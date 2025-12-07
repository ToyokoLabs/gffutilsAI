import os
import argparse
import requests
import time
from importlib import resources
from dotenv import load_dotenv
from importlib.metadata import version

from strands import Agent
from strands.models.ollama import OllamaModel
from strands.models.anthropic import AnthropicModel
from strands.models.gemini import GeminiModel
from strands.models.openai import OpenAIModel


# Import all tool functions from the gff_tools module
from .gff_tools import (
    file_read, file_write, list_files, get_country_or_region,
    get_gff_feature_types, get_gene_lenght, get_gene_attributes, get_multiple_gene_lenght,
    get_all_attributes, get_protein_product_from_gene,
    get_features_in_region, get_features_at_position, get_gene_structure, 
    get_feature_parents, get_features_by_type,
    get_feature_statistics, get_chromosome_summary, get_length_distribution,
    search_features_by_attribute, get_features_with_attribute,
    get_intergenic_regions, get_feature_density, get_strand_distribution,
    export_features_to_csv, get_feature_summary_report, get_genes_and_features_from_attribute,
    get_tools_list, get_organism_info, get_chromosomes_info, 
    search_genes_by_go_function_attribute, extract_genes_to_gff
)


__version__ = version("gffutilsai")

# Global variable to store tool call information for debugging
tool_call_log = []


def main():
    global tool_call_log
    
    # Check if --version is being used (don't load env vars for version check)
    import sys
    is_version_check = "--version" in sys.argv or "-v" in sys.argv
    
    # Parse command line arguments first to get env-file option
    env_file_path = None
    if "--env-file" in sys.argv:
        try:
            env_file_index = sys.argv.index("--env-file")
            if env_file_index + 1 < len(sys.argv):
                env_file_path = sys.argv[env_file_index + 1]
        except (ValueError, IndexError):
            pass
    
    # Load environment variables from .env file (but don't show message for version check)
    if env_file_path:
        load_dotenv(env_file_path)
        if not is_version_check:
            if os.path.exists(env_file_path):
                print(f"🔧 Loaded environment variables from: {env_file_path}")
            else:
                print(f"⚠️  Warning: .env file not found: {env_file_path}")
    else:
        # Try to load from default .env file
        if os.path.exists(".env"):
            load_dotenv()
            if not is_version_check:
                print("🔧 Loaded environment variables from: .env")
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="GFF Analysis Tools - AI Agent for bioinformatics analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  gffai --model llama3.1 --server local
  gffai --model gpt-4 --server cloud
  gffai --model codellama:13b --server local --query "What features are in my GFF file?"
  
  From source:
  uv run gffai --model llama3.1 --server local
  
  Environment variables:
  gffai --env-file my.env --server cloud
  
  Provider examples:
  gffai --anthropic --model claude-3-5-sonnet-latest
  gffai --gemini --model gemini-2.0-flash-exp
  gffai --openai --model gpt-4o
  
  Batch mode (for benchmarking):
  gffai --batch queries.txt --model llama3.1
  
  Note: To use cloud models you need to set the API key as an environment variable. 
  You can use a .env file or export the variables directly. See README.md for more information.
        """
    )
    
    parser.add_argument(
        "--version", "-v",
        action="version",
        version=__version__
    )
    
    parser.add_argument(
        "--model", "-m",
        type=str,
        help="Model to use. Default: llama3.1 for local, gpt-oss:20b-cloud for cloud. Examples: llama3.1, codellama:13b, gpt-4, etc."
    )
    
    parser.add_argument(
        "--server", "-s",
        type=str,
        choices=["local", "cloud"],
        default="local",
        help="Server to use: 'local' for localhost:11434 or 'cloud' for ollama.com (default: local)"
    )
    
    parser.add_argument(
        "--host",
        type=str,
        help="Custom host URL (overrides --server setting)"
    )
    
    parser.add_argument(
        "--query", "-q",
        type=str,
        help="Run a single query and exit"
    )
    
    parser.add_argument(
        "--batch", "-b",
        type=str,
        help="Run queries from a file (one query per line) for benchmarking"
    )
    
    parser.add_argument(
        "--temperature", "-t",
        type=float,
        default=0.1,
        help="Temperature for model responses (0.0-1.0, default: 0.1)"
    )
    
    parser.add_argument(
        "--max-tokens",
        type=int,
        default=4096,
        help="Maximum tokens for model responses (default: 4096)"
    )
    
    parser.add_argument(
        "--system-prompt",
        type=str,
        default="system_prompt.txt",
        help="Path to system prompt file (default: system_prompt.txt)"
    )
    
    parser.add_argument(
        "--anthropic",
        action="store_true",
        help="Use Anthropic Claude model (default: claude-3-5-haiku-latest)"
    )
    
    parser.add_argument(
        "--gemini",
        action="store_true",
        help="Use Google Gemini model (default: gemini-2.0-flash-exp)"
    )
    
    parser.add_argument(
        "--openai",
        action="store_true",
        help="Use OpenAI model (default: gpt-4o-mini)"
    )
    
    parser.add_argument(
        "--env-file",
        type=str,
        help="Path to .env file (default: .env in current directory)"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show detailed debug information including tool calls and parameters"
    )
    
    args = parser.parse_args()
    
    # Set default model based on server/provider if not specified
    if args.model is None:
        if args.anthropic:
            args.model = "claude-3-5-haiku-latest" # claude-sonnet-4-5-20250929
        elif args.gemini:
            args.model = "gemini-2.5-flash"
        elif args.openai:
            # old gpt-4.1-mini
            args.model = "gpt-5-nano"
        elif args.server == "cloud":
            args.model = "gpt-oss:20b-cloud"
        else:
            args.model = "llama3.1"
    
    # Determine host URL
    if args.host:
        host_url = args.host
    elif args.server == "cloud":
        host_url = "https://ollama.com"
    else:  # local
        host_url = "http://localhost:11434"
    
    print(f"🤖 GFF Analysis AI Agent")
    print(f"📊 Model: {args.model}")
    if args.gemini:
        print(f"� Provider: Google Gemini")
    elif args.anthropic:
        print(f"🌐 Provider: Anthropic")
    elif args.openai:
        print(f"🌐 Provider: OpenAI")
    else:
        print(f"🌐 Server: {args.server} ({host_url})")
    if args.model.startswith("gpt-5"):
        print(f"🌡️ Temperature: 1")
    else:
        print(f"🌡️ Temperature: {args.temperature}")
    print("-" * 50)
    
    # Load system prompt from file
    try:
        # First try to open the file as specified (for custom prompts or development)
        with open(args.system_prompt, "r", encoding="utf-8") as f:
            system_prompt = f.read().strip()
        print(f"📝 System prompt loaded from: {args.system_prompt}")
    except FileNotFoundError:
        # If default system_prompt.txt not found, try to load from package resources
        if args.system_prompt == "system_prompt.txt":
            try:
                with resources.open_text("gffutilsAI", "system_prompt.txt") as f:
                    system_prompt = f.read().strip()
                print(f"📝 System prompt loaded from package resources")
            except FileNotFoundError:
                print(f"❌ Error: Could not find system prompt file: {args.system_prompt}")
                print("   Make sure the file exists or use --system-prompt to specify a custom file.")
                return
        else:
            print(f"❌ Error: Could not find system prompt file: {args.system_prompt}")
            return



    # Configure the model based on provider
    if args.anthropic:
        # Use Anthropic Claude model
        a_model = AnthropicModel(
            client_args={
                "api_key": os.environ.get('ANTHROPIC_API_KEY', ""),
            },
            max_tokens=args.max_tokens,
            model_id=args.model,
            params={
                "temperature": args.temperature,
            }
        )
        model_to_use = a_model
        print(f"🤖 Using Anthropic Claude model: {args.model}")
    elif args.gemini:
        # Use Google Gemini model
        gemini_model = GeminiModel(
            client_args={
                "api_key": os.environ.get('GEMINI_API_KEY', ""),
            },
            model_id=args.model,
            params={
                "temperature": args.temperature,
                "max_output_tokens": args.max_tokens,
                "top_p": 0.9,
                "top_k": 40,
            }
        )
        model_to_use = gemini_model
        print(f"🤖 Using Google Gemini model: {args.model}")
    elif args.openai:
        # Use OpenAI model
        if args.model.startswith("gpt-4"):
            openai_model = OpenAIModel(
                client_args={
                    "api_key": os.environ.get('OPENAI_API_KEY', ""),
                },
                model_id=args.model,
                params={
                    "temperature": args.temperature,
                    "max_tokens": args.max_tokens,
                }
            )
        else:
            # for 5 and better
            openai_model = OpenAIModel(
                client_args={
                    "api_key": os.environ.get('OPENAI_API_KEY', ""),
                },
                model_id=args.model,
                params={
                    "temperature": 1,
                    "max_completion_tokens": args.max_tokens,
                }
            )
        model_to_use = openai_model
        print(f"🤖 Using OpenAI model: {args.model}")
    else:
        # Use Ollama model
        ollama_model = OllamaModel(
            model_id=args.model,
            host=host_url,
            max_tokens=args.max_tokens,
            temperature=args.temperature,
            top_p=0.9,
        )
        model_to_use = ollama_model
        
        # Note: For cloud server authentication, you may need to set OLLAMA_API_KEY
        # as an environment variable or configure it differently based on the SDK version

    # Create tools list based on server type
    base_tools = [
        file_write, list_files,
        get_gff_feature_types, get_gene_lenght, get_gene_attributes, get_multiple_gene_lenght,
        get_all_attributes, get_protein_product_from_gene,
        get_features_in_region, get_features_at_position, get_gene_structure, 
        get_feature_parents, get_features_by_type,
        get_feature_statistics, get_chromosome_summary, get_length_distribution,
        search_features_by_attribute, get_features_with_attribute,
        get_intergenic_regions, get_feature_density, get_strand_distribution,
        export_features_to_csv, get_feature_summary_report, get_tools_list, 
        get_genes_and_features_from_attribute, get_organism_info, get_chromosomes_info,
        search_genes_by_go_function_attribute, extract_genes_to_gff
    ]
    
    # Add file_read tool only for local server (security restriction for cloud/anthropic/gemini/openai)
    if args.server == "local" and not args.anthropic and not args.gemini and not args.openai:
        all_tools = [file_read] + base_tools
        print("🔓 Local server: file_read tool enabled")
    else:
        all_tools = base_tools
        if args.anthropic:
            print("🔒 Anthropic: file_read tool disabled for security")
        elif args.gemini:
            print("🔒 Gemini: file_read tool disabled for security")
        elif args.openai:
            print("🔒 OpenAI: file_read tool disabled for security")
        else:
            print("🔒 Cloud server: file_read tool disabled for security")

    local_agent = Agent(
        system_prompt=system_prompt,
        model=model_to_use,
        tools=base_tools,
    )
    
    # Handle single query mode, batch mode, or interactive mode
    if args.query:
        print(f"🔍 Query: {args.query}")
        print("-" * 50)
        try:
            # Execute the query
            result = local_agent(args.query)
            if args.debug:
                print("**DEBUG****************************")
                metrics = result.metrics
                for tool_name, tool_metric in metrics.tool_metrics.items():
                    tool_data = tool_metric.tool     
                    print(f"Tool Name: {tool_data['name']}")
                    print(f"Parameters: {tool_data['input']}")
                    print(f"Execution Time: {tool_metric.total_time:.4f}s")
                    print("-------------------------")

                print("***********************************")
                
        except Exception as e:
            print(f"❌ Error: {str(e)}")
            if args.debug:
                import traceback
                print("\n🔧 DEBUG - Full Error Traceback:")
                print("-" * 40)
                traceback.print_exc()
                print("-" * 40)
    elif args.batch:
        # Batch mode - process queries from file
        print(f"📋 Batch mode - Processing queries from: {args.batch}")
        print("-" * 50)
        
        try:
            with open(args.batch, 'r', encoding='utf-8') as f:
                queries = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
            
            if not queries:
                print("⚠️  Warning: No queries found in file")
                return
            
            print(f"📊 Found {len(queries)} queries to process\n")
            
            results = []
            for i, query in enumerate(queries, 1):
                print(f"\n{'='*60}")
                print(f"Query {i}/{len(queries)}: {query}")
                print('='*60)
                
                try:
                    result = local_agent(query)
                    results.append({
                        'query': query,
                        'result': result,
                        'status': 'success'
                    })
                except Exception as e:
                    error_msg = str(e)
                    print(f"❌ Error: {error_msg}")
                    results.append({
                        'query': query,
                        'error': error_msg,
                        'status': 'error'
                    })
                    if args.debug:
                        import traceback
                        traceback.print_exc()
                # wait time to avoid rate limit error
                time.sleep(2)
            
            # Print summary
            print(f"\n{'='*60}")
            print("📊 Batch Processing Summary")
            print('='*60)
            successful = sum(1 for r in results if r['status'] == 'success')
            failed = sum(1 for r in results if r['status'] == 'error')
            print(f"✅ Successful: {successful}/{len(queries)}")
            print(f"❌ Failed: {failed}/{len(queries)}")
            print('='*60)
            print(f"Version: {__version__}")
            
        except FileNotFoundError:
            print(f"❌ Error: File not found: {args.batch}")
        except Exception as e:
            print(f"❌ Error reading batch file: {str(e)}")
            if args.debug:
                import traceback
                traceback.print_exc()
    else:
        print("💬 Interactive mode - Type your questions about GFF files")
        print("   Type 'quit' or 'exit' to stop")
        print("-" * 50)
        
        while True:
            try:
                user_input = input("\n🧬 GFF Query: ").strip()
                
                if user_input.lower() in ['quit', 'exit', 'q']:
                    print("👋 Goodbye!")
                    print("\n" + "=" * 80)
                    print("To cite this program, use this:")
                    print("gffutilsAI: An AI-Agent for Interactive Genomic Feature Exploration in GFF files")
                    print("Sebastian Bassi, Tristan Yang, Virginia Gonzalez, bioRxiv 2025.12.02.690645;")
                    print("doi: https://doi.org/10.64898/2025.12.02.690645")
                    print("=" * 80)
                    break
                
                if not user_input:
                    continue
                
                print("-" * 30)
                
                # Execute the query
                result = local_agent(user_input)
                
                #print(result)
                
                print("\n" + "-" * 30)
                
            except KeyboardInterrupt:
                print("\n👋 Goodbye!")
                print("\n" + "=" * 80)
                print("To cite this program, use this:")
                print("gffutilsAI: An AI-Agent for Interactive Genomic Feature Exploration in GFF files")
                print("Sebastian Bassi, Tristan Yang, Virginia Gonzalez, bioRxiv 2025.12.02.690645;")
                print("doi: https://doi.org/10.64898/2025.12.02.690645")
                print("=" * 80)
                break
            except Exception as e:
                print(f"❌ Error: {str(e)}")
                print("Please try again or type 'quit' to exit.")


if __name__ == "__main__":
    main()