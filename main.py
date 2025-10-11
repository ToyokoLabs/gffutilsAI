import os
import argparse
import requests

from strands import Agent
from strands.models.ollama import OllamaModel

# Import all tool functions from the gff_tools module
from gff_tools import (
    file_read, file_write, list_directory,
    get_gff_feature_types, get_gene_lenght, get_gene_attributes, get_multiple_gene_lenght,
    get_all_attributes, get_protein_product_from_gene,
    get_features_in_region, get_features_at_position, get_gene_structure, 
    get_feature_parents, get_features_by_type,
    get_feature_statistics, get_chromosome_summary, get_length_distribution,
    search_features_by_attribute, get_features_with_attribute,
    get_intergenic_regions, get_feature_density, get_strand_distribution,
    export_features_to_csv, get_feature_summary_report
)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="GFF Analysis Tools - AI Agent for bioinformatics analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --model llama3.1 --server local
  python main.py --model gpt-4 --server cloud
  python main.py --model codellama:13b --server local --query "What features are in my GFF file?"
        """
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
    
    args = parser.parse_args()
    
    # Set default model based on server if not specified
    if args.model is None:
        if args.server == "cloud":
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
    print(f"🌐 Server: {args.server} ({host_url})")
    print(f"🌡️  Temperature: {args.temperature}")
    print("-" * 50)
    
    # Load system prompt from file
    with open(args.system_prompt, "r", encoding="utf-8") as f:
        system_prompt = f.read().strip()
    print(f"📝 System prompt loaded from: {args.system_prompt}")



    # Configure the model
    ollama_model = OllamaModel(
        model_id=args.model,
        host=host_url,
        params={
            "max_tokens": args.max_tokens,
            "temperature": args.temperature,
            "top_p": 0.9,
            "stream": True,
        },
    )
    
    # Add authorization header for cloud server if needed
    headers = {}
    if args.server == "cloud" and os.environ.get('OLLAMA_API_KEY'):
        headers['Authorization'] = 'Bearer ' + os.environ.get('OLLAMA_API_KEY')
    
    # Update model configuration with headers if needed
    if headers:
        ollama_model = OllamaModel(
            model_id=args.model,
            host=host_url,
            headers=headers,
            params={
                "max_tokens": args.max_tokens,
                "temperature": args.temperature,
                "top_p": 0.9,
                "stream": True,
            },
        )

    # Create tools list based on server type
    base_tools = [
        file_write, list_directory,
        get_gff_feature_types, get_gene_lenght, get_gene_attributes, get_multiple_gene_lenght,
        get_all_attributes, get_protein_product_from_gene,
        get_features_in_region, get_features_at_position, get_gene_structure, 
        get_feature_parents, get_features_by_type,
        get_feature_statistics, get_chromosome_summary, get_length_distribution,
        search_features_by_attribute, get_features_with_attribute,
        get_intergenic_regions, get_feature_density, get_strand_distribution,
        export_features_to_csv, get_feature_summary_report
    ]
    
    # Add file_read tool only for local server (security restriction for cloud)
    if args.server == "local":
        tools_list = [file_read] + base_tools
        print("🔓 Local server: file_read tool enabled")
    else:
        tools_list = base_tools
        print("🔒 Cloud server: file_read tool disabled for security")

    # Create the agent
    local_agent = Agent(
        system_prompt=system_prompt,
        model=ollama_model,
        tools=tools_list,
    )

    # Handle single query mode or interactive mode
    if args.query:
        print(f"🔍 Query: {args.query}")
        print("-" * 50)
        try:
            result = local_agent(args.query)
            print(result)
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    else:
        print("💬 Interactive mode - Type your questions about GFF files")
        print("   Type 'quit' or 'exit' to stop")
        print("-" * 50)
        
        while True:
            try:
                user_input = input("\n🧬 GFF Query: ").strip()
                
                if user_input.lower() in ['quit', 'exit', 'q']:
                    print("👋 Goodbye!")
                    break
                
                if not user_input:
                    continue
                
                print("-" * 30)
                result = local_agent(user_input)
                print(result)
                print("-" * 30)
                
            except KeyboardInterrupt:
                print("\n👋 Goodbye!")
                break
            except Exception as e:
                print(f"❌ Error: {str(e)}")
                print("Please try again or type 'quit' to exit.")


if __name__ == "__main__":
    main()