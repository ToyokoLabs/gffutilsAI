# GFF Analysis Tools - AI Agent

A comprehensive bioinformatics AI agent for analyzing GFF (General Feature Format) files using natural language queries. This project extends a basic GFF analysis agent with advanced querying capabilities, statistical analysis, and data export features.

<a target="_blank" href="https://colab.research.google.com/github/ToyokoLabs/gffutilsAI/blob/main/notebooks/gffai.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>


## Overview

This AI agent provides an intuitive interface for bioinformatics researchers to analyze GFF files without writing complex code. Users can ask questions in natural language, and the agent will select and execute the appropriate analysis tools to provide comprehensive answers.

## Features

### 🧬 Coordinate-based Queries
- Find features by genomic coordinates (regions and specific positions)
- Query features overlapping genomic regions with filtering by type and strand
- Identify features containing specific genomic positions

### 🔗 Relationship and Hierarchy Analysis
- Explore gene structure and organization (get all child features like exons, CDS, UTRs)
- Find parent features of any given feature using upward traversal
- Get all features of specific types with efficient iteration

### 📊 Statistical Analysis
- Calculate comprehensive feature statistics (counts, length distributions per feature type)
- Generate per-chromosome feature summaries and analysis
- Analyze length distributions with detailed statistics (min, max, mean, median, std dev, percentiles)
- Create histogram data for feature length distributions

### 🔍 Attribute-based Searches
- Search features by attribute key-value pairs (exact and partial matching)
- Find features containing specific attribute keys
- Support pattern matching and logical operations for attribute queries

### 📍 Positional Analysis
- Identify intergenic regions (gaps between genes) with filtering options
- Calculate feature density in genomic windows across chromosomes
- Analyze strand distribution of features with counts and percentages
- Support clustering analysis and positional insights

### 📤 Export and Reporting
- Export feature data to CSV format with comprehensive filtering
- Generate human-readable summary reports of GFF file contents
- Provide formatted output for downstream analysis

## Installation

### Option 1: Install from PyPI (Recommended)

```bash
pip install gffutilsai
```

After installation, you can use the command:

```bash
gffai --help
```

### Option 2: Development Installation

#### Prerequisites

- Python 3.12+
- Ollama (for running local LLM models)

#### Dependencies

Install the required Python packages:

```bash
pip install gffutils strands requests
```

### Ollama Setup

1. Install Ollama from [https://ollama.ai](https://ollama.ai)
2. Pull a compatible model (e.g., llama3.1):
   ```bash
   ollama pull llama3.1
   ```
3. Ensure Ollama is running on `http://localhost:11434`

## Usage

### Running the Agent

#### Basic Usage (Interactive Mode)

**If installed via pip:**
```bash
# Use default settings (llama3.1 model on local server)
gffai

# Use cloud server with default model (gpt-oss:20b-cloud)
gffai --server cloud

# Use Anthropic Claude model (default: claude-3-5-haiku-latest)
gffai --anthropic

# Use Google Gemini model (default: gemini-2.0-flash-exp)
gffai --gemini

# Use OpenAI model (default: gpt-4o-mini)
gffai --openai

# Specify custom model and server
gffai --model llama3.1 --server local
gffai --model codellama:13b --server local
gffai --model gpt-4 --server cloud
gffai --anthropic --model claude-3-5-sonnet-latest
```

**If running from source:**
```bash
# Use default settings (llama3.1 model on local server)
uv run gffai

# Use cloud server with default model (gpt-oss:20b-cloud)
uv run gffai --server cloud

# Use Anthropic Claude model (default: claude-3-5-haiku-latest)
uv run gffai --anthropic

# Use Google Gemini model (default: gemini-2.0-flash-exp)
uv run gffai --gemini

# Use OpenAI model (default: gpt-4o-mini)
uv run gffai --openai
```

If you are going to use a cloud model you need to export the api key.

For Ollama:

```
export OLLAMA_API_KEY="XXXXXXXXXXXXXX"
```

For Anthropic:

```
export ANTHROPIC_API_KEY="XXXXXXXXXXXXX"
```


#### Single Query Mode

```bash
# Run a single query and exit
gffai --query "What feature types are in my GFF file?"
gffai --model llama3.1 --query "Find all genes on chromosome 1"

# Or if running from source:
uv run gffai --query "What feature types are in my GFF file?"
uv run gffai --model llama3.1 --query "Find all genes on chromosome 1"
```

#### Batch Mode (for Benchmarking)

Process multiple queries from a file, one query per line. This is useful for benchmarking and automated testing.

```bash
# Create a queries file (one query per line)
cat > queries.txt << EOF
What feature types are available in my GFF file?
How many genes are in the genome?
List all chromosomes in the GFF file
EOF

# Run in batch mode
gffai --batch queries.txt --model llama3.1

# Or if running from source:
uv run gffai --batch queries.txt --model llama3.1

# With different providers
gffai --batch queries.txt --anthropic
gffai --batch queries.txt --openai --model gpt-4o
```

**Batch file format:**
- One query per line
- Lines starting with `#` are treated as comments and ignored
- Empty lines are ignored
- Results are printed for each query with a summary at the end

**Example batch file (`example_queries.txt`):**
```
# Basic information queries
What feature types are available in my GFF file?
How many genes are in the genome?
List all chromosomes in the GFF file

# Statistical queries
Calculate feature statistics for this GFF file
What's the length distribution of genes?
```

#### Command Line Options

- `--model, -m`: Model to use (default: llama3.1 for local, gpt-oss:20b-cloud for cloud)
- `--server, -s`: Server type - 'local' or 'cloud' (default: local)
- `--anthropic`: Use Anthropic Claude model (default: claude-3-5-haiku-latest)
- `--gemini`: Use Google Gemini model (default: gemini-2.0-flash-exp)
- `--openai`: Use OpenAI model (default: gpt-4o-mini)
- `--host`: Custom host URL (overrides --server setting)
- `--query, -q`: Run a single query and exit
- `--batch, -b`: Run queries from a file (one query per line) for benchmarking
- `--temperature, -t`: Temperature for responses (0.0-1.0, default: 0.1)
- `--max-tokens`: Maximum tokens for responses (default: 4096)
- `--system-prompt`: Path to system prompt file (default: system_prompt.txt)
- `--env-file`: Path to .env file (default: .env in current directory)
- `--version, -v`: Show version information
- `--debug`: Show detailed debug information including tool calls and parameters

#### Server Options

**Local Server (default):**
- Uses `http://localhost:11434`
- Requires Ollama running locally
- Free and private

**Cloud Server:**
- Uses `https://ollama.com`
- Requires `OLLAMA_API_KEY` environment variable
- May have usage costs
- **Security restriction**: `file_read` tool is disabled for security

**Anthropic Claude:**
- Uses Anthropic's Claude models via API
- Requires `ANTHROPIC_API_KEY` environment variable
- **Security restriction**: `file_read` tool is disabled for security
- Default model: `claude-3-5-haiku-latest`

**Google Gemini:**
- Uses Google's Gemini models via API
- Requires `GEMINI_API_KEY` environment variable
- **Security restriction**: `file_read` tool is disabled for security
- Default model: `gemini-2.0-flash-exp`

**OpenAI:**
- Uses OpenAI's models via API
- Requires `OPENAI_API_KEY` environment variable
- **Security restriction**: `file_read` tool is disabled for security
- Default model: `gpt-4o-mini`

The agent will start in interactive mode where you can ask questions about your GFF files, or use `--query` for single commands.

### Example Queries

Here are some example questions you can ask the agent:

#### Basic Information
- "What feature types are available in my GFF file?"
- "How many genes are in chromosome 1?"
- "What's the length of gene AT1G01010?"

#### Coordinate-based Queries
- "Find all genes in chromosome 1 between positions 1000-5000"
- "What features are at position 2500 on chromosome 2?"
- "Show me all exons on the positive strand in the region chr1:10000-20000"

#### Gene Structure Analysis
- "Get the structure of gene AT1G01010 including all exons and CDS"
- "What are the parent features of exon AT1G01010.1?"
- "Show me all CDS features in the genome"

#### Statistical Analysis
- "Calculate feature statistics for this GFF file"
- "What's the length distribution of genes?"
- "Give me a summary of features on each chromosome"

#### Attribute Searches
- "Find all features with 'kinase' in their Name attribute"
- "Show me features that have a Note attribute"
- "Search for genes with 'hypothetical' in their description"

#### Positional Analysis
- "Identify intergenic regions longer than 1000bp on chromosome 2"
- "Calculate gene density in 10kb windows across chromosome 1"
- "What's the strand distribution of genes?"

#### Data Export
- "Export all exon features to CSV format"
- "Generate a summary report of this GFF file"
- "Save gene information to a CSV file"

#### Sample session


```
$ uv run gffai -s cloud -m glm-4.6:cloud
🤖 GFF Analysis AI Agent
📊 Model: glm-4.6:cloud
🌐 Server: cloud (https://ollama.com)
🌡️  Temperature: 0.1
--------------------------------------------------
📝 System prompt loaded from: system_prompt.txt
🔒 Cloud server: file_read tool disabled for security
💬 Interactive mode - Type your questions about GFF files
   Type 'quit' or 'exit' to stop
--------------------------------------------------

🧬 GFF Query: list files in src/gffutilsAI/data/
------------------------------

Tool #1: list_directory
The directory `src/gffutilsAI/data/` contains the following files:

1. **GCF_040556925.1_ASM4055692v1_genomic.gff**
2. **subset_4percent.gff**

Both appear to be GFF (General Feature Format) files, which are commonly used for genome annotation data. Would you like me to analyze either of these files or perform any specific operations on them?
------------------------------

🧬 GFF Query: how many chromosomes are there in GCF_040556925.1_ASM4055692v1_genomic.gff
------------------------------

Tool #2: get_chromosomes_info
The GFF file `GCF_040556925.1_ASM4055692v1_genomic.gff` contains **1 chromosome**:

- **NZ_CP159789.1**

This appears to be a bacterial genome assembly, as indicated by the "NZ_CP" prefix which is typically used for RefSeq complete bacterial genomes.
------------------------------

🧬 GFF Query: which species it is?
------------------------------

Tool #3: get_organism_info
The GFF file `GCF_040556925.1_ASM4055692v1_genomic.gff` contains the genome of **Mycoplasmoides genitalium**, which is a bacterial species (taxonomy ID: 2097).

This is consistent with the fact that there's only one chromosome (NZ_CP159789.1) in the file, as bacteria typically have a single circular chromosome.
------------------------------

🧬 GFF Query: how many genes there are?
------------------------------

Tool #4: get_feature_statistics
The Mycoplasmoides genitalium genome contains **545 genes**.

Here's a breakdown of the genomic features in this assembly:

- **Genes**: 545
- **CDS (coding sequences)**: 524
- **Pseudogenes**: 21
- **tRNA**: 36
- **rRNA**: 3
- **Other RNA features**: 4 (RNase_P_RNA, SRP_RNA, tmRNA)
- **Exons**: 42
- **Total features**: 1,175

The genome has a total length of 1,706,785 base pairs (~1.7 Mb), which is typical for Mycoplasma species, which have some of the smallest known bacterial genomes.
------------------------------

🧬 GFF Query: how many membrane genes? 
------------------------------

Tool #5: search_genes_by_go_function_attribute

Tool #6: search_genes_by_go_function_attribute

Tool #7: search_genes_by_go_function_attribute

Tool #8: search_genes_by_go_function_attribute

Tool #9: search_features_by_attribute

Tool #10: search_features_by_attribute
Based on my analysis of the Mycoplasmoides genitalium genome, I found **approximately 35-40 membrane-related genes**. Here's the breakdown:

## Membrane Genes Found:

### Transporters (Primary Category):
- **ABC transporters**: ~20 genes (ATP-binding proteins, permeases, substrate-binding proteins)
- **Energy-coupling factor (ECF) transporters**: ~5 genes
- **MFS transporters**: 1 gene
- **MATE family efflux transporter**: 1 gene
- **PTS transporters**: ~3 genes (phosphotransferase system)

### Other Membrane Proteins:
- **Membrane protein insertase YidC**: 1 gene
- **Lysylphosphatidylglycerol synthase**: 1 gene (membrane domain-containing)
- **Replication initiation and membrane attachment protein**: 1 gene

## Key Findings:
- **Total membrane-related genes**: ~35-40 out of 545 total genes (~6-7% of genome)
- **Transporters dominate**: ~30 genes are involved in various transport processes
- **ABC transporters are most abundant**: ~20 genes for importing nutrients and exporting waste
- **ECF transporters**: Specialized for vitamin and micronutrient uptake

This is typical for Mycoplasma species, which have reduced genomes and rely heavily on transport systems to import nutrients from their host environment since they lack many metabolic pathways.
------------------------------

🧬 GFF Query: can you write a CSV file with gene name, position and description? call it membrane_genes.csv
------------------------------

Tool #11: file_write
I've successfully created the CSV file `membrane_genes.csv` with **37 membrane-related genes** from the Mycoplasmoides genitalium genome. The file contains:

**Columns:**
- `gene_name`: Locus tag (e.g., ABS333_RS00080)
- `chromosome`: Chromosome identifier (NZ_CP159789.1)
- `start`: Start position in the genome
- `end`: End position in the genome  
- `description`: Protein description

**Summary of the data:**
- **37 total membrane genes** included
- **ABC transporters**: 20 genes (most abundant category)
- **ECF transporters**: 5 genes
- **PTS transporters**: 3 genes
- **Other membrane proteins**: 9 genes (including YidC, MFS transporter, MATE efflux transporter, etc.)

The file is now ready for use in downstream analysis, visualization, or further bioinformatics work. You can open it in Excel, R, Python, or any spreadsheet program to work with the membrane gene data.
------------------------------

```

## Configuration

### Model Configuration

#### Command Line Configuration (Recommended)

Configure the model and server using command line arguments:

```bash
# Local server with different models
gffai --model llama3.1 --server local
gffai --model codellama:13b --server local
gffai --model mistral:7b --server local

# Cloud server with default model (gpt-oss:20b-cloud)
export OLLAMA_API_KEY="your-api-key-here"
gffai --server cloud

# Cloud server with custom model
export OLLAMA_API_KEY="your-api-key-here"
gffai --model gpt-4 --server cloud

# Anthropic Claude with default model (claude-3-5-haiku-latest)
export ANTHROPIC_API_KEY="your-anthropic-api-key"
gffai --anthropic

# Anthropic Claude with custom model
export ANTHROPIC_API_KEY="your-anthropic-api-key"
gffai --anthropic --model claude-3-5-sonnet-latest

# Google Gemini with default model (gemini-2.0-flash-exp)
export GEMINI_API_KEY="your-gemini-api-key"
gffai --gemini

# Google Gemini with custom model
export GEMINI_API_KEY="your-gemini-api-key"
gffai --gemini --model gemini-1.5-pro

# OpenAI with default model (gpt-4o-mini)
export OPENAI_API_KEY="your-openai-api-key"
gffai --openai

# OpenAI with custom model
export OPENAI_API_KEY="your-openai-api-key"
gffai --openai --model gpt-4o

# Custom settings
gffai --model llama3.1 --temperature 0.3 --max-tokens 2048

# Use custom system prompt
gffai --system-prompt my_custom_prompt.txt

# Enable debug mode to see tool calls and parameters
gffai --debug --query "What features are in my GFF file?"

# If running from source, prefix with 'uv run':
uv run gffai --model llama3.1 --server local
```

#### Environment Variables

You can set environment variables in three ways:

**1. Using a .env file (Recommended):**

Create a `.env` file in your working directory:
```bash
# .env file
OLLAMA_API_KEY=your_ollama_api_key_here
ANTHROPIC_API_KEY=your_anthropic_api_key_here
GEMINI_API_KEY=your_gemini_api_key_here
OPENAI_API_KEY=your_openai_api_key_here
```

**2. Using a custom .env file:**
```bash
gffai --env-file path/to/your.env --server cloud
```

**3. Using export commands:**
```bash
export OLLAMA_API_KEY="your-ollama-api-key"
export ANTHROPIC_API_KEY="your-anthropic-api-key"
export GEMINI_API_KEY="your-gemini-api-key"
export OPENAI_API_KEY="your-openai-api-key"
```

The application will automatically load variables from a `.env` file in the current directory if it exists. You can also specify a custom .env file path using the `--env-file` parameter.

#### Available Models

**Local Models** (require `ollama pull <model>`):
- `llama3.1` - General purpose, good balance
- `codellama:13b` - Code-focused, good for technical queries
- `mistral:7b` - Faster, lighter model
- `llama2:70b` - Larger, more capable (requires more resources)

**Cloud Models** (via ollama.com):
- `gpt-oss:20b-cloud` - Default cloud model, good balance of capability and speed
- `gpt-4` - Most capable, requires API key
- `gpt-3.5-turbo` - Fast and capable
- Various other models available through the service

**Anthropic Claude Models**:
- `claude-3-5-haiku-latest` - Default Anthropic model, fast and efficient
- `claude-3-5-sonnet-latest` - More capable, balanced performance
- `claude-3-opus-latest` - Most capable Claude model
- Various other Claude models available

**Google Gemini Models**:
- `gemini-2.0-flash-exp` - Default Gemini model, fast and efficient
- `gemini-1.5-pro` - More capable, balanced performance
- `gemini-1.5-flash` - Fast and lightweight
- Various other Gemini models available

**OpenAI Models**:
- `gpt-4o-mini` - Default OpenAI model, fast and cost-effective
- `gpt-4o` - Most capable OpenAI model, balanced performance
- `gpt-4-turbo` - Fast and capable
- `gpt-3.5-turbo` - Lightweight and fast
- Various other OpenAI models available

### Database Management

The agent automatically creates and manages GFF databases:
- First query creates a database file named after your GFF file (e.g., `file.gff` → `file.db`)
- Subsequent queries reuse the existing database for faster performance
- Database files are created in the same directory as the GFF file
- Multiple GFF files can have their own separate database files

## Project Structure

```
├── main.py              # Main application with CLI interface and agent setup
├── gff_tools.py         # All GFF analysis tool functions
├── system_prompt.txt    # Editable system prompt for the AI agent
├── README.md            # This documentation
└── .kiro/specs/         # Development specifications (optional)
```

## Available Tools

The agent has access to 18+ specialized tools for GFF analysis (defined in `gff_tools.py`):

### File Operations
- `file_read` - Read and display file contents (local server only)
- `file_write` - Write content to files
- `list_directory` - List directory contents

### GFF Analysis Tools
- `get_gff_feature_types` - Get all available feature types
- `get_gene_length` - Get length of specific genes
- `get_multiple_gene_length` - Get lengths of multiple genes
- `get_gene_attributes` - Get gene attributes (ID, Name, Note, etc.)
- `get_features_in_region` - Find features in genomic regions
- `get_features_at_position` - Find features at specific positions
- `get_gene_structure` - Get gene structure with child features
- `get_feature_parents` - Find parent features
- `get_features_by_type` - Get all features of a specific type
- `get_feature_statistics` - Calculate comprehensive statistics
- `get_chromosome_summary` - Per-chromosome analysis
- `get_length_distribution` - Length distribution analysis
- `search_features_by_attribute` - Search by attributes
- `get_features_with_attribute` - Find features with specific attributes
- `get_intergenic_regions` - Identify gaps between genes
- `get_feature_density` - Calculate feature density in windows
- `get_strand_distribution` - Analyze strand distribution
- `export_features_to_csv` - Export data to CSV
- `get_feature_summary_report` - Generate summary reports

## File Formats

### Input
- **GFF3 files** - Standard GFF3 format with proper feature hierarchies
- **GFF2 files** - Older GFF format (limited support)

### Output
- **CSV** - Tabular data export with all feature information
- **JSON** - Structured data format (via agent responses)
- **Text Reports** - Human-readable summaries and statistics

## Performance Considerations

- **Database Creation**: First analysis of a GFF file creates a database, which may take time for large files
- **Memory Usage**: Large GFF files may require significant memory for analysis
- **Result Limiting**: Tools support limiting results for very large datasets
- **Database Reuse**: Subsequent queries on the same GFF file are much faster

### Help and Examples

Get help with command line options:

```bash
gffai --help
```

Example commands:

```bash
# Interactive mode with local llama3.1
gffai

# Interactive mode with cloud default model (gpt-oss:20b-cloud)
export OLLAMA_API_KEY="your-key"
gffai --server cloud

# Interactive mode with cloud GPT-4
export OLLAMA_API_KEY="your-key"
gffai --model gpt-4 --server cloud

# Interactive mode with Anthropic Claude
export ANTHROPIC_API_KEY="your-anthropic-key"
gffai --anthropic

# Interactive mode with Google Gemini
export GEMINI_API_KEY="your-gemini-key"
gffai --gemini

# Interactive mode with OpenAI
export OPENAI_API_KEY="your-openai-key"
gffai --openai

# Single query mode
gffai --query "What chromosomes are in my GFF file?" --model llama3.1

# Batch mode for benchmarking
gffai --batch queries.txt --model llama3.1

# Custom temperature for more creative responses
gffai --temperature 0.5 --model codellama:13b

# Custom host
gffai --host "http://my-server:8080" --model custom-model

# If running from source, prefix with 'uv run':
uv run gffai --help
uv run gffai --server cloud
```

## Troubleshooting

### Common Issues

1. **"File not found" errors**
   - Ensure GFF file path is correct
   - Check file permissions

2. **Database creation fails**
   - Verify GFF file format is valid
   - Check available disk space
   - Ensure write permissions in current directory

3. **Ollama connection errors**
   - Verify Ollama is running: `ollama list`
   - Check if model is available: `ollama pull llama3.1`
   - Ensure correct host URL in configuration

4. **Memory issues with large files**
   - Use result limiting parameters where available
   - Consider analyzing subsets of data
   - Increase system memory if possible

### Debug Mode

Use the `--debug` flag to see detailed information about tool execution:

```bash
# Enable debug mode for single query
gffai --debug --query "What features are in my GFF file?"

# Enable debug mode for interactive session
gffai --debug

# Debug with specific model
gffai --debug --anthropic --query "Analyze my GFF file"

# If running from source:
uv run gffai --debug --query "What features are in my GFF file?"
```

Debug mode shows:
- **Model Information**: Which model and parameters were used
- **Tool Calls**: Which tools were executed and with what parameters
- **Tool Results**: Preview of tool outputs (truncated for readability)
- **Performance Metrics**: Token usage and execution time (when available)
- **Error Details**: Full stack traces for troubleshooting

## Development

### Adding New Tools

To add new GFF analysis tools:

1. Add your tool function to `gff_tools.py` with the `@tool` decorator
2. Import the new tool in `main.py`
3. Add it to the `tools` list in the Agent initialization
4. Update `system_prompt.txt` to describe the new capability

### Customizing the System Prompt

The AI agent's behavior is controlled by the system prompt in `system_prompt.txt`. You can:

1. **Edit the default prompt**: Modify `system_prompt.txt` directly
2. **Use a custom prompt file**: `gffai --system-prompt my_custom_prompt.txt` (or `uv run gffai --system-prompt my_custom_prompt.txt` from source)
3. **Customize for specific use cases**: Create different prompt files for different analysis workflows

The system prompt defines:
- The agent's personality and communication style
- Available capabilities and how to describe them
- Example queries users can ask
- Guidelines for tool usage and error handling

### Project Architecture

- **main.py**: Entry point with CLI argument parsing, model configuration, and agent setup
- **gff_tools.py**: All tool functions decorated with `@tool` for the AI agent to use
- **Modular Design**: Tools are separated from the main application logic for better maintainability

## Contributing

This project follows a specification-driven development approach. See the `.kiro/specs/gff-analysis-tools/` directory for:
- `requirements.md` - Feature requirements
- `design.md` - Technical design
- `tasks.md` - Implementation tasks

## License

GPL-3.0

## Acknowledgments

- Built with [gffutils](https://github.com/daler/gffutils) for GFF file parsing
- Powered by [Strands](https://github.com/weaviate/strands) AI agent framework
- Uses [Ollama](https://ollama.ai) for local LLM inference

## Support

For issues and questions:

sebastian@toyoko.io
