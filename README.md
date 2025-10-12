# GFF Analysis Tools - AI Agent

A comprehensive bioinformatics AI agent for analyzing GFF (General Feature Format) files using natural language queries. This project extends a basic GFF analysis agent with advanced querying capabilities, statistical analysis, and data export features.

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

### Prerequisites

- Python 3.8+
- Ollama (for running local LLM models)

### Dependencies

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

```bash
# Use default settings (llama3.1 model on local server)
python main.py

# Use cloud server with default model (gpt-oss:20b-cloud)
python main.py --server cloud

# Use Anthropic Claude model (default: claude-3-5-haiku-latest)
python main.py --anthropic

# Specify custom model and server
python main.py --model llama3.1 --server local
python main.py --model codellama:13b --server local
python main.py --model gpt-4 --server cloud
python main.py --anthropic --model claude-3-5-sonnet-latest
```

#### Single Query Mode

```bash
# Run a single query and exit
python main.py --query "What feature types are in my GFF file?"
python main.py --model llama3.1 --query "Find all genes on chromosome 1"
```

#### Command Line Options

- `--model, -m`: Model to use (default: llama3.1 for local, gpt-oss:20b-cloud for cloud)
- `--server, -s`: Server type - 'local' or 'cloud' (default: local)
- `--anthropic`: Use Anthropic Claude model (default: claude-3-5-haiku-latest)
- `--host`: Custom host URL (overrides --server setting)
- `--query, -q`: Run a single query and exit
- `--temperature, -t`: Temperature for responses (0.0-1.0, default: 0.1)
- `--max-tokens`: Maximum tokens for responses (default: 4096)
- `--system-prompt`: Path to system prompt file (default: system_prompt.txt)
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

## Configuration

### Model Configuration

#### Command Line Configuration (Recommended)

Configure the model and server using command line arguments:

```bash
# Local server with different models
python main.py --model llama3.1 --server local
python main.py --model codellama:13b --server local
python main.py --model mistral:7b --server local

# Cloud server with default model (gpt-oss:20b-cloud)
export OLLAMA_API_KEY="your-api-key-here"
python main.py --server cloud

# Cloud server with custom model
export OLLAMA_API_KEY="your-api-key-here"
python main.py --model gpt-4 --server cloud

# Anthropic Claude with default model (claude-3-5-haiku-latest)
export ANTHROPIC_API_KEY="your-anthropic-api-key"
python main.py --anthropic

# Anthropic Claude with custom model
export ANTHROPIC_API_KEY="your-anthropic-api-key"
python main.py --anthropic --model claude-3-5-sonnet-latest

# Custom settings
python main.py --model llama3.1 --temperature 0.3 --max-tokens 2048

# Use custom system prompt
python main.py --system-prompt my_custom_prompt.txt

# Enable debug mode to see tool calls and parameters
python main.py --debug --query "What features are in my GFF file?"
```

#### Environment Variables

For cloud server usage, set your API key:

```bash
export OLLAMA_API_KEY="your-ollama-api-key"
```

For Anthropic Claude usage, set your API key:

```bash
export ANTHROPIC_API_KEY="your-anthropic-api-key"
```

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
python main.py --help
```

Example commands:

```bash
# Interactive mode with local llama3.1
python main.py

# Interactive mode with cloud default model (gpt-oss:20b-cloud)
export OLLAMA_API_KEY="your-key"
python main.py --server cloud

# Interactive mode with cloud GPT-4
export OLLAMA_API_KEY="your-key"
python main.py --model gpt-4 --server cloud

# Interactive mode with Anthropic Claude
export ANTHROPIC_API_KEY="your-anthropic-key"
python main.py --anthropic

# Single query mode
python main.py --query "What chromosomes are in my GFF file?" --model llama3.1

# Custom temperature for more creative responses
python main.py --temperature 0.5 --model codellama:13b

# Custom host
python main.py --host "http://my-server:8080" --model custom-model
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
python main.py --debug --query "What features are in my GFF file?"

# Enable debug mode for interactive session
python main.py --debug

# Debug with specific model
python main.py --debug --anthropic --query "Analyze my GFF file"
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
2. **Use a custom prompt file**: `python main.py --system-prompt my_custom_prompt.txt`
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

[Add your license information here]

## Acknowledgments

- Built with [gffutils](https://github.com/daler/gffutils) for GFF file parsing
- Powered by [Strands](https://github.com/weaviate/strands) AI agent framework
- Uses [Ollama](https://ollama.ai) for local LLM inference

## Support

For issues and questions:

sebastian@toyoko.io