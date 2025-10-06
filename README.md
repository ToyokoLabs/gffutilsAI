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

```bash
python main.py
```

The agent will start and you can interact with it using natural language queries about your GFF files.

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

You can modify the model settings in `main.py`:

```python
model_id = "llama3.1"  # Change to your preferred model

ollama_model = OllamaModel(
    model_id=model_id,
    host="http://localhost:11434",
    params={
        "max_tokens": 4096,
        "temperature": 0.1,  # Lower for more deterministic responses
        "top_p": 0.9,
        "stream": True,
    },
)
```

### Database Management

The agent automatically creates and manages GFF databases:
- First query creates `annotation.db` from your GFF file
- Subsequent queries reuse the existing database for faster performance
- Database files are created in the current working directory

## Available Tools

The agent has access to 18 specialized tools for GFF analysis:

### File Operations
- `file_read` - Read and display file contents
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

To enable more detailed logging, modify the Ollama parameters:

```python
params={
    "max_tokens": 4096,
    "temperature": 0.1,
    "top_p": 0.9,
    "stream": True,
    "verbose": True,  # Add for debugging
}
```

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