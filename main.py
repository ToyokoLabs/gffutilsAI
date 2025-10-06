

import os

import requests

from strands import Agent, tool
from strands.models.ollama import OllamaModel
import gffutils




@tool
def file_read(file_path: str) -> str:
    """Read a file and return its content.

    Args:
        file_path (str): Path to the file to read

    Returns:
        str: Content of the file

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        with open(file_path, "r") as file:
            return file.read()
    except FileNotFoundError:
        return f"Error: File '{file_path}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


#gffutils.interface.FeatureDB

@tool
def get_gff_feature_types(gffpath: str) -> list:
    """Given the path of a gff file, generate a database and then get the list of all available features types.
    Sample features are:
    'CDS', 'chromosome', 'exon', 'five_prime_UTR',
    'gene', 'mRNA', 'mRNA_TE_gene', 'miRNA', 'ncRNA', 'protein'

    Args:
        gffpath (str): Path to the file to read

    Returns:
        list: All available feature types

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # check if file is there
        if os.path.exists('annotation.db'):
            # is here
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        return list(db.featuretypes())
    except FileNotFoundError:
        return f"Error: File '{file_path}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


@tool
def get_gene_lenght(gffpath: str, gene_id: str) -> list:
    """From a gff file and a gene id, returns the lenght of the gene.

    Args:
        gffpath (str): Path to the file to read
        gene_id (str): The gene name

    Returns:
        list: The lenght of the gene

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # check if file is there
        if os.path.exists('annotation.db'):
            # is here
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        g = db[gene_id]
        return abs(g.start-g.end)
    except FileNotFoundError:
        return f"Error: File '{file_path}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


@tool
def get_gene_attributes(gffpath: str, gene_id: str) -> dict:
    """From a gff file and a gene id, returns gene attributes, these are the gene attributes: ID, Note, Name.

    Args:
        gffpath (str): Path to the file to read
        gene_id (str): The gene name

    Returns:
        dictionary: A dictionary with the gene attributes. For example: 
        {'ID': ['AT1G01183'], 'Note': ['miRNA'], 'Name': ['AT1G01183']}

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # check if file is there
        if os.path.exists('annotation.db'):
            # is here
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        g = db[gene_id]
        return dict(g.attributes.items())
    except FileNotFoundError:
        return f"Error: File '{file_path}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


#print(f"Found gene: {found_gene.id}")
#    print(f"Chromosome: {found_gene.chrom}")
#    print(f"Start: {found_gene.start}")
#    print(f"End: {found_gene.end}")
#    print(f"Strand: {found_gene.strand}")

@tool
def get_multiple_gene_lenght(gffpath: str, gene_ids: list) -> list:
    """From a gff file and a list of gene ids, returns a list with the lenght of all genes.

    Args:
        gffpath (str): Path to the file to read
        gene_ids (list): The gene name

    Returns:
        list: The lenght of the gene

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # check if file is there
        if os.path.exists('annotation.db'):
            # is here
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        out = []
        for gid in gene_ids:
            g = db[gid]
            out.append(abs(g.start-g.end))
        return out
    except FileNotFoundError:
        return f"Error: File '{file_path}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"

@tool
def file_write(file_path: str, content: str) -> str:
    """Write content to a file.

    Args:
        file_path (str): The path to the file
        content (str): The content to write to the file

    Returns:
        str: A message indicating success or failure
    """
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)

        with open(file_path, "w") as file:
            file.write(content)
        return f"File '{file_path}' written successfully."
    except Exception as e:
        return f"Error writing to file: {str(e)}"


@tool
def list_directory(directory_path: str = ".") -> str:
    """List files and directories in the specified path.

    Args:
        directory_path (str): Path to the directory to list

    Returns:
        str: A formatted string listing all files and directories
    """
    try:
        items = os.listdir(directory_path)
        files = []
        directories = []

        for item in items:
            full_path = os.path.join(directory_path, item)
            if os.path.isdir(full_path):
                directories.append(f"Folder: {item}/")
            else:
                files.append(f"File: {item}")

        result = f"Contents of {os.path.abspath(directory_path)}:\n"
        result += (
            "\nDirectories:\n" + "\n".join(sorted(directories))
            if directories
            else "\nNo directories found."
        )
        result += (
            "\n\nFiles:\n" + "\n".join(sorted(files)) if files else "\nNo files found."
        )

        return result
    except Exception as e:
        return f"Error listing directory: {str(e)}"

# Define a comprehensive system prompt for our agent
system_prompt = """
You are a helpful personal assistant with bioinformatics capabilities and can perfom local file actions and simple tasks for the user.

Your key capabilities:
1. Read, understand, and summarize files.
2. Create and write to files.
3. List directory contents and provide information on the files.
4. Summarize text content
5. Get information out of gff files
6. Get the lenght of a gene

When using tools:
- Always verify file paths before operations
- Be careful with system commands
- Provide clear explanations of what you're doing
- If a task cannot be completed, explain why and suggest alternatives

Always be helpful, concise, and focus on addressing the user's needs efficiently.
"""

model_id = (
    "llama3.1"  # You can change this to any model you have pulled with Ollama.
)


ollama_model = OllamaModel(
    model_id=model_id,
    host="http://localhost:11434",
    params={
        "max_tokens": 4096,  # Adjust based on your model's capabilities
        "temperature": 0.1,  # Lower for more deterministic responses, higher for more creative
        "top_p": 0.9,  # Nucleus sampling parameter
        "stream": True,  # Enable streaming responses
    },
)

# Create the agent
local_agent = Agent(
    system_prompt=system_prompt,
    model=ollama_model,
    tools=[file_read, file_write, list_directory, get_gff_feature_types, get_gene_lenght, 
           get_multiple_gene_lenght, get_gene_attributes],
)



def main():
    # Define a comprehensive system prompt for our agent
    system_prompt = """
You are a helpful personal assistant with bioinformatics capabilities and can perfom local file actions and simple tasks for the user.

Your key capabilities:
1. Read, understand, and summarize files.
2. Create and write to files.
3. List directory contents and provide information on the files.
4. Summarize text content
5. Get information out of gff files
6. Get the lenght of a gene

When using tools:
- Always verify file paths before operations
- Be careful with system commands
- Provide clear explanations of what you're doing
- If a task cannot be completed, explain why and suggest alternatives

Always be helpful, concise, and focus on addressing the user's needs efficiently.
"""

    model_id = (
    "llama3.1"  # You can change this to any model you have pulled with Ollama.
)


    ollama_model = OllamaModel(
    model_id=model_id,
    host="http://localhost:11434",
    params={
        "max_tokens": 4096,  # Adjust based on your model's capabilities
        "temperature": 0.1,  # Lower for more deterministic responses, higher for more creative
        "top_p": 0.9,  # Nucleus sampling parameter
        "stream": True,  # Enable streaming responses
    },
)

    # Create the agent
    local_agent = Agent(
    system_prompt=system_prompt,
    model=ollama_model,
    tools=[file_read, file_write, list_directory, get_gff_feature_types, get_gene_lenght, 
           get_multiple_gene_lenght, get_gene_attributes],
)


    print("Hello from gffutilsai!")
    r=local_agent(
    "want to know which features there are in a gff file called ./subset_4percent.gff. Show me the list"
)


if __name__ == "__main__":
    main()
