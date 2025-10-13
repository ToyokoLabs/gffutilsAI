"""
GFF Analysis Tools - Tool Functions

This module contains all the tool functions for GFF file analysis.
These functions are decorated with @tool to be used by the AI agent.
"""

import os
import csv
from strands import tool
import gffutils
from Bio import Entrez


def get_db_filename(gffpath: str) -> str:
    """Generate database filename based on GFF file path.
    
    Args:
        gffpath (str): Path to the GFF file
        
    Returns:
        str: Database filename (e.g., 'file.gff' -> 'file.db')
    """
    base_name = os.path.splitext(gffpath)[0]
    return f"{base_name}.db"


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


@tool
def get_organism_info(accession: str) -> dict:
    """Given an accession id get organism info like specie name

    Args:
        accession (str): Accession id like GCF_036512215.1

    Returns:
        dict: organism and specie information

    """

    def search_refseq_assembly(refseq_accession):
        # Step 1: Search assembly database for this accession
        handle = Entrez.esearch(db="assembly", term=refseq_accession)
        search_results = Entrez.read(handle)
        handle.close()
        
        uid_list = search_results['IdList']
        return uid_list


    def get_assembly_summary(uid):
        # Step 2: Fetch summary info by UID
        handle = Entrez.esummary(db="assembly", id=uid)
        summary = Entrez.read(handle, validate=False)
        handle.close()
        return summary


    Entrez.email = "example@example.com"

    uids = search_refseq_assembly(accession)

    if uids:
        uid = uids[0]  # take the first match
        summary = get_assembly_summary(uid)
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
    else:
        print("No assembly found for accession", refseq_accession)

    organism = docsum.get('Organism', 'N/A')
    taxonomy_id = docsum.get('Taxid', 'N/A')
    species = docsum.get('SpeciesName', 'N/A')  # Often available
    return {"organism": organism, "taxonomy_id": taxonomy_id, "species": species}


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
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        return list(db.featuretypes())
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
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
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        g = db[gene_id]
        return abs(g.start-g.end)
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
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
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        g = db[gene_id]
        return dict(g.attributes.items())
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


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
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        out = []
        for gid in gene_ids:
            g = db[gid]
            out.append(abs(g.start-g.end))
        return out
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


@tool
def get_all_attributes(gffpath: str) -> set:
    """Given the path of a gff file, generate a database and then get the list of all available attributes.
    Sample attributes are:
    'Dbxref', 'ID', 'Is_circular', 'Name', 'Note', 'Ontology_term', 'Parent', 'anticodon', 'collection-date', 'country',
    'end_range', 'gb-synonym', 'gbkey', 'gene', 'gene_biotype', 'gene_synonym',  'genome', 'go_component', 'go_function', 'go_process',
    'inference', 'isolation-source',  'locus_tag',  'mol_type', 'nat-host', 'partial',
     'product', 'protein_id', 'pseudo', 'start_range', 'strain', 'transl_table', 'type-material'

    Args:
        gffpath (str): Path to the file to read

    Returns:
        set: All available attributes

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        attribute_types = set()
        for feature in db.all_features():
            attribute_types.update(feature.attributes.keys())
        return attribute_types
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


@tool
def get_genes_and_features_from_attribute(gffpath: str, attr: str) -> dict:
    """Get all genes that has an attribute.

    Args:
        gffpath (str): Path to the file to read
        attr (str): The attribute content (for example "regulation of sporulation")

    Returns:
        Dict (dict): A dictionary with all genes and features where the attribute is present

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        genes = set()
        features = []
        for feature in db.all_features():    
            if any(attr in val for vals in feature.attributes.values() for val in vals):
                #print(f"{feature.id}\t{feature.featuretype}\t{feature.attributes}")
                features.append(feature.id)
                gene = get_feature_parents(gffpath, feature.id)[0]['id']
                genes.add(gene)

        return {"genes": genes, "features": features}
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"


@tool
def get_protein_product_from_gene(gffpath: str, gene: str) -> list:
    """Get all protein products given a gene name or id.

    Args:
        gffpath (str): Path to the file to read
        gene (str): The gene name or id

    Returns:
        list: List of protein products

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get all children of the gene and look for CDS features
        gene_feature = db[gene]
        products = []
        
        for child in db.children(gene_feature, featuretype='CDS'):
            if 'product' in child.attributes:
                products.extend(child.attributes['product'])
        
        return products
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error reading file: {str(e)}"

@tool
def get_features_in_region(gffpath: str, chrom: str, start: int, end: int, feature_type: str = None, strand: str = None) -> list:
    """Find all features overlapping a genomic region.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str): Chromosome name
        start (int): Start coordinate
        end (int): End coordinate
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon')
        strand (str, optional): Filter by strand ('+', '-', or '.')

    Returns:
        list: List of dictionaries containing feature information

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Query features in the region
        features = []
        for feature in db.region(seqid=chrom, start=start, end=end):
            # Apply feature type filter if specified
            if feature_type and feature.featuretype != feature_type:
                continue
            
            # Apply strand filter if specified
            if strand and feature.strand != strand:
                continue
            
            # Convert feature to dictionary
            feature_dict = {
                'id': feature.id,
                'chrom': feature.chrom,
                'start': feature.start,
                'end': feature.end,
                'strand': feature.strand,
                'feature_type': feature.featuretype,
                'attributes': dict(feature.attributes.items()),
                'length': abs(feature.end - feature.start)
            }
            features.append(feature_dict)
        
        return features
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error querying features: {str(e)}"


@tool
def get_features_at_position(gffpath: str, chrom: str, position: int, feature_type: str = None) -> list:
    """Find features that contain a specific genomic position.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str): Chromosome name
        position (int): Genomic position to query
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon')

    Returns:
        list: List of dictionaries containing feature information

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Query features at the position (using a single position as both start and end)
        features = []
        for feature in db.region(seqid=chrom, start=position, end=position):
            # Apply feature type filter if specified
            if feature_type and feature.featuretype != feature_type:
                continue
            
            # Check if the position is actually within the feature bounds
            if feature.start <= position <= feature.end:
                # Convert feature to dictionary
                feature_dict = {
                    'id': feature.id,
                    'chrom': feature.chrom,
                    'start': feature.start,
                    'end': feature.end,
                    'strand': feature.strand,
                    'feature_type': feature.featuretype,
                    'attributes': dict(feature.attributes.items()),
                    'length': abs(feature.end - feature.start)
                }
                features.append(feature_dict)
        
        return features
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error querying features: {str(e)}"


@tool
def get_gene_structure(gffpath: str, gene_id: str) -> dict:
    """Get all child features of a gene (exons, CDS, UTRs) organized by feature type.

    Args:
        gffpath (str): Path to the GFF file
        gene_id (str): The gene ID to query

    Returns:
        dict: Dictionary with gene structure organized by feature type
              Format: {
                  'gene_info': {...},
                  'children': {
                      'exon': [...],
                      'CDS': [...],
                      'UTR': [...]
                  }
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get the gene feature
        try:
            gene = db[gene_id]
        except gffutils.exceptions.FeatureNotFoundError:
            return f"Error: Gene '{gene_id}' not found in database."
        
        # Get gene information
        gene_info = {
            'id': gene.id,
            'chrom': gene.chrom,
            'start': gene.start,
            'end': gene.end,
            'strand': gene.strand,
            'feature_type': gene.featuretype,
            'attributes': dict(gene.attributes.items()),
            'length': abs(gene.end - gene.start)
        }
        
        # Get all child features and organize by type
        children_by_type = {}
        for child in db.children(gene):
            feature_type = child.featuretype
            
            if feature_type not in children_by_type:
                children_by_type[feature_type] = []
            
            child_dict = {
                'id': child.id,
                'chrom': child.chrom,
                'start': child.start,
                'end': child.end,
                'strand': child.strand,
                'feature_type': child.featuretype,
                'attributes': dict(child.attributes.items()),
                'length': abs(child.end - child.start)
            }
            children_by_type[feature_type].append(child_dict)
        
        # Sort children within each type by start position
        for feature_type in children_by_type:
            children_by_type[feature_type].sort(key=lambda x: x['start'])
        
        return {
            'gene_info': gene_info,
            'children': children_by_type
        }
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error getting gene structure: {str(e)}"


@tool
def get_feature_parents(gffpath: str, feature_id: str) -> list:
    """Find parent features of any given feature using upward traversal.

    Args:
        gffpath (str): Path to the GFF file
        feature_id (str): The feature ID to find parents for

    Returns:
        list: List of dictionaries containing parent feature details

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)

        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get the feature
        try:
            feature = db[feature_id]
        except gffutils.exceptions.FeatureNotFoundError:
            return f"Error: Feature '{feature_id}' not found in database."
        
        # Get all parent features
        parents = []
        for parent in db.parents(feature):
            parent_dict = {
                'id': parent.id,
                'chrom': parent.chrom,
                'start': parent.start,
                'end': parent.end,
                'strand': parent.strand,
                'feature_type': parent.featuretype,
                'attributes': dict(parent.attributes.items()),
                'length': abs(parent.end - parent.start)
            }
            parents.append(parent_dict)
        
        # Sort parents by hierarchical level (larger features typically higher in hierarchy)
        parents.sort(key=lambda x: x['length'], reverse=True)
        
        return parents
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error getting feature parents: {str(e)}"


@tool
def get_features_by_type(gffpath: str, feature_type: str, limit: int = None) -> list:
    """Get all features of a specific type using efficient iteration.

    Args:
        gffpath (str): Path to the GFF file
        feature_type (str): The feature type to query (e.g., 'gene', 'exon', 'CDS')
        limit (int, optional): Maximum number of features to return (for large datasets)

    Returns:
        list: List of dictionaries containing feature details

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)

        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Check if feature type exists in database
        available_types = list(db.featuretypes())
        if feature_type not in available_types:
            return f"Error: Feature type '{feature_type}' not found. Available types: {available_types}"
        
        # Get features of the specified type
        features = []
        count = 0
        
        for feature in db.features_of_type(feature_type):
            if limit and count >= limit:
                break
                
            feature_dict = {
                'id': feature.id,
                'chrom': feature.chrom,
                'start': feature.start,
                'end': feature.end,
                'strand': feature.strand,
                'feature_type': feature.featuretype,
                'attributes': dict(feature.attributes.items()),
                'length': abs(feature.end - feature.start)
            }
            features.append(feature_dict)
            count += 1
        
        # Sort features by chromosome and start position
        features.sort(key=lambda x: (x['chrom'], x['start']))
        
        return features
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error getting features by type: {str(e)}"


@tool
def get_feature_statistics(gffpath: str) -> dict:
    """Calculate comprehensive feature statistics including counts and length statistics per feature type.
    Use this function to get ammount of genes, exons, CDNs, chromosomes, and so on.

    Args:
        gffpath (str): Path to the GFF file

    Returns:
        dict: Dictionary containing comprehensive statistics
              Format: {
                  'total_features': int,
                  'feature_types': {
                      'gene': {'count': int, 'total_length': int, 'avg_length': float, 'min_length': int, 'max_length': int},
                      'exon': {'count': int, 'total_length': int, 'avg_length': float, 'min_length': int, 'max_length': int},
                      ...
                  },
                  'chromosomes': list,
                  'total_genome_length': int
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)

        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get all feature types
        feature_types = list(db.featuretypes())
        
        # Initialize statistics dictionary
        stats = {
            'total_features': 0,
            'feature_types': {},
            'chromosomes': [],
            'total_genome_length': 0
        }
        
        # Get chromosome list
        chromosomes = set()
        
        # Calculate statistics for each feature type
        for feature_type in feature_types:
            type_stats = {
                'count': 0,
                'total_length': 0,
                'lengths': []  # Temporary list to calculate min, max, avg
            }
            
            # Iterate through all features of this type
            for feature in db.features_of_type(feature_type):
                length = abs(feature.end - feature.start)
                type_stats['count'] += 1
                type_stats['total_length'] += length
                type_stats['lengths'].append(length)
                chromosomes.add(feature.chrom)
            
            # Calculate derived statistics
            if type_stats['count'] > 0:
                type_stats['avg_length'] = round(type_stats['total_length'] / type_stats['count'], 2)
                type_stats['min_length'] = min(type_stats['lengths'])
                type_stats['max_length'] = max(type_stats['lengths'])
            else:
                type_stats['avg_length'] = 0
                type_stats['min_length'] = 0
                type_stats['max_length'] = 0
            
            # Remove temporary lengths list
            del type_stats['lengths']
            
            # Add to main stats
            stats['feature_types'][feature_type] = type_stats
            stats['total_features'] += type_stats['count']
            stats['total_genome_length'] += type_stats['total_length']
        
        # Convert chromosomes set to sorted list
        stats['chromosomes'] = sorted(list(chromosomes))
        
        return stats
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating feature statistics: {str(e)}"



@tool
def get_chromosomes_info(gffpath: str) -> list:
    """Get all chromosome names, use this to get number of chromosomes.

    Args:
        gffpath (str): Path to the GFF file

    Returns:
        list: List with the name of all chromosomes

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get all chromosomes or filter to specific one
        all_chromosomes = set()
        for feature in db.all_features():
            all_chromosomes.add(feature.chrom)
        


        
        return list(all_chromosomes)
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating chromosome summary: {str(e)}"


@tool
def get_chromosome_summary(gffpath: str, chrom: str = None) -> dict:
    """Calculate per-chromosome feature analysis with counts and statistics.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str, optional): Specific chromosome to analyze. If None, analyzes all chromosomes.

    Returns:
        dict: Dictionary containing per-chromosome statistics
              Format: {
                  'chromosome_name': {
                      'total_features': int,
                      'feature_types': {
                          'gene': {'count': int, 'total_length': int, 'avg_length': float},
                          'exon': {'count': int, 'total_length': int, 'avg_length': float},
                          ...
                      },
                      'chromosome_length': int,
                      'feature_density': float  # features per kb
                  },
                  ...
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get all chromosomes or filter to specific one
        all_chromosomes = set()
        for feature in db.all_features():
            all_chromosomes.add(feature.chrom)
        
        if chrom:
            if chrom not in all_chromosomes:
                return f"Error: Chromosome '{chrom}' not found. Available chromosomes: {sorted(list(all_chromosomes))}"
            chromosomes_to_analyze = [chrom]
        else:
            chromosomes_to_analyze = sorted(list(all_chromosomes))
        
        # Initialize summary dictionary
        summary = {}
        
        # Analyze each chromosome
        for chromosome in chromosomes_to_analyze:
            chrom_stats = {
                'total_features': 0,
                'feature_types': {},
                'chromosome_length': 0,
                'feature_density': 0.0
            }
            
            # Get all features for this chromosome
            chromosome_features = {}  # feature_type -> list of features
            max_end = 0
            
            for feature in db.region(seqid=chromosome):
                feature_type = feature.featuretype
                
                if feature_type not in chromosome_features:
                    chromosome_features[feature_type] = []
                
                chromosome_features[feature_type].append(feature)
                chrom_stats['total_features'] += 1
                
                # Track chromosome length (maximum end coordinate)
                if feature.end > max_end:
                    max_end = feature.end
            
            chrom_stats['chromosome_length'] = max_end
            
            # Calculate statistics for each feature type on this chromosome
            for feature_type, features in chromosome_features.items():
                type_stats = {
                    'count': len(features),
                    'total_length': 0,
                    'avg_length': 0.0
                }
                
                # Calculate lengths
                lengths = []
                for feature in features:
                    length = abs(feature.end - feature.start)
                    lengths.append(length)
                    type_stats['total_length'] += length
                
                # Calculate average length
                if type_stats['count'] > 0:
                    type_stats['avg_length'] = round(type_stats['total_length'] / type_stats['count'], 2)
                
                chrom_stats['feature_types'][feature_type] = type_stats
            
            # Calculate feature density (features per kb)
            if chrom_stats['chromosome_length'] > 0:
                chrom_stats['feature_density'] = round(chrom_stats['total_features'] / (chrom_stats['chromosome_length'] / 1000), 2)
            
            summary[chromosome] = chrom_stats
        
        return summary
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating chromosome summary: {str(e)}"


@tool
def get_length_distribution(gffpath: str, feature_type: str) -> dict:
    """Calculate length statistics and distribution for a specific feature type.

    Args:
        gffpath (str): Path to the GFF file
        feature_type (str): The feature type to analyze (e.g., 'gene', 'exon', 'CDS')

    Returns:
        dict: Dictionary containing length distribution statistics
              Format: {
                  'feature_type': str,
                  'total_count': int,
                  'statistics': {
                      'min': int,
                      'max': int,
                      'mean': float,
                      'median': float,
                      'std_dev': float,
                      'total_length': int
                  },
                  'histogram': {
                      'bins': list,  # bin edges
                      'counts': list,  # count in each bin
                      'bin_width': int
                  },
                  'percentiles': {
                      '25th': float,
                      '75th': float,
                      '90th': float,
                      '95th': float
                  }
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)

        # Check if database exists, create if not
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Check if feature type exists
        available_types = list(db.featuretypes())
        if feature_type not in available_types:
            return f"Error: Feature type '{feature_type}' not found. Available types: {available_types}"
        
        # Collect all lengths for the specified feature type
        lengths = []
        for feature in db.features_of_type(feature_type):
            length = abs(feature.end - feature.start)
            lengths.append(length)
        
        if not lengths:
            return f"Error: No features of type '{feature_type}' found."
        
        # Sort lengths for percentile calculations
        lengths.sort()
        n = len(lengths)
        
        # Calculate basic statistics
        min_length = min(lengths)
        max_length = max(lengths)
        total_length = sum(lengths)
        mean_length = total_length / n
        
        # Calculate median
        if n % 2 == 0:
            median_length = (lengths[n//2 - 1] + lengths[n//2]) / 2
        else:
            median_length = lengths[n//2]
        
        # Calculate standard deviation
        variance = sum((x - mean_length) ** 2 for x in lengths) / n
        std_dev = variance ** 0.5
        
        # Calculate percentiles
        def percentile(data, p):
            index = int(p * len(data) / 100)
            if index >= len(data):
                index = len(data) - 1
            return data[index]
        
        percentiles = {
            '25th': percentile(lengths, 25),
            '75th': percentile(lengths, 75),
            '90th': percentile(lengths, 90),
            '95th': percentile(lengths, 95)
        }
        
        # Create histogram (10 bins)
        num_bins = min(10, n)  # Use fewer bins if we have very few features
        bin_width = (max_length - min_length) / num_bins if num_bins > 1 else 1
        
        # Initialize bins
        bins = []
        counts = []
        
        for i in range(num_bins):
            bin_start = min_length + i * bin_width
            bin_end = min_length + (i + 1) * bin_width
            bins.append(f"{int(bin_start)}-{int(bin_end)}")
            
            # Count features in this bin
            count = 0
            for length in lengths:
                if i == num_bins - 1:  # Last bin includes the maximum
                    if bin_start <= length <= bin_end:
                        count += 1
                else:
                    if bin_start <= length < bin_end:
                        count += 1
            counts.append(count)
        
        # Prepare result
        result = {
            'feature_type': feature_type,
            'total_count': n,
            'statistics': {
                'min': min_length,
                'max': max_length,
                'mean': round(mean_length, 2),
                'median': round(median_length, 2),
                'std_dev': round(std_dev, 2),
                'total_length': total_length
            },
            'histogram': {
                'bins': bins,
                'counts': counts,
                'bin_width': round(bin_width, 2)
            },
            'percentiles': {
                '25th': round(percentiles['25th'], 2),
                '75th': round(percentiles['75th'], 2),
                '90th': round(percentiles['90th'], 2),
                '95th': round(percentiles['95th'], 2)
            }
        }
        
        return result
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating length distribution: {str(e)}"


@tool
def search_features_by_attribute(gffpath: str, attribute_key: str, attribute_value: str, exact_match: bool = True) -> list:
    """Search features by attribute key-value pairs with exact or partial matching.

    Args:
        gffpath (str): Path to the GFF file
        attribute_key (str): The attribute key to search for (e.g., 'Name', 'ID', 'Note')
        attribute_value (str): The attribute value to match
        exact_match (bool): If True, performs exact matching; if False, performs partial matching

    Returns:
        list: List of dictionaries containing matching features with complete attribute information

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        matching_features = []
        
        # Iterate through all features to search for matching attributes
        for feature in db.all_features():
            # Check if the feature has the specified attribute key
            if attribute_key in feature.attributes:
                # Get the attribute values (attributes are stored as lists)
                attr_values = feature.attributes[attribute_key]
                
                # Check if any of the attribute values match our search criteria
                match_found = False
                for attr_val in attr_values:
                    if exact_match:
                        if attr_val == attribute_value:
                            match_found = True
                            break
                    else:
                        if attribute_value.lower() in attr_val.lower():
                            match_found = True
                            break
                
                if match_found:
                    feature_dict = {
                        'id': feature.id,
                        'chrom': feature.chrom,
                        'start': feature.start,
                        'end': feature.end,
                        'strand': feature.strand,
                        'feature_type': feature.featuretype,
                        'attributes': dict(feature.attributes.items()),
                        'length': abs(feature.end - feature.start)
                    }
                    matching_features.append(feature_dict)
        
        # Sort results by chromosome and start position
        matching_features.sort(key=lambda x: (x['chrom'], x['start']))
        
        return matching_features
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error searching features by attribute: {str(e)}"



#search_features_by_attribute("GCA_050947715.1_ASM5094771v1_genomic.gff", "go_function", "transmembrane", False)


@tool
def search_genes_by_go_function_attribute(gffpath: str, attribute_value: str, exact_match: bool = True) -> list:
    """Search genes by matching go function (GO: Gene Ontology) with exact or partial matching.
    Use this function when asked about genes that encodes a specific protein.

    Args:
        gffpath (str): Path to the GFF file
        attribute_value (str): The attribute value to match
        exact_match (bool): If True, performs exact matching; if False, performs partial matching

    Returns:
        list: List of string with the gene id

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        matching_features = []
        
        # Iterate through all features to search for matching attributes
        for feature in db.all_features():
            # Check if the feature has the specified attribute key
            if "go_function" in feature.attributes:
                # Get the attribute values (attributes are stored as lists)
                attr_values = feature.attributes["go_function"]
                
                # Check if any of the attribute values match our search criteria
                match_found = False
                for attr_val in attr_values:
                    if exact_match:
                        if attr_val == attribute_value:
                            match_found = True
                            break
                    else:
                        if attribute_value.lower() in attr_val.lower():
                            match_found = True
                            break
                
                if match_found:
                    #feature_dict = {
                    #    'id': feature.id,
                    #    'attributes': feature.attributes["go_function"],
                    #}
                    matching_features.append(feature.id)
        
        return matching_features
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error searching features by attribute: {str(e)}"


@tool
def get_features_with_attribute(gffpath: str, attribute_key: str) -> list:
    """Find all features that have a specific attribute key present.

    Args:
        gffpath (str): Path to the GFF file
        attribute_key (str): The attribute key to search for (e.g., 'Name', 'ID', 'Note')

    Returns:
        list: List of dictionaries containing features with the specified attribute

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        matching_features = []
        
        # Iterate through all features to find those with the specified attribute key
        for feature in db.all_features():
            if attribute_key in feature.attributes:
                feature_dict = {
                    'id': feature.id,
                    'chrom': feature.chrom,
                    'start': feature.start,
                    'end': feature.end,
                    'strand': feature.strand,
                    'feature_type': feature.featuretype,
                    'attributes': dict(feature.attributes.items()),
                    'length': abs(feature.end - feature.start)
                }
                matching_features.append(feature_dict)
        
        # Sort results by chromosome and start position
        matching_features.sort(key=lambda x: (x['chrom'], x['start']))
        
        return matching_features
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error finding features with attribute: {str(e)}"


@tool
def get_intergenic_regions(gffpath: str, chrom: str = None, min_length: int = 0) -> list:
    """Identify gaps between genes with filtering by minimum length and chromosome.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str, optional): Specific chromosome to analyze. If None, analyzes all chromosomes.
        min_length (int): Minimum length of intergenic regions to include (default: 0)

    Returns:
        list: List of dictionaries containing intergenic region information
              Format: [
                  {
                      'chrom': str,
                      'start': int,
                      'end': int,
                      'length': int,
                      'upstream_gene': str,
                      'downstream_gene': str
                  },
                  ...
              ]

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get all chromosomes or filter to specific one
        all_chromosomes = set()
        for feature in db.all_features():
            all_chromosomes.add(feature.chrom)
        
        if chrom:
            if chrom not in all_chromosomes:
                return f"Error: Chromosome '{chrom}' not found. Available chromosomes: {sorted(list(all_chromosomes))}"
            chromosomes_to_analyze = [chrom]
        else:
            chromosomes_to_analyze = sorted(list(all_chromosomes))
        
        intergenic_regions = []
        
        # Analyze each chromosome
        for chromosome in chromosomes_to_analyze:
            # Get all genes on this chromosome, sorted by start position
            genes = []
            for feature in db.region(seqid=chromosome):
                if feature.featuretype == 'gene':
                    genes.append(feature)
            
            # Sort genes by start position
            genes.sort(key=lambda x: x.start)
            
            # Find gaps between consecutive genes
            for i in range(len(genes) - 1):
                current_gene = genes[i]
                next_gene = genes[i + 1]
                
                # Calculate intergenic region
                intergenic_start = current_gene.end + 1
                intergenic_end = next_gene.start - 1
                intergenic_length = intergenic_end - intergenic_start + 1
                
                # Only include if it meets minimum length requirement
                if intergenic_length >= min_length and intergenic_length > 0:
                    intergenic_region = {
                        'chrom': chromosome,
                        'start': intergenic_start,
                        'end': intergenic_end,
                        'length': intergenic_length,
                        'upstream_gene': current_gene.id,
                        'downstream_gene': next_gene.id
                    }
                    intergenic_regions.append(intergenic_region)
        
        # Sort by chromosome and start position
        intergenic_regions.sort(key=lambda x: (x['chrom'], x['start']))
        
        return intergenic_regions
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error identifying intergenic regions: {str(e)}"


@tool
def get_feature_density(gffpath: str, chrom: str, window_size: int, feature_type: str = None) -> list:
    """Calculate feature density in genomic windows across a chromosome.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str): Chromosome name to analyze
        window_size (int): Size of genomic windows in base pairs
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon'). If None, counts all features.

    Returns:
        list: List of dictionaries containing density information for each window
              Format: [
                  {
                      'chrom': str,
                      'window_start': int,
                      'window_end': int,
                      'feature_count': int,
                      'density': float  # features per kb
                  },
                  ...
              ]

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Check if chromosome exists
        all_chromosomes = set()
        for feature in db.all_features():
            all_chromosomes.add(feature.chrom)
        
        if chrom not in all_chromosomes:
            return f"Error: Chromosome '{chrom}' not found. Available chromosomes: {sorted(list(all_chromosomes))}"
        
        # Get chromosome length (maximum end coordinate)
        max_end = 0
        for feature in db.region(seqid=chrom):
            if feature.end > max_end:
                max_end = feature.end
        
        if max_end == 0:
            return f"Error: No features found on chromosome '{chrom}'"
        
        # Create windows
        windows = []
        current_start = 1
        
        while current_start <= max_end:
            window_end = min(current_start + window_size - 1, max_end)
            
            # Count features in this window
            feature_count = 0
            for feature in db.region(seqid=chrom, start=current_start, end=window_end):
                # Apply feature type filter if specified
                if feature_type is None or feature.featuretype == feature_type:
                    feature_count += 1
            
            # Calculate density (features per kb)
            actual_window_size = window_end - current_start + 1
            density = (feature_count / actual_window_size) * 1000  # per kb
            
            window_data = {
                'chrom': chrom,
                'window_start': current_start,
                'window_end': window_end,
                'feature_count': feature_count,
                'density': round(density, 4)
            }
            windows.append(window_data)
            
            current_start += window_size
        
        return windows
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating feature density: {str(e)}"


@tool
def get_strand_distribution(gffpath: str, feature_type: str = None) -> dict:
    """Analyze strand distribution of features with counts and percentages.

    Args:
        gffpath (str): Path to the GFF file
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon'). If None, analyzes all features.

    Returns:
        dict: Dictionary containing strand distribution analysis
              Format: {
                  'feature_type': str,
                  'total_features': int,
                  'strand_counts': {'+': int, '-': int, '.': int, 'other': int},
                  'strand_percentages': {'+': float, '-': float, '.': float, 'other': float},
                  'chromosome_breakdown': {
                      'chr1': {
                          'total': int,
                          'strand_counts': {'+': int, '-': int, '.': int, 'other': int},
                          'strand_percentages': {'+': float, '-': float, '.': float, 'other': float}
                      },
                      ...
                  }
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Check if feature type exists (if specified)
        if feature_type:
            available_types = list(db.featuretypes())
            if feature_type not in available_types:
                return f"Error: Feature type '{feature_type}' not found. Available types: {available_types}"
        
        # Initialize counters
        strand_counts = {'+': 0, '-': 0, '.': 0, 'other': 0}
        chromosome_data = {}
        total_features = 0
        
        # Determine which features to analyze
        if feature_type:
            features_to_analyze = db.features_of_type(feature_type)
            analysis_type = feature_type
        else:
            features_to_analyze = db.all_features()
            analysis_type = 'all'
        
        # Count strand distribution
        for feature in features_to_analyze:
            total_features += 1
            strand = feature.strand
            chrom = feature.chrom
            
            # Count overall strand distribution
            if strand in ['+', '-', '.']:
                strand_counts[strand] += 1
            else:
                strand_counts['other'] += 1
            
            # Initialize chromosome data if not exists
            if chrom not in chromosome_data:
                chromosome_data[chrom] = {
                    'total': 0,
                    'strand_counts': {'+': 0, '-': 0, '.': 0, 'other': 0}
                }
            
            # Count per chromosome
            chromosome_data[chrom]['total'] += 1
            if strand in ['+', '-', '.']:
                chromosome_data[chrom]['strand_counts'][strand] += 1
            else:
                chromosome_data[chrom]['strand_counts']['other'] += 1
        
        # Calculate overall percentages
        strand_percentages = {}
        for strand, count in strand_counts.items():
            if total_features > 0:
                strand_percentages[strand] = round((count / total_features) * 100, 2)
            else:
                strand_percentages[strand] = 0.0
        
        # Calculate per-chromosome percentages
        for chrom in chromosome_data:
            chrom_total = chromosome_data[chrom]['total']
            chromosome_data[chrom]['strand_percentages'] = {}
            
            for strand, count in chromosome_data[chrom]['strand_counts'].items():
                if chrom_total > 0:
                    chromosome_data[chrom]['strand_percentages'][strand] = round((count / chrom_total) * 100, 2)
                else:
                    chromosome_data[chrom]['strand_percentages'][strand] = 0.0
        
        # Prepare result
        result = {
            'feature_type': analysis_type,
            'total_features': total_features,
            'strand_counts': strand_counts,
            'strand_percentages': strand_percentages,
            'chromosome_breakdown': chromosome_data
        }
        
        return result
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error analyzing strand distribution: {str(e)}"


@tool
def export_features_to_csv(gffpath: str, output_path: str, feature_type: str = None, chrom: str = None) -> str:
    """Export feature data to CSV format with filtering options.

    Args:
        gffpath (str): Path to the GFF file
        output_path (str): Path for the output CSV file
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon')
        chrom (str, optional): Filter by chromosome

    Returns:
        str: Success message with export details or error message

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Collect features based on filters
        features_to_export = []
        
        if feature_type and chrom:
            # Filter by both feature type and chromosome
            for feature in db.features_of_type(feature_type):
                if feature.chrom == chrom:
                    features_to_export.append(feature)
        elif feature_type:
            # Filter by feature type only
            features_to_export = list(db.features_of_type(feature_type))
        elif chrom:
            # Filter by chromosome only
            features_to_export = list(db.region(seqid=chrom))
        else:
            # No filters - export all features
            features_to_export = list(db.all_features())
        
        if not features_to_export:
            return f"No features found matching the specified criteria."
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        
        # Write to CSV
        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            # Define CSV headers
            fieldnames = ['id', 'chrom', 'start', 'end', 'strand', 'feature_type', 'length']
            
            # Get all unique attribute keys to include as columns
            all_attributes = set()
            for feature in features_to_export:
                all_attributes.update(feature.attributes.keys())
            
            # Add attribute columns to fieldnames
            fieldnames.extend(sorted(all_attributes))
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write feature data
            for feature in features_to_export:
                row = {
                    'id': feature.id,
                    'chrom': feature.chrom,
                    'start': feature.start,
                    'end': feature.end,
                    'strand': feature.strand,
                    'feature_type': feature.featuretype,
                    'length': abs(feature.end - feature.start)
                }
                
                # Add attribute values (join lists with semicolons)
                for attr_key in all_attributes:
                    if attr_key in feature.attributes:
                        attr_values = feature.attributes[attr_key]
                        row[attr_key] = ';'.join(attr_values) if isinstance(attr_values, list) else str(attr_values)
                    else:
                        row[attr_key] = ''
                
                writer.writerow(row)
        
        # Prepare summary message
        filter_info = []
        if feature_type:
            filter_info.append(f"feature_type='{feature_type}'")
        if chrom:
            filter_info.append(f"chromosome='{chrom}'")
        
        filter_str = f" (filtered by: {', '.join(filter_info)})" if filter_info else ""
        
        return f"Successfully exported {len(features_to_export)} features to '{output_path}'{filter_str}."
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error exporting features to CSV: {str(e)}"


@tool
def get_feature_summary_report(gffpath: str) -> str:
    """Generate a human-readable summary report of GFF file contents.

    Args:
        gffpath (str): Path to the GFF file

    Returns:
        str: Formatted text report with comprehensive GFF file summary

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        db_filename = get_db_filename(gffpath)
        if os.path.exists(db_filename):
            db = gffutils.FeatureDB(db_filename)
        else:
            db = gffutils.create_db(gffpath, dbfn=db_filename, force=True, keep_order=True, merge_strategy="create_unique")
        
        # Get basic statistics
        feature_types = list(db.featuretypes())
        total_features = sum(1 for _ in db.all_features())
        
        # Get chromosomes
        chromosomes = set()
        for feature in db.all_features():
            chromosomes.add(feature.chrom)
        chromosomes = sorted(list(chromosomes))
        
        # Count features by type
        feature_counts = {}
        for feature_type in feature_types:
            count = sum(1 for _ in db.features_of_type(feature_type))
            feature_counts[feature_type] = count
        
        # Calculate some basic length statistics for genes (if available)
        gene_stats = ""
        if 'gene' in feature_types:
            gene_lengths = []
            for gene in db.features_of_type('gene'):
                gene_lengths.append(abs(gene.end - gene.start))
            
            if gene_lengths:
                gene_stats = f"""
Gene Length Statistics:
  - Total genes: {len(gene_lengths):,}
  - Average length: {sum(gene_lengths) / len(gene_lengths):,.0f} bp
  - Shortest gene: {min(gene_lengths):,} bp
  - Longest gene: {max(gene_lengths):,} bp
"""
        
        # Build the report
        report = f"""
=== GFF FILE SUMMARY REPORT ===

File: {gffpath}
Total Features: {total_features:,}
Feature Types: {len(feature_types)}
Chromosomes: {len(chromosomes)}

CHROMOSOMES:
{', '.join(chromosomes)}

FEATURE TYPE BREAKDOWN:
"""
        
        # Add feature type counts (sorted by count, descending)
        sorted_features = sorted(feature_counts.items(), key=lambda x: x[1], reverse=True)
        for feature_type, count in sorted_features:
            percentage = (count / total_features) * 100 if total_features > 0 else 0
            report += f"  - {feature_type}: {count:,} ({percentage:.1f}%)\n"
        
        # Add gene statistics if available
        if gene_stats:
            report += gene_stats
        
        # Add chromosome breakdown for genes (if available)
        if 'gene' in feature_types and len(chromosomes) > 1:
            report += "\nGENES PER CHROMOSOME:\n"
            for chrom in chromosomes:
                gene_count = sum(1 for gene in db.features_of_type('gene') if gene.chrom == chrom)
                report += f"  - {chrom}: {gene_count:,} genes\n"
        
        # Add some attribute information
        all_attributes = set()
        sample_size = min(1000, total_features)  # Sample first 1000 features for attributes
        count = 0
        for feature in db.all_features():
            all_attributes.update(feature.attributes.keys())
            count += 1
            if count >= sample_size:
                break
        
        if all_attributes:
            report += f"\nCOMMON ATTRIBUTES ({len(all_attributes)} total):\n"
            # Show most common attributes (up to 10)
            common_attrs = sorted(list(all_attributes))[:10]
            report += f"  {', '.join(common_attrs)}\n"
            if len(all_attributes) > 10:
                report += f"  ... and {len(all_attributes) - 10} more\n"
        
        report += "\n=== END REPORT ===\n"
        
        return report
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error generating summary report: {str(e)}"


@tool
def get_tools_list() -> list:
    """Get a list of all available GFF analysis tools with their descriptions.

    Returns:
        list: List of tuples containing (function_name, description) for each tool
    """
    tools_info = [
        ("file_read", "Read a file and return its content"),
        ("file_write", "Write content to a file"),
        ("list_directory", "List files and directories in the specified path"),
        ("get_gff_feature_types", "Get all available feature types from a GFF file"),
        ("get_gene_lenght", "Get the length of a specific gene"),
        ("get_gene_attributes", "Get gene attributes (ID, Note, Name, etc.)"),
        ("get_multiple_gene_lenght", "Get lengths of multiple genes"),
        ("get_all_attributes", "Get all available attributes from a GFF file"),
        ("get_protein_product_from_gene", "Get protein products for a specific gene"),
        ("get_features_in_region", "Find features overlapping a genomic region"),
        ("get_features_at_position", "Find features at a specific genomic position"),
        ("get_gene_structure", "Get gene structure with child features (exons, CDS, UTRs)"),
        ("get_feature_parents", "Find parent features of a given feature"),
        ("get_features_by_type", "Get all features of a specific type"),
        ("get_feature_statistics", "Calculate comprehensive feature statistics. "
        "Use this function to get ammount of genes, exons, CDNs, chromosomes, and so on"),
        ("get_chromosome_summary", "Generate per-chromosome feature analysis"),
        ("get_length_distribution", "Calculate length statistics for a feature type"),
        ("search_features_by_attribute", "Search features by attribute key-value pairs"),
        ("get_features_with_attribute", "Find features with a specific attribute key"),
        ("get_intergenic_regions", "Identify gaps between genes"),
        ("get_feature_density", "Calculate feature density in genomic windows"),
        ("get_strand_distribution", "Analyze strand distribution of features"),
        ("export_features_to_csv", "Export feature data to CSV format"),
        ("get_feature_summary_report", "Generate human-readable GFF summary report"),
        ("get_tools_list", "Get list of all available tools with descriptions")
    ]
    
    return tools_info