

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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
                
                # If we found a match, add this feature to results
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        matching_features = []
        
        # Iterate through all features to find those with the specified attribute
        for feature in db.all_features():
            # Check if the feature has the specified attribute key
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
        list: List of dictionaries containing intergenic region coordinates and lengths
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
            for gene in db.features_of_type('gene'):
                if gene.chrom == chromosome:
                    genes.append({
                        'id': gene.id,
                        'start': gene.start,
                        'end': gene.end,
                        'strand': gene.strand
                    })
            
            # Sort genes by start position
            genes.sort(key=lambda x: x['start'])
            
            # Find intergenic regions between consecutive genes
            for i in range(len(genes) - 1):
                current_gene = genes[i]
                next_gene = genes[i + 1]
                
                # Calculate intergenic region
                intergenic_start = current_gene['end'] + 1
                intergenic_end = next_gene['start'] - 1
                intergenic_length = intergenic_end - intergenic_start + 1
                
                # Only include if it meets minimum length requirement and is positive
                if intergenic_length >= min_length and intergenic_length > 0:
                    intergenic_region = {
                        'chrom': chromosome,
                        'start': intergenic_start,
                        'end': intergenic_end,
                        'length': intergenic_length,
                        'upstream_gene': current_gene['id'],
                        'downstream_gene': next_gene['id']
                    }
                    intergenic_regions.append(intergenic_region)
        
        # Sort results by chromosome and start position
        intergenic_regions.sort(key=lambda x: (x['chrom'], x['start']))
        
        return intergenic_regions
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error finding intergenic regions: {str(e)}"


@tool
def get_feature_density(gffpath: str, chrom: str, window_size: int, feature_type: str = None) -> list:
    """Calculate feature density in genomic windows across a chromosome.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str): Chromosome name to analyze
        window_size (int): Size of genomic windows in base pairs
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon')

    Returns:
        list: List of dictionaries containing density values across chromosome
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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        # Check if chromosome exists
        all_chromosomes = set()
        for feature in db.all_features():
            all_chromosomes.add(feature.chrom)
        
        if chrom not in all_chromosomes:
            return f"Error: Chromosome '{chrom}' not found. Available chromosomes: {sorted(list(all_chromosomes))}"
        
        # Find the maximum coordinate on this chromosome to determine chromosome length
        max_coord = 0
        for feature in db.region(seqid=chrom):
            if feature.end > max_coord:
                max_coord = feature.end
        
        if max_coord == 0:
            return f"Error: No features found on chromosome '{chrom}'"
        
        # Create windows and calculate density
        density_data = []
        window_start = 1
        
        while window_start <= max_coord:
            window_end = min(window_start + window_size - 1, max_coord)
            
            # Count features in this window
            feature_count = 0
            for feature in db.region(seqid=chrom, start=window_start, end=window_end):
                # Apply feature type filter if specified
                if feature_type is None or feature.featuretype == feature_type:
                    feature_count += 1
            
            # Calculate density (features per kb)
            actual_window_size = window_end - window_start + 1
            density = (feature_count / actual_window_size) * 1000  # features per kb
            
            window_data = {
                'chrom': chrom,
                'window_start': window_start,
                'window_end': window_end,
                'feature_count': feature_count,
                'density': round(density, 4)
            }
            density_data.append(window_data)
            
            window_start += window_size
        
        return density_data
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating feature density: {str(e)}"


@tool
def get_strand_distribution(gffpath: str, feature_type: str = None) -> dict:
    """Analyze strand distribution of features with counts and percentages.

    Args:
        gffpath (str): Path to the GFF file
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon')

    Returns:
        dict: Dictionary containing strand distribution analysis
              Format: {
                  'feature_type': str or 'all',
                  'total_features': int,
                  'strand_counts': {
                      '+': int,
                      '-': int,
                      '.': int,  # unstranded
                      'other': int  # any other strand values
                  },
                  'strand_percentages': {
                      '+': float,
                      '-': float,
                      '.': float,
                      'other': float
                  },
                  'strand_balance': float  # ratio of + to - strands
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        # Initialize strand counts
        strand_counts = {'+': 0, '-': 0, '.': 0, 'other': 0}
        
        # Count features by strand
        if feature_type:
            # Check if feature type exists
            available_types = list(db.featuretypes())
            if feature_type not in available_types:
                return f"Error: Feature type '{feature_type}' not found. Available types: {available_types}"
            
            # Count for specific feature type
            for feature in db.features_of_type(feature_type):
                strand = feature.strand
                if strand in strand_counts:
                    strand_counts[strand] += 1
                else:
                    strand_counts['other'] += 1
        else:
            # Count for all features
            for feature in db.all_features():
                strand = feature.strand
                if strand in strand_counts:
                    strand_counts[strand] += 1
                else:
                    strand_counts['other'] += 1
        
        # Calculate total and percentages
        total_features = sum(strand_counts.values())
        
        if total_features == 0:
            return f"Error: No features found" + (f" of type '{feature_type}'" if feature_type else "")
        
        strand_percentages = {}
        for strand, count in strand_counts.items():
            strand_percentages[strand] = round((count / total_features) * 100, 2)
        
        # Calculate strand balance (+ to - ratio)
        plus_count = strand_counts['+']
        minus_count = strand_counts['-']
        
        if minus_count > 0:
            strand_balance = round(plus_count / minus_count, 2)
        elif plus_count > 0:
            strand_balance = float('inf')  # All plus, no minus
        else:
            strand_balance = 0  # No stranded features
        
        result = {
            'feature_type': feature_type if feature_type else 'all',
            'total_features': total_features,
            'strand_counts': strand_counts,
            'strand_percentages': strand_percentages,
            'strand_balance': strand_balance
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
        output_path (str): Path where the CSV file will be saved
        feature_type (str, optional): Filter by feature type (e.g., 'gene', 'exon')
        chrom (str, optional): Filter by chromosome

    Returns:
        str: Success message with export details or error message

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        import csv
        
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        # Collect features based on filters
        features_to_export = []
        
        if feature_type and chrom:
            # Filter by both feature type and chromosome
            for feature in db.features_of_type(feature_type):
                if feature.chrom == chrom:
                    features_to_export.append(feature)
        elif feature_type:
            # Filter by feature type only
            for feature in db.features_of_type(feature_type):
                features_to_export.append(feature)
        elif chrom:
            # Filter by chromosome only
            for feature in db.region(seqid=chrom):
                features_to_export.append(feature)
        else:
            # Export all features
            for feature in db.all_features():
                features_to_export.append(feature)
        
        if not features_to_export:
            filter_desc = []
            if feature_type:
                filter_desc.append(f"feature type '{feature_type}'")
            if chrom:
                filter_desc.append(f"chromosome '{chrom}'")
            filter_str = " and ".join(filter_desc) if filter_desc else "specified criteria"
            return f"No features found matching {filter_str}"
        
        # Prepare CSV data
        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            # Define CSV columns
            fieldnames = [
                'id', 'chromosome', 'start', 'end', 'strand', 'feature_type', 'length'
            ]
            
            # Collect all unique attribute keys to include as columns
            all_attributes = set()
            for feature in features_to_export:
                all_attributes.update(feature.attributes.keys())
            
            # Add attribute columns (sorted for consistency)
            attribute_columns = sorted(list(all_attributes))
            fieldnames.extend(attribute_columns)
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write feature data
            for feature in features_to_export:
                row = {
                    'id': feature.id,
                    'chromosome': feature.chrom,
                    'start': feature.start,
                    'end': feature.end,
                    'strand': feature.strand,
                    'feature_type': feature.featuretype,
                    'length': abs(feature.end - feature.start)
                }
                
                # Add attribute values (join multiple values with semicolon)
                for attr_key in attribute_columns:
                    if attr_key in feature.attributes:
                        attr_values = feature.attributes[attr_key]
                        row[attr_key] = ';'.join(attr_values) if isinstance(attr_values, list) else str(attr_values)
                    else:
                        row[attr_key] = ''
                
                writer.writerow(row)
        
        # Prepare success message
        filter_desc = []
        if feature_type:
            filter_desc.append(f"feature type '{feature_type}'")
        if chrom:
            filter_desc.append(f"chromosome '{chrom}'")
        
        filter_str = f" (filtered by {' and '.join(filter_desc)})" if filter_desc else ""
        
        return f"Successfully exported {len(features_to_export)} features to '{output_path}'{filter_str}. Columns include: {', '.join(fieldnames)}"
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except PermissionError:
        return f"Error: Permission denied writing to '{output_path}'. Check file permissions and path."
    except Exception as e:
        return f"Error exporting features to CSV: {str(e)}"


@tool
def get_feature_summary_report(gffpath: str) -> str:
    """Generate a human-readable summary report of GFF file contents.

    Args:
        gffpath (str): Path to the GFF file

    Returns:
        str: Formatted text report with key statistics, feature type breakdown, and genome overview

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        # Get comprehensive statistics
        stats = get_feature_statistics(gffpath)
        if isinstance(stats, str) and stats.startswith("Error"):
            return stats
        
        # Get chromosome information
        chromosomes = stats['chromosomes']
        
        # Start building the report
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("GFF FILE ANALYSIS SUMMARY REPORT")
        report_lines.append("=" * 80)
        report_lines.append(f"File: {gffpath}")
        report_lines.append(f"Generated: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append("")
        
        # Overall Statistics
        report_lines.append("OVERALL STATISTICS")
        report_lines.append("-" * 40)
        report_lines.append(f"Total Features: {stats['total_features']:,}")
        report_lines.append(f"Total Genome Length: {stats['total_genome_length']:,} bp")
        report_lines.append(f"Number of Chromosomes: {len(chromosomes)}")
        report_lines.append(f"Feature Types: {len(stats['feature_types'])}")
        report_lines.append("")
        
        # Chromosome Overview
        report_lines.append("CHROMOSOME OVERVIEW")
        report_lines.append("-" * 40)
        if len(chromosomes) <= 10:
            report_lines.append("Chromosomes: " + ", ".join(chromosomes))
        else:
            report_lines.append(f"Chromosomes: {', '.join(chromosomes[:5])}, ... and {len(chromosomes)-5} more")
        report_lines.append("")
        
        # Feature Type Breakdown
        report_lines.append("FEATURE TYPE BREAKDOWN")
        report_lines.append("-" * 40)
        report_lines.append(f"{'Feature Type':<20} {'Count':<10} {'Avg Length':<12} {'Total Length':<15}")
        report_lines.append("-" * 57)
        
        # Sort feature types by count (descending)
        sorted_types = sorted(stats['feature_types'].items(), key=lambda x: x[1]['count'], reverse=True)
        
        for feature_type, type_stats in sorted_types:
            count = type_stats['count']
            avg_length = type_stats['avg_length']
            total_length = type_stats['total_length']
            
            report_lines.append(f"{feature_type:<20} {count:<10,} {avg_length:<12.1f} {total_length:<15,}")
        
        report_lines.append("")
        
        # Top Feature Types (most abundant)
        report_lines.append("TOP 5 MOST ABUNDANT FEATURE TYPES")
        report_lines.append("-" * 40)
        top_5_types = sorted_types[:5]
        for i, (feature_type, type_stats) in enumerate(top_5_types, 1):
            percentage = (type_stats['count'] / stats['total_features']) * 100
            report_lines.append(f"{i}. {feature_type}: {type_stats['count']:,} features ({percentage:.1f}%)")
        report_lines.append("")
        
        # Length Statistics for Key Feature Types
        key_types = ['gene', 'exon', 'CDS', 'mRNA']
        available_key_types = [ft for ft in key_types if ft in stats['feature_types']]
        
        if available_key_types:
            report_lines.append("LENGTH STATISTICS FOR KEY FEATURE TYPES")
            report_lines.append("-" * 40)
            report_lines.append(f"{'Type':<10} {'Min':<8} {'Max':<10} {'Average':<10} {'Total':<12}")
            report_lines.append("-" * 50)
            
            for feature_type in available_key_types:
                type_stats = stats['feature_types'][feature_type]
                report_lines.append(f"{feature_type:<10} {type_stats['min_length']:<8,} {type_stats['max_length']:<10,} "
                                  f"{type_stats['avg_length']:<10.1f} {type_stats['total_length']:<12,}")
            report_lines.append("")
        
        # Genome Composition Analysis
        if 'gene' in stats['feature_types']:
            gene_stats = stats['feature_types']['gene']
            gene_coverage = (gene_stats['total_length'] / stats['total_genome_length']) * 100
            
            report_lines.append("GENOME COMPOSITION ANALYSIS")
            report_lines.append("-" * 40)
            report_lines.append(f"Gene Coverage: {gene_coverage:.2f}% of genome")
            report_lines.append(f"Average Gene Density: {gene_stats['count'] / (stats['total_genome_length'] / 1000000):.1f} genes per Mb")
            
            if 'exon' in stats['feature_types']:
                exon_stats = stats['feature_types']['exon']
                exon_coverage = (exon_stats['total_length'] / stats['total_genome_length']) * 100
                report_lines.append(f"Exon Coverage: {exon_coverage:.2f}% of genome")
                
                if gene_stats['total_length'] > 0:
                    exon_gene_ratio = (exon_stats['total_length'] / gene_stats['total_length']) * 100
                    report_lines.append(f"Exonic Content of Genes: {exon_gene_ratio:.2f}%")
            
            report_lines.append("")
        
        # Data Quality Indicators
        report_lines.append("DATA QUALITY INDICATORS")
        report_lines.append("-" * 40)
        
        # Check for common feature types
        expected_types = ['gene', 'mRNA', 'exon', 'CDS']
        present_types = [ft for ft in expected_types if ft in stats['feature_types']]
        missing_types = [ft for ft in expected_types if ft not in stats['feature_types']]
        
        report_lines.append(f"Standard Feature Types Present: {', '.join(present_types) if present_types else 'None'}")
        if missing_types:
            report_lines.append(f"Standard Feature Types Missing: {', '.join(missing_types)}")
        
        # Check for features with attributes
        features_with_id = 0
        features_with_name = 0
        total_checked = 0
        
        # Sample first 1000 features to check attribute presence
        for feature in db.all_features():
            if total_checked >= 1000:
                break
            if 'ID' in feature.attributes:
                features_with_id += 1
            if 'Name' in feature.attributes:
                features_with_name += 1
            total_checked += 1
        
        if total_checked > 0:
            id_percentage = (features_with_id / total_checked) * 100
            name_percentage = (features_with_name / total_checked) * 100
            report_lines.append(f"Features with ID attribute: {id_percentage:.1f}% (sampled {total_checked} features)")
            report_lines.append(f"Features with Name attribute: {name_percentage:.1f}% (sampled {total_checked} features)")
        
        report_lines.append("")
        
        # Footer
        report_lines.append("=" * 80)
        report_lines.append("End of Summary Report")
        report_lines.append("=" * 80)
        
        return "\n".join(report_lines)
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error generating summary report: {str(e)}"


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
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
    """Identify gaps between genes (intergenic regions) with filtering options.

    Args:
        gffpath (str): Path to the GFF file
        chrom (str, optional): Specific chromosome to analyze. If None, analyzes all chromosomes.
        min_length (int): Minimum length threshold for intergenic regions (default: 0)

    Returns:
        list: List of dictionaries containing intergenic region information
              Format: [
                  {
                      'chrom': str,
                      'start': int,
                      'end': int,
                      'length': int,
                      'upstream_gene': str,  # ID of gene before the region
                      'downstream_gene': str  # ID of gene after the region
                  },
                  ...
              ]

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
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
            for gene in db.features_of_type('gene'):
                if gene.chrom == chromosome:
                    genes.append({
                        'id': gene.id,
                        'start': gene.start,
                        'end': gene.end,
                        'strand': gene.strand
                    })
            
            # Sort genes by start position
            genes.sort(key=lambda x: x['start'])
            
            # Find intergenic regions between consecutive genes
            for i in range(len(genes) - 1):
                current_gene = genes[i]
                next_gene = genes[i + 1]
                
                # Calculate intergenic region coordinates
                # Start after current gene ends, end before next gene starts
                intergenic_start = current_gene['end'] + 1
                intergenic_end = next_gene['start'] - 1
                
                # Only include if there's actually a gap
                if intergenic_end >= intergenic_start:
                    intergenic_length = intergenic_end - intergenic_start + 1
                    
                    # Apply minimum length filter
                    if intergenic_length >= min_length:
                        intergenic_region = {
                            'chrom': chromosome,
                            'start': intergenic_start,
                            'end': intergenic_end,
                            'length': intergenic_length,
                            'upstream_gene': current_gene['id'],
                            'downstream_gene': next_gene['id']
                        }
                        intergenic_regions.append(intergenic_region)
        
        # Sort results by chromosome and start position
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
        feature_type (str, optional): Filter by specific feature type (e.g., 'gene', 'exon')

    Returns:
        list: List of dictionaries containing density information for each window
              Format: [
                  {
                      'chrom': str,
                      'window_start': int,
                      'window_end': int,
                      'window_number': int,
                      'feature_count': int,
                      'density': float,  # features per kb
                      'features': list  # list of feature IDs in this window
                  },
                  ...
              ]

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        # Validate chromosome exists
        all_chromosomes = set()
        for feature in db.all_features():
            all_chromosomes.add(feature.chrom)
        
        if chrom not in all_chromosomes:
            return f"Error: Chromosome '{chrom}' not found. Available chromosomes: {sorted(list(all_chromosomes))}"
        
        # Find the maximum coordinate on this chromosome to determine chromosome length
        max_coord = 0
        for feature in db.region(seqid=chrom):
            if feature.end > max_coord:
                max_coord = feature.end
        
        if max_coord == 0:
            return f"Error: No features found on chromosome '{chrom}'"
        
        # Calculate number of windows needed
        num_windows = (max_coord // window_size) + (1 if max_coord % window_size > 0 else 0)
        
        density_data = []
        
        # Analyze each window
        for window_num in range(num_windows):
            window_start = window_num * window_size + 1
            window_end = min((window_num + 1) * window_size, max_coord)
            
            # Count features in this window
            features_in_window = []
            feature_count = 0
            
            for feature in db.region(seqid=chrom, start=window_start, end=window_end):
                # Apply feature type filter if specified
                if feature_type and feature.featuretype != feature_type:
                    continue
                
                # Check if feature overlaps with window
                if not (feature.end < window_start or feature.start > window_end):
                    features_in_window.append(feature.id)
                    feature_count += 1
            
            # Calculate density (features per kb)
            window_length_kb = (window_end - window_start + 1) / 1000
            density = feature_count / window_length_kb if window_length_kb > 0 else 0
            
            window_data = {
                'chrom': chrom,
                'window_start': window_start,
                'window_end': window_end,
                'window_number': window_num + 1,
                'feature_count': feature_count,
                'density': round(density, 3),
                'features': features_in_window
            }
            
            density_data.append(window_data)
        
        return density_data
        
    except FileNotFoundError:
        return f"Error: File '{gffpath}' not found."
    except Exception as e:
        return f"Error calculating feature density: {str(e)}"


@tool
def get_strand_distribution(gffpath: str, feature_type: str = None) -> dict:
    """Analyze strand distribution of features with counts and percentages.

    Args:
        gffpath (str): Path to the GFF file
        feature_type (str, optional): Filter by specific feature type (e.g., 'gene', 'exon')

    Returns:
        dict: Dictionary containing strand distribution analysis
              Format: {
                  'feature_type': str,  # 'all' if no filter applied
                  'total_features': int,
                  'strand_counts': {
                      '+': int,
                      '-': int,
                      '.': int,  # unstranded
                      'other': int  # any other strand values
                  },
                  'strand_percentages': {
                      '+': float,
                      '-': float,
                      '.': float,
                      'other': float
                  },
                  'chromosome_breakdown': {
                      'chrom1': {
                          'total': int,
                          'strand_counts': {...},
                          'strand_percentages': {...}
                      },
                      ...
                  }
              }

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        # Check if database exists, create if not
        if os.path.exists('annotation.db'):
            db = gffutils.FeatureDB('annotation.db')
        else:
            db = gffutils.create_db(gffpath, dbfn='annotation.db', force=True, keep_order=True)
        
        # Validate feature type if specified
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
You are a helpful personal assistant with comprehensive bioinformatics capabilities for GFF file analysis and can perform local file actions and simple tasks for the user.

Your key capabilities:

## File Operations:
1. Read, understand, and summarize files
2. Create and write to files
3. List directory contents and provide information on files
4. Export data in various formats (CSV, JSON, TSV)

## GFF Analysis - Coordinate-based Queries:
5. Find features by genomic coordinates (regions and specific positions)
6. Query features overlapping genomic regions with filtering by type and strand
7. Identify features containing specific genomic positions

## GFF Analysis - Relationship and Hierarchy Queries:
8. Explore gene structure and organization (get all child features like exons, CDS, UTRs)
9. Find parent features of any given feature using upward traversal
10. Get all features of specific types with efficient iteration

## GFF Analysis - Statistical Analysis:
11. Calculate comprehensive feature statistics (counts, length distributions per feature type)
12. Generate per-chromosome feature summaries and analysis
13. Analyze length distributions with detailed statistics (min, max, mean, median, std dev, percentiles)
14. Create histogram data for feature length distributions

## GFF Analysis - Attribute-based Searches:
15. Search features by attribute key-value pairs (exact and partial matching)
16. Find features containing specific attribute keys
17. Support pattern matching and logical operations for attribute queries

## GFF Analysis - Positional Analysis:
18. Identify intergenic regions (gaps between genes) with filtering options
19. Calculate feature density in genomic windows across chromosomes
20. Analyze strand distribution of features with counts and percentages
21. Support clustering analysis and positional insights

## GFF Analysis - Export and Reporting:
22. Export feature data to CSV format with comprehensive filtering
23. Generate human-readable summary reports of GFF file contents
24. Provide formatted output for downstream analysis

Example queries you can handle:
- "Find all genes in chromosome 1 between positions 1000-5000"
- "Get the structure of gene AT1G01010 including all exons and CDS"
- "Calculate feature statistics for this GFF file"
- "Find all features with 'kinase' in their Name attribute"
- "Identify intergenic regions longer than 1000bp on chromosome 2"
- "Calculate gene density in 10kb windows across chromosome 1"
- "Export all exon features to CSV format"

When using tools:
- Always verify file paths before operations
- Provide clear explanations of what you're doing
- If a task cannot be completed, explain why and suggest alternatives
- Use appropriate filtering and limiting for large datasets
- Handle database creation and reuse efficiently

Always be helpful, concise, and focus on addressing the user's bioinformatics analysis needs efficiently.
"""

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
           get_multiple_gene_lenght, get_gene_attributes, get_features_in_region, get_features_at_position,
           get_gene_structure, get_feature_parents, get_features_by_type, get_feature_statistics,
           get_chromosome_summary, get_length_distribution, search_features_by_attribute, 
           get_features_with_attribute, get_intergenic_regions, get_feature_density, get_strand_distribution,
           export_features_to_csv, get_feature_summary_report],
)
"""



def main():
    # Define a comprehensive system prompt for our agent
    system_prompt = """
You are a helpful personal assistant with comprehensive bioinformatics capabilities for GFF file analysis and can perform local file actions and simple tasks for the user.

Your key capabilities:

## File Operations:
1. Read, understand, and summarize files
2. Create and write to files
3. List directory contents and provide information on files
4. Export data in various formats (CSV, JSON, TSV)

## GFF Analysis - Coordinate-based Queries:
5. Find features by genomic coordinates (regions and specific positions)
6. Query features overlapping genomic regions with filtering by type and strand
7. Identify features containing specific genomic positions

## GFF Analysis - Relationship and Hierarchy Queries:
8. Explore gene structure and organization (get all child features like exons, CDS, UTRs)
9. Find parent features of any given feature using upward traversal
10. Get all features of specific types with efficient iteration

## GFF Analysis - Statistical Analysis:
11. Calculate comprehensive feature statistics (counts, length distributions per feature type)
12. Generate per-chromosome feature summaries and analysis
13. Analyze length distributions with detailed statistics (min, max, mean, median, std dev, percentiles)
14. Create histogram data for feature length distributions

## GFF Analysis - Attribute-based Searches:
15. Search features by attribute key-value pairs (exact and partial matching)
16. Find features containing specific attribute keys
17. Support pattern matching and logical operations for attribute queries

## GFF Analysis - Positional Analysis:
18. Identify intergenic regions (gaps between genes) with filtering options
19. Calculate feature density in genomic windows across chromosomes
20. Analyze strand distribution of features with counts and percentages
21. Support clustering analysis and positional insights

## GFF Analysis - Export and Reporting:
22. Export feature data to CSV format with comprehensive filtering
23. Generate human-readable summary reports of GFF file contents
24. Provide formatted output for downstream analysis

Example queries you can handle:
- "Find all genes in chromosome 1 between positions 1000-5000"
- "Get the structure of gene AT1G01010 including all exons and CDS"
- "Calculate feature statistics for this GFF file"
- "Find all features with 'kinase' in their Name attribute"
- "Identify intergenic regions longer than 1000bp on chromosome 2"
- "Calculate gene density in 10kb windows across chromosome 1"
- "Export all exon features to CSV format"

When using tools:
- Always verify file paths before operations
- Provide clear explanations of what you're doing
- If a task cannot be completed, explain why and suggest alternatives
- Use appropriate filtering and limiting for large datasets
- Handle database creation and reuse efficiently

Always be helpful, concise, and focus on addressing the user's bioinformatics analysis needs efficiently.
"""

    model_id = (
    "gpt-oss:20b-cloud" 
)

    #host="https://ollama.com",
    #headers={'Authorization': 'Bearer ' + os.environ.get('OLLAMA_API_KEY')}



    ollama_model = OllamaModel(
    model_id=model_id,
    #host="http://localhost:11434",
    host="https://ollama.com",
    headers={'Authorization': 'Bearer ' + os.environ.get('OLLAMA_API_KEY')},
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
           get_multiple_gene_lenght, get_gene_attributes, get_features_in_region, get_features_at_position,
           get_gene_structure, get_feature_parents, get_features_by_type, get_feature_statistics,
           get_chromosome_summary, get_length_distribution, search_features_by_attribute, 
           get_features_with_attribute, get_intergenic_regions, get_feature_density, get_strand_distribution,
           export_features_to_csv, get_feature_summary_report],
)


    print("Hello from gffutilsai!")
    #r=local_agent(
    #"want to know which features there are in a gff file called ./subset_4percent.gff. Show me the list"
    #)
    #r=local_agent(
    #"want to know which chromosomes there are in the gff file called ./subset_4percent.gff. Show me the list"
    #)
    r=local_agent(
    "want to know strand distribution for genes in the gff file called ./subset_4percent.gff."
    )

if __name__ == "__main__":
    main()
