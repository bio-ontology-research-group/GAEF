import ast
#### PROCESS COHERENCE ####
def parse_has_part(file_path):
    """
    Parse a file where each line contains two space-separated GO terms:
       GO:XXXXXXX GO:YYYYYYY
    Returns a dictionary mapping each GO term to a set of GO terms it has-part.
    """
    has_part_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            parts = line.split()
            if len(parts) == 2:
                key = parts[1]
                value = parts[0]
                if key in has_part_dict:
                    has_part_dict[key].add(value)
                else:
                    has_part_dict[key] = {value}
    return has_part_dict


def check_completeness(protein_go_terms, has_part_dict):
    """
    Check if the required 'has-part' GO terms are included in the genome's GO terms.
    Returns the percentage of missing 'has-part' relations at the genome level.
    """
    genome_go_terms = set()
    for go_terms in protein_go_terms.values():
        genome_go_terms.update(go_terms)

    missing_relations_count = 0
    total_relations_count = 0

    for go_term in genome_go_terms:
        if go_term in has_part_dict:
            required_parts = has_part_dict[go_term]
            missing = required_parts - genome_go_terms
            missing_relations_count += len(missing)
            total_relations_count += len(required_parts)

    if total_relations_count == 0:
        return 0
    missing_percentage = (missing_relations_count / total_relations_count) * 100
    return (100 - missing_percentage)

#### PATHWAY COHERENCE ####
def parse_ec2go(filename):
    ec2go = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('EC:'):
                parts = line.strip().split(' > ')
                ec_number = parts[0].split(':')[1]
                go_term = parts[1].split('; ')[1]
                if ec_number in ec2go:
                    ec2go[ec_number].append(go_term)
                else:
                    ec2go[ec_number] = [go_term]
    return ec2go


def map_pathways_to_go_terms(pathway_file, ec2go):
    pathway_to_go = {}
    with open(pathway_file, 'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split('\t')
            original_go_term, pathway, ec_combinations = parts[0], parts[1], parts[2]
            if ec_combinations.strip() == '':
                go_terms_sets = []
            else:
                ec_combination_lists = ast.literal_eval(ec_combinations)
                go_terms_sets = []
                for ec_list in ec_combination_lists:
                    go_terms_set = set()
                    for ec in ec_list:
                        if ec in ec2go:
                            go_terms_set.update(ec2go[ec])
                    if go_terms_set:
                        go_terms_sets.append(go_terms_set)
            pathway_to_go[pathway] = (original_go_term, go_terms_sets)
    return pathway_to_go


def analyze_genome(protein_go_terms, pathway_to_go):
    completeness_results = {}
    completed_pathways = []
    annotated_pathways = set()
    
    genome_go_set = set()
    for go_terms in protein_go_terms.values():
        genome_go_set.update(go_terms)

    for pathway, (original_go_term, go_terms_sets) in pathway_to_go.items():
        if original_go_term in genome_go_set:
            annotated_pathways.add(pathway)
            if not go_terms_sets:
                pathway_complete = True
            else:
                pathway_complete = any(all(go_term in genome_go_set for go_term in go_terms) for go_terms in go_terms_sets)
        else:
            pathway_complete = False

        completeness_results[pathway] = pathway_complete
        if pathway_complete:
            completed_pathways.append(pathway)

    return completeness_results, completed_pathways, annotated_pathways

#### PATHWAY COHERENCE ####
