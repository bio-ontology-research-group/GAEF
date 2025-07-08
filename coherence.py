import ast
import re
from collections import defaultdict

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


def check_has_part(protein_go_terms, has_part_dict):
    """
    Check if the required 'has-part' GO terms are included in the genome's GO terms.
    Returns the percentage of missing 'has-part' relations at the genome level.
    """
    genome_go_terms = set()
    for go_terms in protein_go_terms.values():
        genome_go_terms.update(go_terms)

    missing_relations_count = 0
    total_relations_count = 0
    protein_details = []
    for protein_id, go_terms in protein_go_terms.items():
        for go_term in genome_go_terms:
            if go_term in has_part_dict:
                required_parts = has_part_dict[go_term]
                missing = sorted(required_parts - genome_go_terms)
                missing_relations_count += len(missing)
                total_relations_count += len(required_parts)
                protein_details.append({
						  "protein_id": protein_id,
						  "annotated_term": go_term,
						  "missing_parts": missing
					 })

    if total_relations_count == 0:
        return 0
    missing_percentage = (missing_relations_count / total_relations_count) * 100
    process_coherence = (100 - missing_percentage)
    return process_coherence, protein_details

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


def analyze_genome(protein_go_terms, pathway_to_go, ec2go_mapping):
    completeness_results = {}
    completed_pathways = []
    annotated_pathways = set()
    pathway_details = {}

    genome_go_set = set(go_term for go_terms in protein_go_terms.values() for go_term in go_terms)

    for pathway, (original_go_term, go_terms_sets) in pathway_to_go.items():
        # Only proceed if the original GO term is present
        if original_go_term not in genome_go_set:
            continue

        annotated_pathways.add(pathway)
        missing_components = []

        if not go_terms_sets:
            pathway_complete = True
        else:
            pathway_complete = any(
                all(go_term in genome_go_set for go_term in combo)
                for combo in go_terms_sets
            )
            if not pathway_complete:
                for combo in go_terms_sets:
                    missing_ec_go_terms = [go_term for go_term in combo if go_term not in genome_go_set]
                    if missing_ec_go_terms:
                        missing_components.append(missing_ec_go_terms)

        completeness_results[pathway] = pathway_complete
        if pathway_complete:
            completed_pathways.append(pathway)
            
        if ec2go_mapping and missing_components:
            mapped = []
            for go_set in missing_components:
                ecs = set()
                for go in go_set:
                    for ec, gos in ec2go_mapping.items():
                        if go in gos:
                            ecs.add(ec)
                if ecs:
                    mapped.append(list(ecs))
            missing_components = mapped

        pathway_details[pathway] = {
            "complete": pathway_complete,
            "original_go_term": original_go_term,
            "missing_components": missing_components
        }

    return completeness_results, completed_pathways, annotated_pathways, pathway_details

#### COMPLEX COHERENCE ####
MACROMOLECULAR_COMPLEX = "GO:0032991"
HOMODIMERIZATION = "GO:0042803"

def classify_complexes(protein_go_terms, complex_child_terms, homodimer_terms=None):
    complex_term_to_proteins = defaultdict(set)
    for protein_id, terms in protein_go_terms.items():
        if MACROMOLECULAR_COMPLEX in terms:
            protein_complex_terms = complex_child_terms.intersection(terms)
            for term in protein_complex_terms:
                complex_term_to_proteins[term].add(protein_id)

    complex_classifications = {}
    for term, proteins in complex_term_to_proteins.items():
        if homodimer_terms and term in homodimer_terms:
            complex_classifications[term] = "coherent"
        elif len(proteins) > 1:
            complex_classifications[term] = "coherent"
        else:
            complex_classifications[term] = "incoherent"

    return complex_classifications, complex_term_to_proteins


def count_complexes(complex_classifications):
    coherent_count = sum(1 for c in complex_classifications.values() if c == "coherent")
    incoherent_count = sum(1 for c in complex_classifications.values() if c == "incoherent")
    return coherent_count, incoherent_count


def parse_go_ontology(file_path):
    term_to_children = defaultdict(set)
    term_to_parents = defaultdict(set)
    term_to_definition = {}

    current_term = None
    term_id = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '[Term]':
                current_term = True
                term_id = None
                continue
            if line == '' and current_term:
                current_term = False
                continue
            if not current_term:
                continue
            if line.startswith('id:'):
                term_id = line.split('id:')[1].strip()
                continue
            if line.startswith('def:'):
                if term_id:
                    definition = line.split('def:')[1].strip()
                    match = re.search(r'"([^"]*)"', definition)
                    term_to_definition[term_id] = match.group(1) if match else definition
                continue
            if line.startswith('is_a:'):
                if term_id:
                    parent_id = line.split('is_a:')[1].split('!')[0].strip()
                    term_to_parents[term_id].add(parent_id)
                    term_to_children[parent_id].add(term_id)
                continue
            if line.startswith('relationship: part_of'):
                if term_id:
                    parent_id = line.split('part_of')[1].split('!')[0].strip()
                    term_to_parents[term_id].add(parent_id)
                    term_to_children[parent_id].add(term_id)

    return term_to_children, term_to_parents, term_to_definition


def get_all_child_terms(term_id, term_to_children, recursive=True):
    if not recursive:
        return term_to_children.get(term_id, set())

    all_children = set()
    to_process = list(term_to_children.get(term_id, set()))

    while to_process:
        child = to_process.pop()
        if child not in all_children:
            all_children.add(child)
            to_process.extend(term_to_children.get(child, set()))

    return all_children
