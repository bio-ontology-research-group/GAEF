import csv
from information_content import calculate_ic_depth_breadth
from completeness import read_essential_terms, count_essential_terms
from utils import get_ancestors, get_specific
from consistency import check_consistency
from coherence import (
    parse_has_part,
    check_completeness,
    parse_ec2go,
    map_pathways_to_go_terms,
    analyze_genome,
    parse_go_ontology,
    get_all_child_terms,
    classify_complexes,
    count_complexes,
    MACROMOLECULAR_COMPLEX
)

def load_annotations(annotation_file):
    protein_go_terms = {}
    with open(annotation_file, 'r') as f:
        tsv_reader = csv.reader(f, delimiter='\t')
        for row in tsv_reader:
            if row:
                protein_id = row[0]
                go_terms = set(row[1:])  # Skip the protein ID
                protein_go_terms[protein_id] = go_terms
    return protein_go_terms

if __name__ == "__main__":
    annotation_file = "GCF_000007085.1_ASM708v1_IPscan_GO.tsv"
    term_file = "constraints/essential_terms.tsv"
    has_part_file = "constraints/has_part_relations.txt"
    ec2go_file = "constraints/ec2go_v2025-03-16"
    pathway_file = "constraints/metacyc_GO_v2025-03-16_with_EC.tsv"
    ontology_file = "data/go-basic.obo"
    taxa_constraints_file = "constraints/taxon_constraints.tsv"

    # Load annotations
    protein_go_terms = load_annotations(annotation_file)
    protein_go_terms_ancestors = get_ancestors(protein_go_terms)
    protein_go_terms_specific = get_specific(protein_go_terms)
    
    # Essential terms
    core_terms, peripheral_terms = read_essential_terms(term_file)
    core_presence = count_essential_terms(protein_go_terms, core_terms)
    peripheral_presence = count_essential_terms(protein_go_terms, peripheral_terms)

    # Coherence
    has_part_dict = parse_has_part(has_part_file)
    process_coherence = check_completeness(protein_go_terms_ancestors, has_part_dict)

    # Pathway coherence
    ec2go_mapping = parse_ec2go(ec2go_file)
    pathway_to_go = map_pathways_to_go_terms(pathway_file, ec2go_mapping)
    completeness_results, completed_pathways, annotated_pathways = analyze_genome(protein_go_terms_ancestors, pathway_to_go)

    total_complete_pathways = sum(completeness_results.values())
    annotated_pathways_count = len(annotated_pathways)
    proportion_completed = (total_complete_pathways / annotated_pathways_count) * 100 if annotated_pathways_count > 0 else 0

    # Complex coherence
    term_to_children, _, _ = parse_go_ontology(ontology_file)
    complex_child_terms = get_all_child_terms(MACROMOLECULAR_COMPLEX, term_to_children)
    complex_child_terms.add(MACROMOLECULAR_COMPLEX)

    complex_classifications, _ = classify_complexes(protein_go_terms_ancestors, complex_child_terms)
    coherent_count, incoherent_count = count_complexes(complex_classifications)

    check_consistency(protein_go_terms, taxa_constraints_file, output_file="GCF_000005845.2_ASM584v2_consistency.tsv")

    # Output
    print("Core GO term presence:")
    for term, present in core_presence.items():
        print(f"{term}\t{present}")

    print("\nPeripheral GO term presence:")
    for term, present in peripheral_presence.items():
        print(f"{term}\t{present}")

    print(f"\nProcess Coherence: {process_coherence:.2f}%")

    print("\nPathway Coherence:")
    print(f"Total completed pathways: {total_complete_pathways}")
    print(f"Annotated pathways: {annotated_pathways_count}")
    print(f"Proportion completed: {proportion_completed:.2f}%")

    print("\nComplex Coherence:")
    print(f"Coherent complexes: {coherent_count}")
    print(f"Incoherent complexes: {incoherent_count}")
    
    ic_stats = calculate_ic_depth_breadth("GCF_000007085.1_ASM708v1_IPscan_IC_nn")
    print(ic_stats)
