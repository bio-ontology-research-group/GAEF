import csv
from essential_terms import read_essential_terms, count_essential_terms
from coherence import (
    parse_has_part,
    check_completeness,
    parse_ec2go,
    map_pathways_to_go_terms,
    analyze_genome
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
    annotation_file = "GCF_000005845.2_ASM584v2_IPscan_GO_ancestors.tsv"
    term_file = "constraints/essential_terms.tsv"
    has_part_file = "constraints/has_part_relations.txt"
    ec2go_file = "constraints/ec2go_v2025-03-16"
    pathway_file = "constraints/metacyc_GO_v2025-03-16_with_EC.tsv"

    # Load annotations
    protein_go_terms = load_annotations(annotation_file)

    # Essential terms
    core_terms, peripheral_terms = read_essential_terms(term_file)
    core_presence = count_essential_terms(protein_go_terms, core_terms)
    peripheral_presence = count_essential_terms(protein_go_terms, peripheral_terms)

    # Coherence: process
    has_part_dict = parse_has_part(has_part_file)
    process_coherence = check_completeness(protein_go_terms, has_part_dict)

    # Coherence: pathways
    ec2go_mapping = parse_ec2go(ec2go_file)
    pathway_to_go = map_pathways_to_go_terms(pathway_file, ec2go_mapping)
    completeness_results, completed_pathways, annotated_pathways = analyze_genome(protein_go_terms, pathway_to_go)

    total_complete_pathways = sum(completeness_results.values())
    annotated_pathways_count = len(annotated_pathways)
    proportion_completed = (total_complete_pathways / annotated_pathways_count) * 100 if annotated_pathways_count > 0 else 0


    # Output results
    print("Core GO term presence:")
    for term, present in core_presence.items():
        print(f"{term}\t{present}")

    print("\nPeripheral GO term presence:")
    for term, present in peripheral_presence.items():
        print(f"{term}\t{present}")

    print(f"\nProcess Coherence: {process_coherence:.2f}%")

    print(f"\nPathway Completeness:")
    print(f"Total completed pathways: {total_complete_pathways}")
    print(f"Annotated pathways: {annotated_pathways_count}")
    print(f"Proportion completed: {proportion_completed:.2f}%")
