import pandas as pd
from collections import defaultdict

def check_consistency(protein_go_terms, constraints_file, output_file):
    """
    Check taxonomic consistency of GO annotations and save results to a file.

    Parameters:
    - protein_go_terms (dict): Mapping from protein ID to set of GO terms
    - constraints_file (str): Path to the constraints file
    - output_file (str): Path to write the output table
    """

    # === Process constraints ===
    constraints = pd.read_csv(constraints_file, sep="\t", dtype=str)
    only_map = defaultdict(set)
    never_map = defaultdict(set)

    for _, row in constraints.iterrows():
        go_id = row["GO_ID"]
        lineage = row["Taxon_ID"]
        if row["Constraint_Type"] == "only_in_taxon":
            only_map[go_id].add(lineage)
        elif row["Constraint_Type"] == "never_in_taxon":
            never_map[go_id].add(lineage)

    # === Write output ===
    with open(output_file, 'w') as out:
        out.write("protein_name\tGO_ID\tnever_in_taxon\tonly_in_taxon\n")
        for protein_id, go_terms in protein_go_terms.items():
            for go in go_terms:
                only_taxa = list(only_map[go])
                never_taxa = list(never_map[go])
                if only_taxa or never_taxa:
                    row = [
                        protein_id,
                        go,
                        ",".join(never_taxa),
                        ",".join(only_taxa)
                    ]
                    out.write("\t".join(row) + "\n")
