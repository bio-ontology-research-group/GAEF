import os
import gzip
from collections import defaultdict

def extract_go_terms(input_dir="gaf_downloads", output_dir="tsv_outputs"):
    """
    Reads .gaf or .gaf.gz files, extracts Protein Accessions and GO IDs,
    and writes them to .tsv files grouped by protein.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all GAF files in the input directory
    if not os.path.exists(input_dir):
        print(f"Error: The directory '{input_dir}' does not exist.")
        return

    files_to_process = [f for f in os.listdir(input_dir) if f.endswith('.gaf') or f.endswith('.gaf.gz')]
    
    if not files_to_process:
        print(f"No .gaf or .gaf.gz files found in '{input_dir}'.")
        return

    print(f"Found {len(files_to_process)} files. Starting extraction...\n")

    for filename in files_to_process:
        input_path = os.path.join(input_dir, filename)
        
        # Change extension to .tsv
        out_name = filename.replace('.gaf.gz', '.tsv').replace('.gaf', '.tsv')
        output_path = os.path.join(output_dir, out_name)

        # Dictionary to hold protein accession as key, and a set of GO IDs as values
        # Use a set to automatically remove any duplicate GO terms for the same protein
        protein_to_go = defaultdict(set)

        open_func = gzip.open if filename.endswith('.gz') else open
        
        print(f"Processing: {filename}")
        
        try:
            with open_func(input_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    # Skip metadata/comment lines
                    if line.startswith('!'):
                        continue
                    
                    # Split the tab-separated line
                    columns = line.strip('\n').split('\t')
                    
                    # Ensure the line has enough columns to prevent index errors
                    if len(columns) > 4:
                        protein_accession = columns[1]
                        go_id = columns[4]
                        
                        # Basic validation to ensure it's a GO ID
                        if go_id.startswith('GO:'):
                            protein_to_go[protein_accession].add(go_id)
            
            # Write the grouped data to a TSV file
            with open(output_path, 'w', encoding='utf-8') as out_f:
                for protein, go_ids in protein_to_go.items():
                    sorted_go_ids = sorted(list(go_ids))
                    
                    # Join the protein and its GO IDs with tabs
                    row_data = [protein] + sorted_go_ids
                    out_f.write('\t'.join(row_data) + '\n')
                    
            print(f"  -> Saved to {out_name}")

        except Exception as e:
            print(f"  -> Error processing {filename}: {e}")

if __name__ == "__main__":
    INPUT_FOLDER = "gaf_downloads"
    OUTPUT_FOLDER = "tsv_outputs"
    
    extract_go_terms(INPUT_FOLDER, OUTPUT_FOLDER)
