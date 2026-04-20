import os
import urllib.request
import urllib.error

def download_gaf_files(input_file, output_dir="gaf_downloads"):
    """
    Reads genome assembly names from a text file and downloads their GAF annotations.
    """
    # Create the output directory if it doesn't already exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the text file, ignoring empty lines
    with open(input_file, 'r') as f:
        genomes = [line.strip() for line in f if line.strip()]

    print(f"Found {len(genomes)} genomes in {input_file}. Starting downloads...\n")

    for genome in genomes:
        # Example genome: GCF_000007085.1_ASM708v1
        try:
            # Extract the necessary parts to navigate NCBI's FTP folder structure
            prefix = genome[0:3]    # e.g., 'GCF' or 'GCA'
            digits = genome[4:13]   # e.g., '000007085'
            
            # NCBI splits the 9-digit ID into chunks of 3 for their folder tree
            part1 = digits[0:3]
            part2 = digits[3:6]
            part3 = digits[6:9]

            # Construct the exact URL
            base_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{part1}/{part2}/{part3}/{genome}"
            filename = f"{genome}_gene_ontology.gaf.gz"
            download_url = f"{base_url}/{filename}"

            output_path = os.path.join(output_dir, filename)

            print(f"Fetching: {filename}")
            
            # Download the file
            urllib.request.urlretrieve(download_url, output_path)
            print("  -> Success!")

        except urllib.error.HTTPError as e:
            # A 404 error usually means NCBI didn't generate a GAF file for this specific assembly or they suppressed it
            print(f"  -> Failed: HTTP Error {e.code}. (The GAF file likely doesn't exist for this genome)")
        except urllib.error.URLError as e:
            print(f"  -> Failed: Connection Error - {e.reason}")
        except Exception as e:
            print(f"  -> Error processing {genome}: {e}")

if __name__ == "__main__":
    INPUT_TXT = "../evaluation_dataset.tsv" 
    OUTPUT_FOLDER = "gaf_downloads"
    
    download_gaf_files(INPUT_TXT, OUTPUT_FOLDER)
