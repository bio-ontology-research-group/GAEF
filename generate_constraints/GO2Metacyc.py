import sys

go_version = sys.argv[1]
go_obo_file = sys.argv[2]

def parse_obo_file(obo_file, output_file):
    with open(obo_file, 'r') as file:
        lines = file.readlines()

    go_term = None
    meta_cyc = []
    data = []

    for line in lines:
        if line.startswith('[Term]'):
            if go_term and meta_cyc:
                for pathway in meta_cyc:
                    data.append(f"{go_term}\t{pathway}\n")
            go_term = None
            meta_cyc = []
        elif line.startswith('id: '):
            go_term = line.strip().split(' ')[1]
        elif line.startswith('xref: MetaCyc:') and 'PWY' in line:
            pathway = line.strip().split(' ')[1].split(':')[1]
            meta_cyc.append(pathway)

    # to capture the last term in file
    if go_term and meta_cyc:
        for pathway in meta_cyc:
            data.append(f"{go_term}\t{pathway}\n")

    with open(output_file, 'w') as out_file:
        out_file.write('GO_Term\tMetaCyc_Pathway\n')
        out_file.writelines(data)

# Use the function with your specific file paths
parse_obo_file(go_obo_file, f'metacyc_GO_v{go_version}.tsv')

