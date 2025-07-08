import sys

# Read essential terms from file
def read_essential_terms(term_file):
    core_entries = []
    periph_entries = []
    with open(term_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            go_term, category, function = parts
            entry = {'term': go_term, 'category': category, 'Function': function}
            if category == 'Core':
                core_entries.append(entry)
            else:
                periph_entries.append(entry)
    return core_entries, periph_entries

# Count presence/absence of GO terms
def count_essential_terms(protein_go_terms, target_go_terms):
    go_term_presence = {term: 0 for term in target_go_terms}

    for go_terms in protein_go_terms.values():
        for go_term in go_terms:
            if go_term in target_go_terms:
                go_term_presence[go_term] = 1

    return go_term_presence