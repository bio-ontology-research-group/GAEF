import csv
import plotly.express as px
from flask import Flask, render_template
from collections import defaultdict
from information_content import calculate_ic_depth_breadth
from completeness import read_essential_terms, count_essential_terms
from utils import get_ancestors, get_specific, Ontology
from consistency import check_consistency
from coherence import (
    parse_has_part,
    check_has_part,
    parse_ec2go,
    map_pathways_to_go_terms,
    analyze_genome,
    parse_go_ontology,
    get_all_child_terms,
    classify_complexes,
    count_complexes
)

app = Flask(__name__, template_folder='templates', static_folder='static')

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

annotation_file = "GCF_000007085.1_ASM708v1_IPscan_GO.tsv"
term_file = "constraints/essential_terms.tsv"
has_part_file = "constraints/has_part_relations.txt"
ec2go_file = "constraints/ec2go_v2025-03-16"
pathway_file = "constraints/metacyc_GO_v2025-03-16_with_EC.tsv"
ontology_file = "data/go-basic.obo"
taxa_constraints_file = "constraints/taxon_constraints.tsv"
MACROMOLECULAR_COMPLEX = "GO:0032991"
HOMODIMERIZATION = "GO:0042803"

@app.route('/')
def evaluation():
	# Load annotations
	protein_go_terms = load_annotations(annotation_file)
	protein_go_terms_ancestors = get_ancestors(protein_go_terms)
	protein_go_terms_specific = get_specific(protein_go_terms)
	go = Ontology(ontology_file)

	### COMPLETENESS ###
	# Essential terms
	core_terms, periph_terms = read_essential_terms(term_file)
	core_ids = {e['term'] for e in core_terms}
	periph_ids = {e['term'] for e in periph_terms}
	core_labels = [e['Function'] for e in core_terms]
	periph_labels = [e['Function'] for e in periph_terms]

	core_count = count_essential_terms(protein_go_terms_ancestors, core_ids)
	periph_count = count_essential_terms(protein_go_terms_ancestors, periph_ids)
	essential_percentage = (sum(core_count.values())/len(core_ids)) * 100 if core_ids else 0	

	# Completeness plots
	core_vals   = [core_count.get(e['term'], 0) for e in core_terms]
	fig_core    = px.bar(x=core_labels, y=core_vals,
		labels={'x': 'Core Function', 'y': 'Presence'})
	plot_core_html = fig_core.to_html(full_html=False, include_plotlyjs='cdn')

	periph_vals   = [periph_count.get(e['term'], 0) for e in periph_terms]
	fig_periph    = px.bar(x=periph_labels, y=periph_vals,
		labels={'x': 'Peripheral Function', 'y': 'Presence'})
	plot_periph_html = fig_periph.to_html(full_html=False, include_plotlyjs=False)

	# Tables grouping
	go_core   = {'Core': core_terms}
	go_periph = defaultdict(list)
	for e in periph_terms:
		go_periph[e['category']].append(e)

	# Found terms set
	found_terms = {term for term, pres in core_count.items() if pres} | {term for term, pres in periph_count.items() if pres}

	### COHERENCE ###
   # Pathway coherence
	ec2go_mapping      = parse_ec2go(ec2go_file)
	pathway_to_go      = map_pathways_to_go_terms(pathway_file, ec2go_mapping)
	_, metacyc_completed, metacyc_annotated, pathway_details = analyze_genome(protein_go_terms_ancestors, pathway_to_go, ec2go_mapping)
	metacyc_pct       = (len(metacyc_completed) / len(metacyc_annotated)) * 100
	total_completed   = len(metacyc_completed)
	total_annotated   = len(metacyc_annotated)
	total_incomplete  = total_annotated - total_completed

	# Process coherence
	has_part_dict = parse_has_part(has_part_file)
	process_coherence, has_part_protein_details = check_has_part(protein_go_terms, has_part_dict)

	# Protein complex coherence
	term_to_children, _, _ = parse_go_ontology(ontology_file)
	complex_child_terms = get_all_child_terms(MACROMOLECULAR_COMPLEX, term_to_children)
	complex_child_terms.add(MACROMOLECULAR_COMPLEX)
	complex_classifications, _ = classify_complexes(protein_go_terms_ancestors, complex_child_terms)
	coherent_count, incoherent_count = count_complexes(complex_classifications)
	complex_coherence = (coherent_count / (coherent_count + incoherent_count)) * 100 if (coherent_count + incoherent_count) > 0 else 0
	term_names = {t: go.get_term(t)['name'] for t in complex_classifications}

	### CONSISTENCY ###
	# Taxonomic consistency
	check_consistency(protein_go_terms, taxa_constraints_file, output_file="GCF_000005845.2_ASM584v2_consistency.tsv")

	return render_template(
        'html_output_template.html',
        essential_percentage=round(essential_percentage, 2),
        plot_core_html=plot_core_html,
        plot_periph_html=plot_periph_html,
        go_categories_core=go_core,
        go_categories_periph=go_periph,
        found_terms=found_terms,
        ec2go_mapping=ec2go_mapping,
        complete_has_part_percentage=round(process_coherence, 2),
        has_part_data=has_part_protein_details,
        complex_coherence=round(complex_coherence, 2),
        complex_classifications=complex_classifications,
        term_names=term_names,
        satisfiable=True,
        metacyc_complete_percentage=round(metacyc_pct, 2),
        metacyc_completed=total_completed,
        metacyc_annotated=total_annotated,
        metacyc_incomplete=total_incomplete,
        pathway_details=pathway_details
    )

if __name__ == '__main__':
    app.run(debug=True)