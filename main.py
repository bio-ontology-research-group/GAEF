import csv
import json
import os
import argparse
import subprocess
import plotly.express as px
from flask import Flask, render_template
from collections import defaultdict
from information_content import calculate_ic_depth_breadth
from completeness import read_essential_terms, count_essential_terms
from utils import get_ancestors, get_specific, Ontology
from consistency import check_consistency
import plots
import coherence

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

term_file = "constraints/essential_terms.tsv"
has_part_file = "constraints/has_part_relations.txt"
ec2go_file = "constraints/ec2go_v2025-03-16"
pathway_file = "constraints/metacyc_GO_v2025-03-16_with_EC.tsv"
ontology_file = "data/go-basic.obo"
taxa_constraints_file = "constraints/taxon_constraints.tsv"
MACROMOLECULAR_COMPLEX = "GO:0032991"
HOMODIMERIZATION = "GO:0042803"

# @app.route('/')
def evaluation(assembly_name, annotation_file, groovy_flag=False):
	# Load annotations
	protein_go_terms = load_annotations(annotation_file)
	protein_go_terms_ancestors = get_ancestors(protein_go_terms)
	protein_go_terms_specific = get_specific(protein_go_terms)
	go = Ontology(ontology_file)
	specific_terms_file = f"{assembly_name}_specific_GO_terms.tsv"

	with open(specific_terms_file, "w") as f:
		for protein, terms in protein_go_terms_specific.items():
			if terms:  # skip proteins with no terms
				f.write(protein + "\t" + "\t".join(sorted(terms)) + "\n")


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
	ec2go_mapping      = coherence.parse_ec2go(ec2go_file)
	pathway_to_go      = coherence.map_pathways_to_go_terms(pathway_file, ec2go_mapping)
	_, metacyc_completed, metacyc_annotated, pathway_details = coherence.analyze_genome(protein_go_terms_ancestors, pathway_to_go, ec2go_mapping)
	metacyc_pct       = (len(metacyc_completed) / len(metacyc_annotated)) * 100
	total_completed   = len(metacyc_completed)
	total_annotated   = len(metacyc_annotated)
	total_incomplete  = total_annotated - total_completed

	# Process coherence
	has_part_dict = coherence.parse_has_part(has_part_file)
	process_coherence, has_part_protein_details = coherence.check_has_part(protein_go_terms, has_part_dict)

	# Protein complex coherence
	term_to_children, _, _ = coherence.parse_go_ontology(ontology_file)
	complex_child_terms = coherence.get_all_child_terms(MACROMOLECULAR_COMPLEX, term_to_children)
	complex_child_terms.add(MACROMOLECULAR_COMPLEX)
	complex_classifications, _ = coherence.classify_complexes(protein_go_terms_ancestors, complex_child_terms)
	coherent_count, incoherent_count = coherence.count_complexes(complex_classifications)
	complex_coherence = (coherent_count / (coherent_count + incoherent_count)) * 100 if (coherent_count + incoherent_count) > 0 else 0
	term_names = {t: go.get_term(t)['name'] for t in complex_classifications}

	### CONSISTENCY ###
	# Taxonomic consistency
	consistency_file =  check_consistency(protein_go_terms_ancestors, taxa_constraints_file, output_file = assembly_name + "_consistency.tsv")
	if groovy_flag:
		print("Evaluating taxonomic consistency with Groovy.")
		groovy_script = "groovy_scripts/taxon_consistency.groovy"
		groovy_output_file = f"{assembly_name}_taxon_explanations.tsv"
		try:
			result = subprocess.run(["groovy", groovy_script, consistency_file, groovy_output_file], capture_output=True, text=True, check=True)
			print("Groovy script output:")
			print(result.stdout)
			if result.stderr:
				print("Groovy script error output:")
				print(result.stderr)
		except subprocess.CalledProcessError as e:
			print("Error running Groovy script:")
			print(e.stderr)
            
	### OVERVIEW ###
	completeness_data = [
		essential_percentage,
		metacyc_pct,
		process_coherence, 
		complex_coherence]
     
	gauge_html = plots.create_completeness_gauge_html(completeness_data)
	
	context = {
        'assembly_name': assembly_name,
        'essential_percentage': round(essential_percentage, 2),
        'plot_core_html': plot_core_html,
        'plot_periph_html': plot_periph_html,
        'go_categories_core': go_core,
        'go_categories_periph': go_periph,
        'found_terms': found_terms,
        'ec2go_mapping': ec2go_mapping,
        'complete_has_part_percentage': round(process_coherence, 2),
        'has_part_data': has_part_protein_details,
        'complex_coherence': round(complex_coherence, 2),
        'complex_classifications': complex_classifications,
        'term_names': term_names,
        'metacyc_complete_percentage': round(metacyc_pct, 2),
        'metacyc_completed': total_completed,
        'metacyc_annotated': total_annotated,
        'metacyc_incomplete': total_incomplete,
        'pathway_details': pathway_details,
        'gauge_html': gauge_html
    }
   
	if groovy_flag:
		satisfiable = True
		if os.path.exists(groovy_output_file):
			with open(groovy_output_file, encoding="utf-8") as f:
				next(f)  # skip header
				for line in f:
					parts = line.strip().split("\t")
					if len(parts) >= 2 and parts[1].lower() == "false":
						satisfiable = False
						break
			context["satisfiable"] = satisfiable

	if groovy_flag:
		print("Calculating Information Content (IC) with Groovy.")
		groovy_IC = "groovy_scripts/ICVector.groovy"
		groovy_IC_output_file = f"{assembly_name}_IC.tsv"
		try:
			result = subprocess.run(["groovy", groovy_IC, specific_terms_file, groovy_IC_output_file], capture_output=True, text=True, check=True)
			print("Groovy script output:")
			print(result.stdout)
			if result.stderr:
				print("Groovy script error output:")
				print(result.stderr)
		except subprocess.CalledProcessError as e:
			print("Error running Groovy script:")
			print(e.stderr)

	if groovy_flag:
		# Calculate IC depth and breadth
		depth, breadth, normalized_breadth = calculate_ic_depth_breadth(groovy_IC_output_file)
		context["ic_depth"] = depth
		context["ic_breadth"] = breadth
		context["normalized_ic_breadth"] = normalized_breadth

	return context

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate GO annotation completeness/coherence/consistency.')
    parser.add_argument('--assembly_name', required=True, help='Assembly name (e.g., GCF_000007085.1_ASM708v1)')
    parser.add_argument('--annotation_file', required=True, help='Path to tab-separated GO annotation file (protein_id	GO:term1	GO:term2)')
    parser.add_argument('--groovy', action='store_true', help='Run Groovy scripts for IC calculation and taxonomic consistency')

    args = parser.parse_args()
    
    with app.app_context():
        context = evaluation(args.assembly_name, args.annotation_file, args.groovy)

        # Save HTML using full context
        html = render_template("html_output_template.html", **context)
        with open(args.assembly_name + "_report.html", "w", encoding="utf-8") as f:
            f.write(html)

        # Save JSON with selected fields removed
        json_context = context.copy()
        json_context.pop("plot_core_html", None)
        json_context.pop("plot_periph_html", None)
        json_context.pop("gauge_html", None)
        json_context.pop("ec2go_mapping", None)
        json_context.pop("term_names", None)

        with open(args.assembly_name + "_report.json", "w") as f:
            json.dump(json_context, f, indent=2, default=str)

        print(f"Saved {args.assembly_name}_report.html and {args.assembly_name}_report.json")
        
