# Genome Annotation Evaluation Framework

This tool generates a static HTML and JSON report evaluating Gene Ontology (GO) class annotations using the following metrics: **completeness**, **coherence**, and **consistency**. It integrates multiple constraint sources (e.g., EC2GO mappings, GO structure, NCBITaxonomy).

[![DOI](https://zenodo.org/badge/1014841170.svg)](https://doi.org/10.5281/zenodo.15902372)


---

## Evaluation Metrics

- **Completeness**: presence of essential functions for bacterial genomes
- **Coherence**:
  - *Process coherence*: Evaluates presence of dependencies of annotated classes using `has_part` relations
  - *Pathway coherence*: Checks presence of MetaCyc pathway components using EC numbers and their GO class mappings
  - *Complex coherence*: Identifies multi-protein complexes and checks the presence of their components
- **Consistency**: Uses GO taxon constraints and NCBI taxonomy to ensure annotations are taxonomically consistent

---

## Requirements

- Python 3.8+
- Groovy (Optional; required for taxonomic consistency evaluation and Information Content (IC) calculation)
- `Flask`
- `plotly`
- `pandas`
- `numpy`


For generating constriants:
- `pythoncyc`

---

## Setup

1. Clone the repository:

```https://github.com/bio-ontology-research-group/GAEF.git```  
```cd GAEF```

2. Download Gene Ontology (GO) files:

```mkdir data```  
```cd data```  
```wget https://release.geneontology.org/2025-03-16/ontology/go-basic.obo```  
```wget https://release.geneontology.org/2025-03-16/ontology/go-basic.owl```

3. Run framework:

```python main.py --assembly_name {assembly name} --annotation_file {annotation file} --groovy (optional)```

---
## Input

Input files should be TAB-separated, with the protein ID in the first column, followed by TAB-separated GO:XXXXXXX classes.  

*Examples of input file: [examples](examples/input)*

---

## Outputs

| Filename                                | Description                                                                                   |
|-----------------------------------------|-----------------------------------------------------------------------------------------------|
| `{assembly_name}_report.html`           | HTML-formatted report with detailed metrics, figures, and searchable tables                  |
| `{assembly_name}_report.json`           | JSON-formatted report with detailed metrics                                                  |
| `{assembly_name}_consistency.tsv`       | 'Never in taxon' and 'only in taxon' constraints for each protein and GO annotation          |
| `{assembly_name}_taxon_explanations.tsv`| Taxonomic consistency satisfiability and HSTExplanationGenerator explanations (if any)       |
| `{assembly_name}_IC.tsv`                | Information Content (IC) for each GO class                                                   |
| `{assembly_name}_specific_GO_terms.tsv` | Most specific GO classes retained for each protein                                           |

*Examples of output files: [examples](examples/outputs)*

---

### Generating constraint files

Constraint files in the repository were generated using GO release 2025-03-16. To generate the constraint files for a different GO version:


#### Pathway constraints

1. Download go.obo:  
```wget https://release.geneontology.org/{GO_VERSION}/ontology/go.obo```  

2. Map GO classes to MetaCyc pathway IDs:  
```python generate_constraints/GO2Metacyc.py {GO_VERSION} {/PATH/TO/go.obo}```  

3. Extract pathway EC numbers:  
```python generate_constraints/extract_ECs.py {GO_VERSION}```  

4. Download EC2GO:  
```wget http://www.geneontology.org/external2go/ec2go```  

#### Taxon constraints

1. Download GO constraint files and NCBI Taxonomy:  
```wget https://release.geneontology.org/{VERSION}/ontology/imports/go-taxon-groupings.owl```  
```wget https://release.geneontology.org/{VERSION}/ontology/imports/go-computed-taxon-constraints.obo```  
```wget http://purl.obolibrary.org/obo/ncbitaxon/{VERSION}/ncbitaxon.owl```  

2. Create constraint files:  
```groovy generate_constraints/extract_constraints.groovy```  
```groovy generate_constraints/add_taxon_disjointness.groovy```
