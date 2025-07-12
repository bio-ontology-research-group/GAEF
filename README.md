# Genome Annotation Evaluation Framework

This tool generates a static HTML and JSON report evaluating Gene Ontology (GO) class annotations using the following metrics: **completeness**, **coherence**, and **consistency**. It integrates multiple constraint sources (e.g., EC2GO mappings, GO structure, NCBITaxonomy).

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

## Outputs

| Filename                                | Description                                                                                   |
|-----------------------------------------|-----------------------------------------------------------------------------------------------|
| `{assembly_name}_report.html`           | HTML-formatted report with detailed metrics, figures, and searchable tables                  |
| `{assembly_name}_report.json`           | JSON-formatted report with detailed metrics                                                  |
| `{assembly_name}_consistency.tsv`       | 'Never in taxon' and 'only in taxon' constraints for each protein and GO annotation          |
| `{assembly_name}_taxon_explanations.tsv`| Taxonomic consistency satisfiability and HSTExplanationGenerator explanations (if any)       |
| `{assembly_name}_IC.tsv`                | Information Content (IC) for each GO class                                                   |
| `{assembly_name}_specific_GO_terms.tsv` | Most specific GO classes retained for each protein                                           |


