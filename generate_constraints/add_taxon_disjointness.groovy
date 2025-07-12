@GrabResolver(name='mavenCentral', root='https://repo1.maven.org/maven2/')
@GrabConfig(systemClassLoader=true)
@Grab(group='net.sourceforge.owlapi', module='owlapi-distribution', version='5.5.0')
@Grab(group='org.slf4j', module='slf4j-simple', version='2.0.9')
@Grab(group='io.github.liveontologies', module='elk-owlapi', version='0.6.0')

import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.util.DefaultPrefixManager
import org.semanticweb.owlapi.formats.RDFXMLDocumentFormat
import org.semanticweb.owlapi.search.EntitySearcher
import org.semanticweb.elk.owlapi.ElkReasonerFactory
import org.semanticweb.owlapi.reasoner.InferenceType

// Parse command line arguments
def cli = new CliBuilder(usage: 'groovy add_taxon_disjointness.groovy [options]')
cli.with {
    i(longOpt: 'input', args: 1, argName: 'file', 'Input TSV file with taxon constraints (default: taxon_constraints.tsv)')
    o(longOpt: 'output', args: 1, argName: 'file', 'Output OWL file (default: ncbitaxon_with_disjointness.owl)')
    t(longOpt: 'taxonomy', args: 1, argName: 'file', 'NCBI Taxonomy OWL file (default: http://purl.obolibrary.org/obo/ncbitaxon.owl)')
    h(longOpt: 'help', 'Show usage information')
}

def options = cli.parse(args)
if (options.h) {
    cli.usage()
    return
}

// Define file paths
def inputFilePath = options.i ?: "go_taxon_constraints.tsv"
def outputFilePath = options.o ?: "ncbitaxon_with_disjointness.owl"
def ncbitaxonPath = options.t ?: "http://purl.obolibrary.org/obo/ncbitaxon.owl"

// Create IRI for ncbitaxon
def ncbitaxonIRI
if (ncbitaxonPath.startsWith("http://") || ncbitaxonPath.startsWith("https://") || ncbitaxonPath.startsWith("file:")) {
    ncbitaxonIRI = IRI.create(ncbitaxonPath)
} else {
    // Assume it's a local file
    def ncbitaxonFile = new File(ncbitaxonPath)
    if (!ncbitaxonFile.exists()) {
        println "Error: Taxonomy file not found: $ncbitaxonPath"
        return
    }
    ncbitaxonIRI = IRI.create(ncbitaxonFile.toURI())
}

println "Using taxonomy: $ncbitaxonIRI"
println "Input file: $inputFilePath"
println "Output file: $outputFilePath"

// Create OWL ontology manager and data factory
def manager = OWLManager.createOWLOntologyManager()
def dataFactory = manager.getOWLDataFactory()

// Create a new ontology that will import ncbitaxon.owl
def ontologyIRI = IRI.create("http://example.org/ncbitaxon-disjoint.owl")
def ontology = manager.createOntology(ontologyIRI)

// Add import declaration for ncbitaxon.owl
def importDeclaration = dataFactory.getOWLImportsDeclaration(ncbitaxonIRI)
manager.applyChange(new AddImport(ontology, importDeclaration))

// Load ncbitaxon.owl to get taxonomy information
println "Loading ncbitaxon.owl..."
def ncbitaxon = manager.loadOntology(ncbitaxonIRI)

// Set up prefix manager
def prefixManager = new DefaultPrefixManager()
prefixManager.setPrefix("obo:", "http://purl.obolibrary.org/obo/")
prefixManager.setPrefix("NCBITaxon:", "http://purl.obolibrary.org/obo/NCBITaxon_")

// Initialize ELK reasoner
println "Initializing ELK reasoner..."
def reasonerFactory = new ElkReasonerFactory()
def reasoner = reasonerFactory.createReasoner(ncbitaxon)
println "Precomputing inferences (this may take a while)..."
try {
    reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY)
} catch (Exception e) {
    println "Warning: Error during precomputation: ${e.message}"
    println "Continuing without precomputation..."
}

// Parse input file and collect taxon IDs
def taxonIds = new HashSet<String>()
def inputFile = new File(inputFilePath)

if (inputFile.exists()) {
    println "Reading input file: $inputFilePath"
    def lineCount = 0
    inputFile.eachLine { line ->
        lineCount++
        if (lineCount > 1) { // Skip header line
            def parts = line.split('\t')
            if (parts.size() >= 3) {
                taxonIds.add(parts[2])
            }
        }
    }
    println "Taxon IDs: $taxonIds"
} else {
    // Create example data from the provided sample
    println "Input file not found. Creating example data from sample."
    def sampleData = """GO_ID\tConstraint_Type\tTaxon_ID
GO:0000001\tonly_in_taxon\t2759
GO:0000011\tnever_in_taxon\t4896
GO:0000015\tonly_in_taxon\t131567
GO:0000022\tonly_in_taxon\t2759
GO:0000027\tonly_in_taxon\t131567
GO:0000028\tonly_in_taxon\t131567
GO:0000045\tonly_in_taxon\t2759
GO:0000054\tonly_in_taxon\t2759
GO:0000055\tonly_in_taxon\t2759
GO:0000056\tonly_in_taxon\t2759
GO:0000070\tonly_in_taxon\t2759
GO:0000073\tonly_in_taxon\t4751
GO:0000082\tonly_in_taxon\t2759
GO:0000086\tonly_in_taxon\t2759"""

    def lineCount = 0
    sampleData.eachLine { line ->
        lineCount++
        if (lineCount > 1) { // Skip header line
            def parts = line.split('\t')
            if (parts.size() >= 3) {
                taxonIds.add(parts[2])
            }
        }
    }
}

println "Found ${taxonIds.size()} unique taxon IDs"

// Function to get an OWL class for a taxon ID
def getTaxonClass = { String taxonId ->
    // Use the full IRI to avoid potential prefix issues
    def taxonIRI = IRI.create("http://purl.obolibrary.org/obo/NCBITaxon_" + taxonId)
    def taxonClass = dataFactory.getOWLClass(taxonIRI)
    // Check if the class exists in the ontology
    if (!ncbitaxon.containsEntityInSignature(taxonClass)) {
        println "Warning: Taxon ${taxonIRI.getShortForm()} not found in the ontology"
    }
    return taxonClass
}

// Function to get direct children of a taxon
def getDirectChildren = { OWLClass taxonClass ->
    try {
        return reasoner.getSubClasses(taxonClass, true).getFlattened().findAll {
            it != dataFactory.getOWLNothing()
        }
    } catch (Exception e) {
        println "Error getting children for ${taxonClass}: ${e.message}"
        return []
    }
}

// Function to get direct parent of a taxon
def getDirectParent = { OWLClass taxonClass ->
    try {
        return reasoner.getSuperClasses(taxonClass, true).getFlattened().findAll {
            it != dataFactory.getOWLThing()
        }
    } catch (Exception e) {
        println "Error getting parents for ${taxonClass}: ${e.message}"
        return []
    }
}

// Function to get siblings of a taxon (taxa with the same parent)
def getSiblings = { OWLClass taxonClass ->
    def parents = getDirectParent(taxonClass)
    def siblings = new HashSet<OWLClass>()

    parents.each { parent ->
        def children = getDirectChildren(parent)
        siblings.addAll(children)
    }

    siblings.remove(taxonClass)
    return siblings
}

// Add disjointness axioms for each taxon in the input file
println "Adding disjointness axioms..."
def addedSiblingAxioms = 0
def addedNegatedAxioms = 0

taxonIds.each { taxonId ->
    try {
        def taxonClass = getTaxonClass(taxonId)
        def taxonIRI = taxonClass.getIRI()

        // --- Create negated class and disjointness axiom ---
        def negatedIRI = IRI.create(taxonIRI.toString() + "_neg")
        def negatedClass = dataFactory.getOWLClass(negatedIRI)

        // Add declaration for the new class
        def declarationAxiom = dataFactory.getOWLDeclarationAxiom(negatedClass)
        manager.addAxiom(ontology, declarationAxiom)

        // Add disjointness axiom between original and negated class
        def disjointWithNegatedAxiom = dataFactory.getOWLDisjointClassesAxiom(taxonClass, negatedClass)
        manager.addAxiom(ontology, disjointWithNegatedAxiom)
        addedNegatedAxioms++
        println "Added disjointness: ${taxonIRI.getShortForm()} DisjointWith ${negatedIRI.getShortForm()}"
        // --- End negated class ---


        // --- Add disjointness axioms between siblings (Existing logic) ---
        def siblings = getSiblings(taxonClass)

        if (siblings.isEmpty()) {
            //println "No siblings found for ${taxonIRI.getShortForm()}" // Reduce verbosity
        } else {
            //println "Found ${siblings.size()} siblings for ${taxonIRI.getShortForm()}" // Reduce verbosity
        }

        // Add disjointness axioms between this taxon and its siblings
        siblings.each { sibling ->
            // Ensure we don't add disjointness with self if reasoner includes it somehow
            if (taxonClass != sibling) {
                def disjointAxiom = dataFactory.getOWLDisjointClassesAxiom(taxonClass, sibling)
                // Check if the axiom already exists (e.g., added from the sibling's perspective)
                // Note: ontology.containsAxiom might be slow for large ontologies.
                // A set could be used to track added pairs if performance is an issue.
                if (!ontology.containsAxiom(disjointAxiom)) {
                     manager.addAxiom(ontology, disjointAxiom)
                     addedSiblingAxioms++
                     //println "Added disjointness: ${taxonIRI.getShortForm()} DisjointWith ${sibling.getIRI().getShortForm()}" // Reduce verbosity
                }
            }
        }
        // --- End sibling disjointness ---

    } catch (Exception e) {
        println "Error processing taxon $taxonId: ${e.message}"
    }
}

println "Added $addedNegatedAxioms negated class disjointness axioms"
println "Added $addedSiblingAxioms sibling disjointness axioms"

// Save the ontology
def outputFile = new File(outputFilePath)
def documentFormat = new RDFXMLDocumentFormat()
// Ensure the prefix manager is available for saving
documentFormat.setPrefixManager(prefixManager)
// Add prefix for the base OBO namespace if not already present
if (prefixManager.getPrefix("obo:") == null) {
    prefixManager.setPrefix("obo:", "http://purl.obolibrary.org/obo/")
}


println "Saving ontology to $outputFilePath..."
manager.saveOntology(ontology, documentFormat, IRI.create(outputFile.toURI()))
println "Ontology saved to $outputFilePath"

// Clean up
reasoner.dispose()
