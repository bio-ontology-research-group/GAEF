// Updated Grapes section
@Grapes([
  @Grab(group='net.sourceforge.owlapi', module='owlapi-distribution', version='5.5.0'),
  @Grab(group='io.github.liveontologies', module='elk-owlapi', version='0.6.0'),
  @GrabConfig(systemClassLoader=true)
])

import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model.*
import org.semanticweb.elk.owlapi.ElkReasonerFactory
import org.semanticweb.owlapi.reasoner.OWLReasoner
import com.clarkparsia.owlapi.explanation.BlackBoxExplanation
import com.clarkparsia.owlapi.explanation.HSTExplanationGenerator
import groovy.transform.Immutable

@Immutable
class GenomeAnnotation {
    String genomeName
    String proteinName
    String goId
    Set<String> neverInTaxon
    Set<String> onlyInTaxon
}

@Immutable
class Result {
    boolean isSatisfiable
    String explanation = ""
}

if (args.length != 2) {
    println "Usage: groovy check_constraints.groovy <input_file> <OUTPUT_FILE>"
    System.exit(1)
}

def manager = OWLManager.createOWLOntologyManager()
def df = manager.getOWLDataFactory()
def baseIRI = "http://purl.obolibrary.org/obo/"
def config = manager.getOntologyLoaderConfiguration()
    .setMissingImportHandlingStrategy(MissingImportHandlingStrategy.SILENT)
    .setLoadAnnotationAxioms(false)

println "Loading ontologies..."
def taxonFile = new File("ncbitaxon_with_disjointness.owl")
def goTaxonFile = new File("go-taxon-groupings.owl")

if (!taxonFile.exists() || !goTaxonFile.exists()) {
    println "Error: Ontology file(s) not found"
    System.exit(1)
}

def taxonOntology = manager.loadOntologyFromOntologyDocument(new org.semanticweb.owlapi.io.FileDocumentSource(taxonFile), config)
def goTaxonOntology = manager.loadOntologyFromOntologyDocument(new org.semanticweb.owlapi.io.FileDocumentSource(goTaxonFile), config)
def mergedIRI = IRI.create("http://merged.ontology/taxon-go")
def ontology = manager.createOntology(mergedIRI)
manager.addAxioms(ontology, taxonOntology.getAxioms())
manager.addAxioms(ontology, goTaxonOntology.getAxioms())

println "Ontologies loaded and merged successfully."

def inputFile = new File(args[0])
def outputFile = new File(args[1])

if (!inputFile.exists()) {
    println "Error: Input file ${args[0]} does not exist"
    System.exit(1)
}

outputFile.withWriter { writer ->
    writer.println("Genome\tIsSatisfiable\tExplanation")
}

Map<String, List<GenomeAnnotation>> genomeAnnotationsMap = [:]
boolean isFirstLine = true
inputFile.eachLine { line ->
    if (isFirstLine) {
        isFirstLine = false
        return
    }
    def parts = line.split("\t", -1)
    if (parts.length < 5) return
    def genomeName = parts[0].trim().replaceAll(".tsv", "")
    def annotation = new GenomeAnnotation(
        genomeName: genomeName,
        proteinName: parts[1].trim(),
        goId: parts[2].trim(),
        neverInTaxon: parts[3].trim() ? parts[3].split(",").collect { it.trim() } as Set : [] as Set,
        onlyInTaxon: parts[4].trim() ? parts[4].split(",").collect { it.trim() } as Set : [] as Set
    )
    genomeAnnotationsMap.computeIfAbsent(genomeName) { [] } << annotation
}

println "Found ${genomeAnnotationsMap.size()} genomes."

def taxonIRI = { tid -> IRI.create("${baseIRI}${tid}") }
def negatedTaxonIRI = { tid -> IRI.create("${baseIRI}${tid}_neg") }
def getTaxonIdFromIRI = { iri -> iri.getFragment()?.replace("NCBITaxon_", "").replace("NCBITaxon_Union_", "").replace("_neg", "") ?: iri.toString() }
def describeRelationship = { axiom, id1, id2 ->
    if (axiom instanceof OWLDisjointClassesAxiom) return "T${id1} and T${id2} are disjoint"
    if (axiom instanceof OWLSubClassOfAxiom) {
        def subClass = getTaxonIdFromIRI(axiom.getSubClass().asOWLClass().getIRI())
        def superClass = getTaxonIdFromIRI(axiom.getSuperClass().asOWLClass().getIRI())
        if (subClass == id1 && superClass == id2) return "T${id1} is a subclass of T${id2}"
        if (subClass == id2 && superClass == id1) return "T${id2} is a subclass of T${id1}"
    }
    return "Relationship between T${id1} and T${id2}: ${axiom.getAxiomType()}"
}

def normalizeTaxonId = { tid -> tid.replace("NCBITaxon_", "").replace("NCBITaxon_Union_", "") }

def reasonerFactory = new ElkReasonerFactory()
def addedAxioms = []
def reasoner = reasonerFactory.createReasoner(ontology)

try {
    genomeAnnotationsMap.each { genomeName, annotations ->
        println "\nProcessing genome: ${genomeName}"
        def onlySet = annotations.collectMany { it.onlyInTaxon.collect(normalizeTaxonId) } as Set
        def neverSet = annotations.collectMany { it.neverInTaxon.collect(normalizeTaxonId) } as Set

        def onlyClasses = onlySet.collect { df.getOWLClass(taxonIRI("NCBITaxon_" + it)) }
        def neverClasses = neverSet.collect { df.getOWLClass(negatedTaxonIRI("NCBITaxon_" + it)) }

        def expr = df.getOWLObjectIntersectionOf((onlyClasses + neverClasses) ?: [df.getOWLThing()])
        def conceptIRI = IRI.create("http://example.org#${genomeName}_Taxon_Concept")
        def conceptClass = df.getOWLClass(conceptIRI)
        def axiom = df.getOWLEquivalentClassesAxiom(conceptClass, expr)

        manager.addAxiom(ontology, axiom)
        addedAxioms << axiom

        reasoner.flush()
        boolean isSatisfiable = reasoner.isSatisfiable(conceptClass)
        String explanationText = ""

        if (!isSatisfiable) {
            try {
                def blackBox = new BlackBoxExplanation(ontology, reasonerFactory, reasoner)
                def hstGenerator = new HSTExplanationGenerator(blackBox)
                def explanation = hstGenerator.getExplanation(conceptClass)

                if (explanation && !explanation.isEmpty()) {
                    def conflictingTaxa = [:]
                    def relationshipAxiom = null
                    explanation.each { ax ->
                        ax.getClassesInSignature().each { cls ->
                            def iri = cls.getIRI()
                            def taxonId = getTaxonIdFromIRI(iri)
                            boolean isNeg = iri.toString().contains("_neg")
                            def type = isNeg ? "never" : "only"
                            if ((type == "only" && onlySet.contains(taxonId)) || (type == "never" && neverSet.contains(taxonId))) {
                                conflictingTaxa[taxonId] = conflictingTaxa.getOrDefault(taxonId, [only: [] as Set, never: [] as Set])
                                def annots = annotations.findAll {
                                    (type == "only" && it.onlyInTaxon.collect(normalizeTaxonId).contains(taxonId)) ||
                                    (type == "never" && it.neverInTaxon.collect(normalizeTaxonId).contains(taxonId))
                                }
                                conflictingTaxa[taxonId][type].addAll(annots)
                            }
                        }
                        if (!relationshipAxiom && (ax instanceof OWLDisjointClassesAxiom || ax instanceof OWLSubClassOfAxiom)) {
                            relationshipAxiom = ax
                        }
                    }
                    def parts = []
                    conflictingTaxa.each { taxonId, data ->
                        data.each { type, annots ->
                            if (annots) parts << annots.collect { "Protein ${it.proteinName} (${it.goId})" }.join(", ") + " requires ${type} in taxon T${taxonId}"
                        }
                    }
                    explanationText = parts.join("; ")
                    if (relationshipAxiom) {
                        def ids = relationshipAxiom.getClassesInSignature().collect { getTaxonIdFromIRI(it.getIRI()) }
                        explanationText += ", and ${describeRelationship(relationshipAxiom, ids[0], ids[1])}"
                    }
                } else {
                    explanationText = "Explanation could not be generated."
                }
            } catch (Exception ex) {
                explanationText = "Error during explanation: ${ex.message}"
            }
        }

        def result = new Result(isSatisfiable: isSatisfiable, explanation: explanationText)
        outputFile.append("${genomeName}\t${result.isSatisfiable}\t${result.explanation ?: ''}\n")

        manager.removeAxiom(ontology, axiom)
        addedAxioms.remove(axiom)
    }
} finally {
    if (!addedAxioms.isEmpty()) manager.removeAxioms(ontology, addedAxioms as Set)
    if (reasoner) reasoner.dispose()
    println "\nFinished processing genomes."
}
