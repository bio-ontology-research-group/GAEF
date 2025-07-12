@Grab(group='net.sourceforge.owlapi', module='owlapi-distribution', version='5.1.20')
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model.*
import java.nio.file.*

def parseGoConstraints(String oboFile) {
    def constraints = []
    def currentTerm = null

    new File(oboFile).eachLine { line ->
        line = line.trim()

        if (line.startsWith('id: ')) {
            currentTerm = line.substring(4)
        }
        else if (line.startsWith('property_value: RO:0002161') && currentTerm) {
            def taxon = line =~ /(NCBITaxon(?:_Union)?):([\d]+)/
            if (taxon) {
                constraints << [currentTerm, "never_in_taxon", "${taxon[0][1]}_${taxon[0][2]}"]
            }
        }
        else if (line.startsWith('relationship: RO:0002162') && currentTerm) {
            def taxon = line =~ /(NCBITaxon(?:_Union)?):([\d]+)/
            if (taxon) {
                constraints << [currentTerm, "only_in_taxon", "${taxon[0][1]}_${taxon[0][2]}"]
            }
        }
    }
    return constraints
}

def parseGoHierarchy(String oboFile) {
    def parents = [:].withDefault { [] }
    def children = [:].withDefault { [] }
    def currentTerm = null

    new File(oboFile).eachLine { line ->
        line = line.trim()
        if (line.startsWith('id: ')) {
            currentTerm = line.substring(4)
        } else if (line.startsWith('is_a: ') && currentTerm) {
            def parentId = (line =~ /GO:\d+/)[0]
            parents[currentTerm] << parentId
            children[parentId] << currentTerm
        }
        // Optional: handle part_of
        else if (line.startsWith('relationship: part_of') && currentTerm) {
            def parentId = (line =~ /GO:\d+/)[0]
            parents[currentTerm] << parentId
            children[parentId] << currentTerm
        }
    }
    return [parents: parents, children: children]
}

def propagateConstraints(constraints, hierarchy) {
    def propagated = [] as Set

    propagated.addAll(constraints)

    constraints.findAll { it[1] == "never_in_taxon" }.each { constraint ->
        def queue = [constraint[0]]
        def visited = [] as Set

        while (queue) {
            def term = queue.pop()
            if (!visited.contains(term)) {
                visited.add(term)
                propagated.add([term, constraint[1], constraint[2]])
                queue.addAll(hierarchy.parents[term] ?: [])
            }
        }
    }

    constraints.findAll { it[1] == "only_in_taxon" }.each { constraint ->
        def queue = [constraint[0]]
        def visited = [] as Set

        while (queue) {
            def term = queue.pop()
            if (!visited.contains(term)) {
                visited.add(term)
                propagated.add([term, constraint[1], constraint[2]])
                queue.addAll(hierarchy.children[term] ?: [])
            }
        }
    }

    return propagated.toList()
}

// Main pipeline

def directConstraints = parseGoConstraints("go-computed-taxon-constraints.obo")
def hierarchy = parseGoHierarchy("go-computed-taxon-constraints.obo")
def propagatedConstraints = propagateConstraints(directConstraints, hierarchy)

def output = new File("go_taxon_constraints_updated.tsv")
output.withPrintWriter { writer ->
    writer.println("GO_ID\tConstraint_Type\tTaxon_ID")
    propagatedConstraints.each { row ->
        writer.println("${row[0]}\t${row[1]}\t${row[2]}")
    }
}

println "Direct constraints: ${directConstraints.size()}"
println "After propagation: ${propagatedConstraints.size()}"

