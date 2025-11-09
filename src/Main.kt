import kotlin.math.max

fun <A, B> cartesianProduct(a: List<A>, b: List<B>): List<Pair<A, B>> = a.flatMap { x -> b.map { y -> x to y } }

fun <A, B> cartesianProductSeq(a: List<A>, b: List<B>): Sequence<Pair<A, B>> =
  a.asSequence().flatMap { x -> b.asSequence().map { y -> x to y } }

fun main() {
  val (bigGraph, smallGraph) = readTwoAdjacencyMatrices("/Users/ernestmolczan/IdeaProjects/ApproxMinimalExtensionGraph/src/sample_graphs5.txt")

  val n = bigGraph.size
  val bigGraphVertices = (0 until n).toList()
  val k = smallGraph.size
  val smallGraphVertices = (0 until k).toList()
  val m = 3

  val allMappingsOfPVerticesToGVertices = cartesianProduct(smallGraphVertices, bigGraphVertices)
  val allMappingsOfPVerticesToGVerticesSeq = cartesianProductSeq(smallGraphVertices, bigGraphVertices)

  val missingEdgesMatrix: Array<Array<MutableList<Pair<Int, Int>>>> = Array(k) {
    Array(n) {
      mutableListOf()
    }
  }

  val fullMappings: Array<Array<MutableList<Pair<MutableList<Int>, MutableList<Int>>>>> = Array(k) {
    Array(n) {
      mutableListOf()
    }
  }

  // For each mapping of a single vertex of smallGraph to a single vertex of bigGraph
  allMappingsOfPVerticesToGVertices.forEach { globalMapping ->
    val mappedVerticesSmall = mutableListOf<Int>()
    mappedVerticesSmall.add(globalMapping.first)

    val mappedVerticesBig = mutableListOf<Int>()
    mappedVerticesBig.add(globalMapping.second)

    // Greedily extend the mapping
    while (mappedVerticesSmall.size < k) {
      val currentlyPossibleBestNeighborsSmall = smallGraphVertices.filter { v -> !mappedVerticesSmall.contains(v) }
      val currentlyPossibleBestNeighborsBig = bigGraphVertices.filter { v -> !mappedVerticesBig.contains(v) }
      val allMappingsOfNeighbors =
        cartesianProduct(currentlyPossibleBestNeighborsSmall, currentlyPossibleBestNeighborsBig)

      var bestMappingOfNeighbor = allMappingsOfNeighbors[0]
      var bestOutgoingCost = max(
        0,
        smallGraph[mappedVerticesSmall.last()][bestMappingOfNeighbor.first] - bigGraph[mappedVerticesBig.last()][bestMappingOfNeighbor.second]
      )
      var bestIncomingCost = max(
        0,
        smallGraph[bestMappingOfNeighbor.first][mappedVerticesSmall.last()] - bigGraph[bestMappingOfNeighbor.second][mappedVerticesBig.last()]
      )
      var bestCostForNeighborMapping = bestOutgoingCost + bestIncomingCost

      // Find the best neighbor mapping
      allMappingsOfNeighbors.forEach { mapping ->
        val currentOutgoingCost = max(
          0, smallGraph[mappedVerticesSmall.last()][mapping.first] - bigGraph[mappedVerticesBig.last()][mapping.second]
        )
        val currentIncomingCost = max(
          0, smallGraph[mapping.first][mappedVerticesSmall.last()] - bigGraph[mapping.second][mappedVerticesBig.last()]
        )
        val currentCostForNeighborMapping = currentOutgoingCost + currentIncomingCost

        if (currentCostForNeighborMapping < bestCostForNeighborMapping) {
          bestOutgoingCost = currentOutgoingCost
          bestIncomingCost = currentIncomingCost
          bestCostForNeighborMapping = currentCostForNeighborMapping
          bestMappingOfNeighbor = mapping
        }
      }

      // Update the mapping and missing edges
      mappedVerticesSmall.add(bestMappingOfNeighbor.first)
      mappedVerticesBig.add(bestMappingOfNeighbor.second)
      repeat(bestOutgoingCost) {
        missingEdgesMatrix[globalMapping.first][globalMapping.second].add(
          Pair(
            mappedVerticesBig[mappedVerticesBig.size - 2],
            mappedVerticesBig[mappedVerticesBig.size - 1]
          )
        )
      }
      repeat(bestIncomingCost) {
        missingEdgesMatrix[globalMapping.first][globalMapping.second].add(
          Pair(
            mappedVerticesBig[mappedVerticesBig.size - 1],
            mappedVerticesBig[mappedVerticesBig.size - 2]
          )
        )
      }
    }
    // Store the full mapping
    fullMappings[globalMapping.first][globalMapping.second].add(Pair(mappedVerticesSmall, mappedVerticesBig))
  }

  // For now, we have only one best mapping per cell

  // For each mapping, let's find missing edges
  val missingEdgesMatrixForRealNow: Array<Array<MutableList<Pair<Int, Int>>>> = Array(k) {
    Array(n) {
      mutableListOf()
    }
  }

  fullMappings.forEach { row ->
    row.forEach { mapping ->
      mapping.forEach { mappingPair ->
        val mappedVerticesSmall = mappingPair.first
        val mappedVerticesBig = mappingPair.second

        // Create the trimmed big graph adjacency matrix for the current combination
        val trimmedBigMatrix = Array(k) { Array(k) { 0 } }
        for (row in 0 until k) {
          for (col in 0 until k) {
            trimmedBigMatrix[row][col] = bigGraph[mappedVerticesBig[row]][mappedVerticesBig[col]]
          }
        }

        // Create the small graph adjacency matrix for the current mapping
        val permutatedSmallMatrix = Array(k) { Array(k) { 0 } }
        for (row in 0 until k) {
          for (col in 0 until k) {
            permutatedSmallMatrix[row][col] = smallGraph[mappedVerticesSmall[row]][mappedVerticesSmall[col]]
          }
        }

        for (i in 0 until k) {
          for (j in 0 until k) {
            // Calculate missing edges, no negative values
            val numberOfMissingEdges = max(0, permutatedSmallMatrix[i][j] - trimmedBigMatrix[i][j])
            if (numberOfMissingEdges > 0) {
              repeat(numberOfMissingEdges) {
                missingEdgesMatrixForRealNow[fullMappings.indexOf(row)][row.indexOf(mapping)].add(
                  Pair(
                    mappedVerticesBig[i],
                    mappedVerticesBig[j]
                  )
                )
              }
            }
          }
        }
      }
    }
  }

  var minimalNumberOfMissingEdges = missingEdgesMatrixForRealNow[0][0].size
  var minRow = 0
  var minCol = 0
  missingEdgesMatrixForRealNow.forEachIndexed { indexOfRow, row ->
    row.forEachIndexed { indexOfCell, cell ->
      if (cell.size < minimalNumberOfMissingEdges) {
        minimalNumberOfMissingEdges = cell.size
        minRow = indexOfRow
        minCol = indexOfCell
      }
    }
  }

  println("Minimal number of missing edges to add: $minimalNumberOfMissingEdges")
  println("Edges to add: ${missingEdgesMatrixForRealNow[minRow][minCol]}")
  println("Mapping of small graph vertices to big graph vertices: ${fullMappings[minRow][minCol][0]}")

  // now let's choose best m mappings

}