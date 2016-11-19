#include "graphutils.h"
#include "GraphHelper.h"
#include "verifier.h"
#include <iomanip> // Pretty debug output
#include <set>
#include <vector>
#include <queue>

Graph GraphHelper::generateRandomGraph(int numNodes, float density) {
	Graph g(numNodes);
	srand(static_cast<unsigned int>(time(0)));
	int randNodes = rand() % numNodes + 1;
	int numEdges = static_cast<int>(randNodes * (numNodes - randNodes) * density);

	if (numEdges == 0)
		std::cerr << "Warning! Graph has no edges." << std::endl;

	for (int i = 0; i < numEdges; ++i) {
		bool inserted = false;
		while (!inserted) {
			int randLeft = rand() % randNodes;
			int randRight = rand() % (numNodes - randNodes) + randNodes;

			if (!boost::edge(randLeft, randRight, g).second) {
				boost::add_edge(randLeft, randRight, g);
				inserted = true;
			}
		}
	}

	return g;
}

void GraphHelper::readGraphFromFile(Graph& g, Vertex& first_right, const std::string& filePath) {
	std::ifstream inFile(filePath);
	if (inFile.fail()) std::cerr << "Error while reading from file" << std::endl;

	int numNodes, u, v;

	inFile >> numNodes;
	inFile >> first_right;
	g = Graph(numNodes);

	while (inFile >> u >> v) boost::add_edge(u, v, g);
}

void GraphHelper::writeGraphToFile(const std::string& filePath, const Graph& g) {
	std::ofstream outFile(filePath);
	outFile << num_vertices(g) << std::endl;

	for (EdgeIterator e = boost::edges(g).first; e != boost::edges(g).second; e++)
		outFile << source(*e, g) << " " << target(*e, g) << std::endl;

	outFile.flush();
	if (outFile.fail())
		std::cerr << "Error while writing to file." << std::endl;
}

bool GraphHelper::isMaximumMatching(const VertexVector& matching, const Graph& g) {
	try {
		vertex_size_t n = num_vertices(g);
		VertexVector solution_mates(n);
		boost::edmonds_maximum_cardinality_matching(g, &solution_mates[0]);
		vertex_size_t matching_size_solution = boost::matching_size(g, &solution_mates[0]);
		verify_matching(g, matching, matching_size_solution);
	}
	catch (std::string error) {
		std::cout << error << std::endl;
		return false;
	}
	return true;
}

VertexVector GraphHelper::greedyMatching(const Graph& g) {

	VertexVector matching(num_vertices(g));

	Vertex null_vertex = g.null_vertex();

	// set all mates to null vector
	std::fill(matching.begin(), matching.end(), null_vertex);

	// do greedy matching over all edges
	EdgeIterator start, end;
	for (std::tie(start, end) = boost::edges(g); start != end; ++start) {
		Edge e = *start;
		Vertex u = boost::source(e, g);
		Vertex v = boost::target(e, g);
		if (matching[u] == null_vertex && matching[v] == null_vertex) {
			matching[u] = v;
			matching[v] = u;
		}
	}

	return matching;
}

void markAdjacentEdges(Vertex v, const Graph& g, std::set<Edge>& _edges, std::vector<size_t>& degree) {
	auto startOutEdges = boost::out_edges(v, g).first;
	auto endOutEdges = boost::out_edges(v, g).second;
	for (auto o = startOutEdges; o != endOutEdges; ++o) {
		if (_edges.find(*o) == _edges.end()) continue;

		Vertex a = source(*o, g);
		Vertex b = target(*o, g);

		_edges.erase(_edges.find(*o));
		--degree[a];
		--degree[b];
	}
}

VertexVector GraphHelper::karpSipser(const Graph& g) {
	// Initialize matching with null
	size_t numVertices = num_vertices(g);
	size_t numEdges = num_edges(g);
	VertexVector matching(numVertices);
	std::fill(matching.begin(), matching.end(), g.null_vertex());

	std::set<Edge> _edges;
	auto es = boost::edges(g);
	for (auto eit = es.first; eit != es.second; ++eit) {
		_edges.insert(*eit);
	}
	
	std::vector<size_t> degree(numVertices);
	for (Vertex v = 0; v < numVertices; ++v)
		degree[v] = boost::degree(v, g);

	while (!_edges.empty()) {
		// Search for edge with a vertex of degree 1
		bool foundEdge = false;

		for (auto e = boost::edges(g).first; e != boost::edges(g).second; ++e) {
			if (_edges.find(*e) == _edges.end()) {
				continue;
			}

			Vertex u = source(*e, g);
			Vertex v = target(*e, g);

			if (degree[v] == 1 || degree[u] == 1) {
				// Add vertices to matching
				matching[u] = v;
				matching[v] = u;

				// Remove out_edges from u,v
				markAdjacentEdges(u, g, _edges, degree);
				markAdjacentEdges(v, g, _edges, degree);

				foundEdge = true;
				break;
			}
		}

		if (foundEdge) continue;

		// If no edge found, select a random one
		int randEdge = rand() % _edges.size();
		auto e = _edges.begin();
		for (int i = 0; i < randEdge; ++i)
			++e;

		Vertex u = source(*e, g);
		Vertex v = target(*e, g);

		// Add vertices to matching
		matching[u] = v;
		matching[v] = u;

		// Remove out_edges from u,v
		markAdjacentEdges(u, g, _edges, degree);
		markAdjacentEdges(v, g, _edges, degree);
	}
	return matching;
}

void GraphHelper::printOutput(const BenchmarkResult& result) {
    std::cout << result.timeStamp << ",";
    std::cout << result.graphName << ",";
    std::cout << result.numVertices << ",";
    std::cout << result.numEdges << ",";
    std::cout << result.algorithm << ",";
    std::cout << result.numThreads << ",";    
    
    for (double d : result.durations) 
        std::cout << "," << d;
    
    std::cout << std::endl;
}



VertexVector GraphHelper::ks(const Graph& g) {

	auto n = boost::num_vertices(g);
	auto null_vertex = g.null_vertex();

	VertexVector matching(n);
	std::vector<Vertex> deg(n);
	for (Vertex v = 0; v < n; v++) {
		matching[v] = null_vertex;
		deg[v] = boost::degree(v, g);
	}

	bool found_edges;
    do {

		// match all vertices with degree one
		bool found_deg_one;
		do {
			found_deg_one = false;

			for (Vertex v = 0; v < n; v++) {
				if (deg[v] != 1 || is_matched(v, g, matching)) continue;

				AdjVertexIterator s, e;
				for (std::tie(s, e) = boost::adjacent_vertices(v, g); s != e; s++) {
					Vertex u = *s;
					if (is_matched(u, g, matching)) continue;

					matching[v] = u;
					matching[u] = v;
					deg[u]--;
					deg[v]--;
					found_deg_one = true;
				}
			}

		} while (found_deg_one);

		found_edges = false;
		// take first unmatched edge
		EdgeIterator estart, eend;
		for (std::tie(estart, eend) = boost::edges(g); estart != eend; estart++) {
			Edge e = *estart;
			Vertex u = boost::source(e, g);
			Vertex v = boost::target(e, g);

			if (is_matched(u,g,matching) || is_matched(v,g,matching)) continue;

			matching[v] = u;
			matching[u] = v;
			deg[u]--;
			deg[v]--;
			found_edges = true;
		}

	} while (found_edges);

	return matching;

}
