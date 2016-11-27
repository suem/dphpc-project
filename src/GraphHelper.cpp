#include "graphutils.h"
#include "GraphHelper.h"
#include "verifier.h"
#include <iomanip> // Pretty debug output
#include <set>
#include <vector>
#include <queue>

Graph GraphHelper::generateRandomGraph(int numNodes, float density, Vertex& first_right) {
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
	
	first_right = randNodes;
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

void GraphHelper::printOutput(const BenchmarkResult& resultStruct, const std::string& outPath) {
	// generate file name: algorithm name, graphName, and timestamp
	std::string fileName = resultStruct.algorithm + '_' + resultStruct.graphName + '_' + resultStruct.timeStamp + ".csv";
	std::string filePath = outPath + fileName;
	std::replace(filePath.begin(), filePath.end(), ' ', '_');
	// create csv file
	std::ofstream csv(filePath);


    // Header
    csv << "TimeStamp" << ",";
    csv << "GraphName" << ",";
    csv << "NumVertices" << ",";
    csv << "NumEdges" << ",";
    csv << "Algorithm" << std::endl;
	 
	csv << resultStruct.timeStamp << ",";
	csv << resultStruct.graphName << ",";
	csv << resultStruct.numVertices << ",";
	csv << resultStruct.numEdges << ",";
	csv << resultStruct.algorithm << std::endl;
    
    // Data
	std::string temp;
	for (int i : resultStruct.numThreads)
		temp += std::to_string(i) + ",";

	csv << temp.substr(0, temp.size() - 1) << std::endl;

    for (int i = 0; i < resultStruct.iter; ++i) {
        int pos = 0;
		temp.clear();
        for (int nThreads : resultStruct.numThreads) {
			temp += std::to_string(resultStruct.durations[pos][i]) + ",";
            ++pos;
        }
		csv << temp.substr(0, temp.size() - 1) << std::endl;
    } 
	std::cout << "CSV Filepath: " << filePath<< std::endl;
}

std::string GraphHelper::getGraphNameFromPath(const std::string& graphPath) {
	// This helperfunction takes a path to a graph file and returns the name of the file
	std::istringstream ss(graphPath);
	std::string token;
	std::string last;
	char delimiter;

	#ifdef _WIN32
		delimiter = '\\';
	#else
		delimiter = '/';
	#endif

	while (std::getline(ss, token, delimiter)) {
		last = token;
	}

	std::istringstream ss2(last);

	std::getline(ss2, token, '.');
	return token;
}

void markAdjacentEdges(Vertex v, const Graph& g, std::vector<size_t>& degree) {
    auto startOutEdges = boost::out_edges(v, g).first;
    auto endOutEdges = boost::out_edges(v, g).second;
    degree[v] = 0;
    for (auto o = startOutEdges; o != endOutEdges; ++o) {
        Vertex b = target(*o, g);
        --degree[b];
    }
}

VertexVector GraphHelper::karpSipser(const Graph& g) {
    auto const n = boost::num_vertices(g);
    auto const null_vertex = Graph::null_vertex();

    VertexVector matching(n);
    std::vector<Vertex> deg(n);
    for (Vertex v = 0; v < n; v++) {
        matching[v] = null_vertex;
        deg[v] = boost::degree(v, g);
    }

    bool found_edges;

    do {
        found_edges = false;

        for (auto e = boost::edges(g).first; e != boost::edges(g).second; ++e) {
            Vertex u = boost::source(*e, g);
            Vertex v = boost::target(*e, g);

            if (is_matched(u, matching) || is_matched(v, matching)) {
                continue;
            }

            if (deg[v] == 1 || deg[u] == 1) {
                matching[v] = u;
                matching[u] = v;
                markAdjacentEdges(v, g, deg);
                markAdjacentEdges(u, g, deg);
                found_edges = true;
                break;
            }
        }

        if (found_edges) continue;

        // If no edge found, select a random unmatched one
        EdgeIterator eStart, eEnd, eRandom, eIt;
        std::tie(eStart, eEnd) = boost::edges(g);
        int randEdge = rand() % num_edges(g);
        eRandom = eStart;
        for (int i = 0; i < randEdge; i++)
            ++eRandom;

        eIt = eRandom;
        do {
            Edge e = *eIt;
            Vertex u = boost::source(e, g);
            Vertex v = boost::target(e, g);

            if (!is_matched(u, matching) && !is_matched(v, matching)) {
                matching[v] = u;
                matching[u] = v;
                markAdjacentEdges(v, g, deg);
                markAdjacentEdges(u, g, deg);
                found_edges = true;
                break;
            }

            ++eIt;
            if (eIt == eEnd) {
                eIt = eStart;
            }

        } while (eIt != eRandom);

    } while (found_edges);

    return matching;
}

VertexVector GraphHelper::ks(const Graph& g) {

    auto const n = boost::num_vertices(g);
    auto const null_vertex = Graph::null_vertex();

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
                if (deg[v] != 1 || is_matched(v, matching)) continue;

                AdjVertexIterator s, e;
                for (std::tie(s, e) = boost::adjacent_vertices(v, g); s != e; s++) {
                    Vertex u = *s;
                    if (is_matched(u, matching)) continue;

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

            if (is_matched(u, matching) || is_matched(v, matching)) continue;

            matching[v] = u;
            matching[u] = v;
            deg[u]--;
            deg[v]--;
            found_edges = true;
        }

    } while (found_edges);

    return matching;
}


void matchAndUpdate(
        Vertex u, 
        const Graph& g, 
        VertexVector& matching, 
        std::atomic_int* deg,
        std::atomic_flag* visited) {

    if (visited[u].test_and_set()) return;

    AdjVertexIterator start, end;
    std::tie(start, end) = boost::adjacent_vertices(u, g); 
    Vertex v;
    while (start != end) {
        v = *start;
        if (!visited[v].test_and_set()) {
            break;
        }
        start++;
    }

    if (start == end) {
        return;
    }

    // add (u,v) to matching
    matching[u] = v;
    matching[v] = u;

    for (std::tie(start, end) = boost::adjacent_vertices(v, g); start != end; start++) {
        Vertex w = *start;
        if (std::atomic_fetch_add(&deg[w], -1) == 2) {
            matchAndUpdate(w, g, matching, deg, visited);
        }
    }

}

VertexVector GraphHelper::parallelKarpSipser(const Graph& g, const Vertex first_right) {
    vertex_size_t n = boost::num_vertices(g);
    VertexVector matching(n, g.null_vertex());

    std::atomic_flag* visited = new std::atomic_flag[n];
    memset(visited, 0, sizeof(std::atomic_flag) * n);

    std::atomic_int* deg = new std::atomic_int[first_right];

#pragma omp parallel for
    for (int x = 0; x < first_right; x++) {
        deg[x] = boost::degree(x, g);
    }

#pragma omp parallel for
    for (int x = 0; x < first_right; x++) 
        if (deg[x] == 1) matchAndUpdate(x, g, matching, deg, visited);

#pragma omp parallel for
    for (int x = 0; x < first_right; x++) 
        if (deg[x] > 1) matchAndUpdate(x, g, matching, deg, visited);

    return matching;
}


