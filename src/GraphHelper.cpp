#include "GraphHelper.h"
#include "verifier.h"
#include <iomanip> // Pretty debug output

Graph GraphHelper::generateRandomGraph(int numNodes, float density) { 
    Graph g(numNodes);
    srand(time(0));
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
    catch(std::string error) {
        std::cout << error << std::endl;
        return false;
    }
    return true;
}

VertexVector GraphHelper::greedyMatching(Graph g) {

    VertexVector matching(num_vertices(g));

    Vertex null_vertex = g.null_vertex();

    // set all mates to null vector
    std::fill(matching.begin(), matching.end(), null_vertex);

    // do greedy matching over all edges
    EdgeIterator start, end;
    for (std::tie(start, end) = boost::edges(g); start != end; start++) {
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

void removeAdjacentEdges(Vertex v, Graph& g) {
    while (boost::degree(v, g) > 0) {
        auto startOutEdges = boost::out_edges(v, g).first;
        auto endOutEdges = boost::out_edges(v, g).second;
        for (auto o = startOutEdges; o != endOutEdges; ++o) {
            Vertex a = source(*o, g);
            Vertex b = target(*o, g);
            if (boost::edge(a, b, g).second) {
                boost::remove_edge(*o, g);
            }          
        }
    }
}

VertexVector GraphHelper::karpSipser(Graph g) {
    // Initialize matching with null
    VertexVector matching(num_vertices(g));
    std::fill(matching.begin(), matching.end(), g.null_vertex());

    while (num_edges(g) > 0) {
        // Search for edge with a vertex of degree 1
        bool foundEdge = false;

        for (auto e = boost::edges(g).first; e != boost::edges(g).second; ++e) {
            Vertex u = source(*e, g);
            Vertex v = target(*e, g);
            
            if(boost::degree(u, g) == 1 || boost::degree(v, g) == 1) {
                // Add vertices to matching
                matching[u] = v;
                matching[v] = u;
                
                // Remove out_edges from u,v
                removeAdjacentEdges(u, g);
                removeAdjacentEdges(v, g);

                foundEdge = true;
                break;
            }
        }

        if (foundEdge) continue;

        // If no edge found, select a random one
        int randEdge = rand() % num_edges(g);
        auto e = boost::edges(g).first;
        for (int i = 0; i < randEdge; i++)
            ++e;

        Vertex u = source(*e, g);
        Vertex v = target(*e, g);

        // Add vertices to matching
        matching[u] = v;
        matching[v] = u;
        
        // Remove out_edges from u,v
        removeAdjacentEdges(u, g);
        removeAdjacentEdges(v, g);
    }
    return matching;
}

void GraphHelper::printOutput(const std::string& algorithm, int numThreads, const std::vector<double>& durations) {
    std::cout << algorithm << "," << numThreads;
    for (double d : durations) 
        std::cout << "," << d;
    
    std::cout << std::endl;
}









