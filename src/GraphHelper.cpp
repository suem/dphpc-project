#include "GraphHelper.h"
#include "verifier.h"

Graph GraphHelper::generateRandomGraph(int numNodes, float density) { 
    Graph g(numNodes);
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

VertexVector GraphHelper::karpSipser(Graph g) {
    VertexVector matching(num_vertices(g));

    while (num_edges(g) > 0) {
        // Search for edge with a vertex of degree 1
        bool foundEdge = false;

        auto startEdges = boost::edges(g).first;
        auto endEdges = boost::edges(g).second;
        for (auto e = startEdges; e != endEdges; ++e) {
            if(boost::degree(source(*e, g), g) == 1 || boost::degree(target(*e, g), g) == 1) {
                matching[source(*e, g)] = target(*e, g);
                matching[target(*e, g)] = source(*e, g);

                // Remove out_edges from u,v
                auto startOutEdges = boost::out_edges(source(*e, g), g).first;
                auto endOutEdges = boost::out_edges(source(*e, g), g).second;
                for (auto o = startOutEdges; o != endOutEdges; ++o) 
                    boost::remove_edge(*o, g);                

                startOutEdges = boost::out_edges(target(*e, g), g).first;
                endOutEdges = boost::out_edges(target(*e, g), g).second;
                for (auto o = startOutEdges; o != endOutEdges; ++o) 
                    boost::remove_edge(*o, g);

                foundEdge = true;
                break;
            }
        }

        if (foundEdge) continue;

        // If no edge found, select the first one
        auto e = boost::edges(g).first;
        matching[source(*e, g)] = target(*e, g);
        matching[target(*e, g)] = source(*e, g);

        auto startOutEdges = boost::out_edges(source(*e, g), g).first;
        auto endOutEdges = boost::out_edges(source(*e, g), g).second;
        for (auto o = startOutEdges; o != endOutEdges; ++o) 
            boost::remove_edge(*o, g);                
    }
    return matching;
}
