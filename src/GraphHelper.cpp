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

Graph GraphHelper::readGraphFromFile(const std::string& filePath) {
    std::ifstream inFile(filePath);
    if (inFile.fail())
        std::cerr << "Error while reading from file" << std::endl;
    
    int numNodes, u, v;

    inFile >> numNodes;
    Graph g(numNodes);

    while (inFile >> u >> v) 
        boost::add_edge(u, v, g);
       
    return g;
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
		verify_matching(g, matching);
	}
	catch(std::string error) {
		std::cout << error << std::endl;
		return false;
	}
	return true;
}
