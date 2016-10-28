#include "GraphHelper.h"

Graph GraphHelper::generateRandomGraph(int numNodes, float density) { 
    Graph g(numNodes);
    int randNodes = rand() % numNodes + 1;
    int numEdges = randNodes * (numNodes - randNodes) * density;
    
    if (numEdges == 0)
        std::cerr << "Warning! Graph has no edges." << std::endl;
    
    for (int i = 0; i < numEdges; i++) {
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
    int numNodes, u, v;

    inFile >> numNodes;
    Graph g(numNodes);

    while (inFile >> u >> v) 
        boost::add_edge(u, v, g);
    
    if (inFile.fail())
        std::cerr << "Error while reading from file" << std::endl;

    return g;
}

void GraphHelper::writeGraphToFile(const std::string& filePath, const Graph& g) {
    std::ofstream outFile(filePath);
    outFile << num_vertices(g) << std::endl;
    
    for (auto e = edges(g).first; e != edges(g).second; e++)
        outFile << source((*e).first, g) << " " << target((*e).first, g) << std::endl;
    
    outFile.flush();
    if (outFile.fail())
        std::cerr << "Error while writing to file." << std::endl;        
}
