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
    //int randNodes = rand() % numNodes + 1;
    int randNodes = numNodes / 2;
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

void GraphHelper::readGraphFromFile(Graph& g, Vertex& first_right, unsigned long& matching_size_soluton, const std::string& filePath) {
    std::ifstream inFile(filePath);
    if (inFile.fail()) std::cerr << "Error while reading from file" << std::endl;

    int numNodes, u, v;

    inFile >> numNodes;
    inFile >> first_right;
	inFile >> matching_size_soluton;
    g = Graph(numNodes);

    while (inFile >> u >> v) boost::add_edge(u, v, g);
}

void GraphHelper::writeGraphToFile(const std::string& filePath, const Graph& g, const Vertex& first_right) {
    std::ofstream outFile(filePath);
    outFile << num_vertices(g) << std::endl;
    outFile << first_right << std::endl;

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
        unsigned long matching_size_solution = boost::matching_size(g, &solution_mates[0]);
        verify_matching(g, matching, matching_size_solution);
    }
    catch (std::string error) {
        std::cout << error << std::endl;
        return false;
    }
    return true;
}

void GraphHelper::printOutput(const BenchmarkResult& resultStruct, const std::string& outPath) {
	// generate file name: algorithm name, graphName, and timestamp
	std::string fileName = resultStruct.algorithm + '_' + resultStruct.graphName + '_' + resultStruct.timeStamp + ".csv";
	std::string filePath = outPath + fileName;
	std::replace(filePath.begin(), filePath.end(), ' ', '_');
	// create csv file
	std::ofstream csv(filePath);


    // Header
    std::cout << "TimeStamp" << ",";
    std::cout << "GraphName" << ",";
    std::cout << "NumVertices" << ",";
    std::cout << "NumEdges" << ",";
    std::cout << "Algorithm" << std::endl;

	std::cout << resultStruct.timeStamp << ",";
	std::cout << resultStruct.graphName << ",";
	std::cout << resultStruct.numVertices << ",";
	std::cout << resultStruct.numEdges << ",";
	std::cout << resultStruct.algorithm << std::endl;

	std::cout << "Duration,NumThreads,Algorithm" << std::endl;
    
	/*New format
	for (int nOfThreads : resultStruct.numThreads) {
		for (int i : resultStruct.durations[nOfThreads]) {
			std::cout << i << nOfThreads << resultStruct.algorithm << std::endl;
		}
	}
    */

	/* Old format */
    // Data
	std::string temp;
	for (int i : resultStruct.numThreads)
		temp += std::to_string(i) + ",";

	std::cout << temp.substr(0, temp.size() - 1) << std::endl;

    for (int i = 0; i < resultStruct.iter; ++i) {
		// loop over iterations
        int pos = 0;
		temp.clear();
        for (int nThreads : resultStruct.numThreads) {
			// loop over number of threads
			temp += std::to_string(resultStruct.durations[pos][i]) + ",";
            ++pos;
        }
		std::cout << temp.substr(0, temp.size() - 1) << std::endl;
    } 
	//std::cout << "CSV Filepath: " << filePath<< std::endl;
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

float GraphHelper::getAverageVertexDegree(const Graph& g) {
    float result = 0.0;
    for (auto v = boost::vertices(g).first; v != boost::vertices(g).second; ++v) {
        result += boost::degree(*v, g);
    }
    return result / boost::num_vertices(g);    
}
