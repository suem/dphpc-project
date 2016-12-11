#pragma once

#include "graphtypes.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>

struct BenchmarkResult {
    // Graph information
    std::string graphName;
    size_t numEdges;
    size_t numVertices;

    // Benchmark information
    int iter;
    std::string timeStamp;
    std::string algorithm;
    std::vector<std::vector<double>> durations;
    std::vector<int> numThreads;
};

class GraphHelper {
    public:
        static Graph generateRandomGraph(int numNodes, float density, Vertex& first_right);
        static void readGraphFromFile(Graph& g, Vertex& first_right, unsigned long& matching_size_solution, const std::string& filePath);
        static void writeGraphToFile(const std::string& filePath, const Graph& g, const Vertex& first_right);
        static bool isMaximumMatching(const VertexVector& matching, const Graph& g); 
        static void printOutput(const BenchmarkResult& result, const std::string& outPath);
	static std::string getGraphNameFromPath(const std::string& graphPath);
        static float getAverageVertexDegree(const Graph& g);
};
