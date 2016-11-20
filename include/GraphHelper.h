#pragma once

#include "graphtypes.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>

struct BenchmarkResult {
    // Graph information
    std::string graphName;
    size_t numEdges;
    size_t numVertices;

    // Benchmark information
    std::string timeStamp;
    std::string algorithm;
    std::vector<double> durations;
    int numThreads;
};

class GraphHelper {
    public:
        static Graph generateRandomGraph(int numNodes, float density);
        static void readGraphFromFile(Graph& g, Vertex& first_right, const std::string& filePath);
        static void writeGraphToFile(const std::string& filePath, const Graph& g);
        static bool isMaximumMatching(const VertexVector& matching, const Graph& g); 
		// karp sipser using samuels matching method 
        static VertexVector karpSipser(const Graph& g);
		// kind of karp sipser very fast but no break
        static VertexVector ks(const Graph& g);
        static VertexVector greedyMatching(const Graph& g);
        static void printOutput(const BenchmarkResult& result);
};
