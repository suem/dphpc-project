#pragma once

#include "graphtypes.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>

typedef struct benchmarkResult {
    // Graph information
    std::string graphName;
    int numEdges;
    int numVertices;

    // Benchmark information
    std::string timeStamp;
    std::string algorithm;
    std::vector<double> durations;
    int numThreads;
} benchmarkResult;

class GraphHelper {
    public:
        static Graph generateRandomGraph(int numNodes, float density);
        static void readGraphFromFile(Graph& g, Vertex& first_right, const std::string& filePath);
        static void writeGraphToFile(const std::string& filePath, const Graph& g);
        static bool isMaximumMatching(const VertexVector& matching, const Graph& g); 
        static VertexVector karpSipser(Graph g);
        static VertexVector greedyMatching(Graph g);
        static void printOutput(const benchmarkResult& result);
};
