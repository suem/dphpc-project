#pragma once

#include "graphtypes.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>

class GraphHelper {
    public:
        static Graph generateRandomGraph(int numNodes, float density);
        static void readGraphFromFile(Graph& g, Vertex& first_right, const std::string& filePath);
        static void writeGraphToFile(const std::string& filePath, const Graph& g);
        static bool isMaximumMatching(const VertexVector& matching, const Graph& g); 
};
