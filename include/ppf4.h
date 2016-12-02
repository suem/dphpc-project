
#pragma once

#include <atomic>
#include <vector>

#include "graphtypes.h"
#include "graphutils.h"

struct MateVisited {
    Vertex mate;
    std::atomic_size_t iteration;
};

void ppf4(const Graph& g, Vertex first_right, std::vector<MateVisited>& matching, int numThreads);
