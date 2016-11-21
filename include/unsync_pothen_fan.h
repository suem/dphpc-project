
#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

void unsync_parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);
