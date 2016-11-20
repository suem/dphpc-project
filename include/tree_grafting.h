
#pragma once

#include <atomic>
#include "graphtypes.h"

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);

