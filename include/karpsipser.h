#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

void karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching);

void parallel_karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching, const int numThreads);
