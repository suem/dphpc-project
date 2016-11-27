
#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

void ppf2(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);
