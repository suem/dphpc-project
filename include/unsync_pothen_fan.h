
#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

void unsync_parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);

bool dfs_la(
		const Vertex x0,
		const Graph& g,const Vertex first_right, VertexVector& mate,
//		std::atomic_flag* visited,
		std::atomic<unsigned int>* visited, unsigned int iteration,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
		std::vector<PathElement>& stack);

