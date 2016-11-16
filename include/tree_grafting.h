
#pragma once

#include <atomic>
#include "graphtypes.h"

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);

inline bool is_right(const Vertex& v, const Vertex first_right) {
	return v >= first_right;
}

inline bool is_left(const Vertex& v, const Vertex first_right) {
	return v < first_right;
}

inline bool is_matched(const Vertex& v, const Graph& g, const VertexVector& mate) {
	return mate[v] != g.null_vertex();
}

inline bool is_unmatched(const Vertex& v, const Graph& g, const VertexVector& mate) {
	return mate[v] == g.null_vertex();
}

