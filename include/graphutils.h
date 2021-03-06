// Collection of graph matching related utility functions

#pragma once

#include "graphtypes.h"
#include <atomic>

struct PathElement {
	Vertex x0;
	AdjVertexIterator yiter;
	AdjVertexIterator yiter_end;
};

struct FindPathElement {
	Vertex x0;
	Vertex y;
	AdjVertexIterator start;
	AdjVertexIterator end;
	bool found = false;
};

struct UnmatchedVertex {
    UnmatchedVertex(Vertex _x) : x(_x) { }
    Vertex x;
    bool unmatched = true;
};

inline bool is_right(const Vertex& v, const Vertex first_right) {
    return v >= first_right;
}

inline bool is_left(const Vertex& v, const Vertex first_right) {
    return v < first_right;
}

inline bool is_matched(const Vertex& v, const VertexVector& mate) {
    return mate[v] != Graph::null_vertex();
}

inline bool is_unmatched(const Vertex& v, const VertexVector& mate) {
    return mate[v] == Graph::null_vertex();
}

inline bool claim_vertex(
        Vertex y, Vertex first_right,
        std::atomic<unsigned int>* visited, unsigned int iteration) {
    return std::atomic_exchange(&visited[y - first_right], iteration) != iteration;
}
