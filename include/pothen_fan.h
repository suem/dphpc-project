//
// Created by suem on 10/27/16.
//

#pragma once

#include "graphtypes.h"

struct FindPathElement {
	Vertex x0;
	Vertex y;
	AdjVertexIterator start;
	AdjVertexIterator end;
	bool found = false;
};

void pothen_fan(const Graph& g, VertexVector& mate);
void match_one(const Graph& g, VertexVector& mate);
void match_greedy(const Graph& g, VertexVector& mate);
bool find_path(const Vertex x0, const Graph& g, VertexVector& mate, bool* visited);
bool find_path_recursive(const Vertex x0, const Graph& g, VertexVector& mate, bool* visited);

inline bool is_right(const Vertex& v) {
	return v % 2 == 1;
}

inline bool is_left(const Vertex& v) {
	return v % 2 == 0;
}
