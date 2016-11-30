
#pragma once

#include <atomic>
#include "graphtypes.h"
#include "graphutils.h"

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);

VertexVector top_down_bfs(
	const Graph& g, 
	VertexVector& F,
	bool* visited,
	VertexVector& parent, 
	VertexVector& root, 
	VertexVector& leaf, 
	VertexVector& mate,
	int numThreads);

VertexVector bottom_up_bfs(
	const Graph& g,
	VertexVector& R,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate,
	int numThreads);

// rebuild frontier (queue) for the next phase
VertexVector graft(
	const Graph& g,
	const Vertex first_right,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate,
	int numThreads);

// update pointers in BSF traversals
void updatePointers(
	Vertex x,
	Vertex y,
	bool* visited,
	VertexVector& queue,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate);

// whether the Vertex v is in an active tree
inline bool is_active_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	return root[v] != Graph::null_vertex() && leaf[root[v]] == Graph::null_vertex();
}

// whether the Vertex v is in a renewable tree
inline bool is_renewable_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	return root[v] != Graph::null_vertex() && leaf[root[v]] != Graph::null_vertex();
}

bool find_path_tg(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited);
