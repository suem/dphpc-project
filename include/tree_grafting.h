
#pragma once

#include <atomic>
#include <queue>
#include "graphtypes.h"

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);

std::queue<Vertex> top_down_bfs(
	const Graph& g, 
	std::set<Vertex>& unmatchedX, 
	bool* visited, 
	VertexVector& parent, 
	VertexVector& root, 
	VertexVector& leaf, 
	VertexVector& mate);

// update pointers in BSF traversals
void updatePointers(
	Vertex x,
	Vertex y,
	std::queue<Vertex>& queue,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate);

// whether the Vertex v is in an active tree
inline bool is_active_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	root[v] != Graph::null_vertex() && leaf[root[v]] == Graph::null_vertex();
}