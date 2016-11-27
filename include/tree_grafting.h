
#pragma once

#include <atomic>
#include <queue>
#include "graphtypes.h"
#include "graphutils.h"

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);

std::deque<Vertex> top_down_bfs(
	const Graph& g, 
	std::deque<Vertex>& F,
	bool* visited,
	VertexVector& parent, 
	VertexVector& root, 
	VertexVector& leaf, 
	VertexVector& mate);

std::deque<Vertex> bottom_up_bfs(
	const Graph& g,
	std::deque<Vertex>& R,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate);

// rebuild frontier (queue) for the next phase
std::deque<Vertex> graft(
	const Graph& g,
	const Vertex first_right,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate);

// update pointers in BSF traversals
void updatePointers(
	Vertex x,
	Vertex y,
	bool* visited,
	std::deque<Vertex>& queue,
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
