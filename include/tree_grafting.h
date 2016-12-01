
#pragma once

#include <atomic>
#include "graphtypes.h"
#include "graphutils.h"

class TreeGrafting {
public:
	TreeGrafting(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);
	~TreeGrafting();

private:
	// initializes member arrays
	void init();
	// execute the tree grafting
	void ms_bfs_graft();

	// breadth first search from above
	VertexVector top_down_bfs(VertexVector& F);
	// breadth first search from below
	VertexVector bottom_up_bfs(VertexVector& R);
	// rebuild frontier (queue) for the next phase
	VertexVector graft();
	// update pointers in bfs traversals
	void updatePointers(const Vertex x, const Vertex y, VertexVector& queue);
	// find and augment a path
	bool find_path_tg(const Vertex x0);

private:
	int m_numThreads;

	const Graph& m_graph;
	VertexVector& m_mate;
	const Vertex m_firstRight;
	Vertex n_right;
	size_t n;

	VertexVector m_parent;
	VertexVector m_leaf;
	VertexVector m_root;

	bool* m_visited = nullptr;
	std::atomic_flag* m_augmentVisited = nullptr;

	static const float ALPHA;
};

