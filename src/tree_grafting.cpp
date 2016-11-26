#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>
#include <set>

#include "graphtypes.h"
#include "tree_grafting.h"
#include "graphutils.h"

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const int NO_THREADS = numThreads; 

	const vertex_size_t n = boost::num_vertices(g);
	const vertex_size_t n_right = n - first_right;
	const Vertex null_vertex = Graph::null_vertex();

	const int nt = std::min(static_cast<int>(n), NO_THREADS);

	volatile bool path_found;

	bool* visited = new bool[n_right];

	// initialize the root vector
	VertexVector root(n);
	std::fill(root.begin(), root.end(), null_vertex);
	// initialize the parent vector
	VertexVector parent(n_right);
	std::fill(parent.begin(), parent.end(), null_vertex);
	// initialize the leaf vector
	VertexVector leaf(first_right);
	std::fill(leaf.begin(), leaf.end(), null_vertex);
	// initialize the visited array
	memset(visited, 0, sizeof(bool) * n_right);

	// unmatchedX <- all unmatched X vertices
	// for all those unmatched vertices, set the root to itself
	std::set<Vertex> unmatchedX;
	for (Vertex v = first_right; v < n; ++v) {
		if (is_unmatched(v, mate)) {
			unmatchedX.insert(v);
			root[v] = v;
		}
	}

	// initialize numUnvisited and alpha
	size_t numUnvisitedY = first_right; // all in Y are unvisited
	float alpha = 5.f; // value taken from the paper

	do {
		path_found = false;

		// construct alternating BFS Forest
		while (!unmatchedX.empty()) {
			if (unmatchedX.size() < numUnvisitedY / alpha) {
				// TODO
				// top_down(g, unmachedX, parent, root, leaf, mate); 
			}
			else {
				// fill unvisited Y vertices
				std::set<Vertex> unvisitedY;
				for (Vertex v = 0; v < first_right; ++v) {
					if (!visited[v]) {
						unvisitedY.insert(v);
					}
				}
				// bottomUp(g, unvisitedY, visited, parent, root, leaf, mate);
			}
		}

		// step 2: augment matching
		for (Vertex v = first_right; v < n; ++v) {
			if (is_unmatched(v, mate)) {
				// if an augmenting path P from x0 is found then
				// Augment matching by P
			}
		}

		// step 3: construct frontier for the next phase
		// unmatchedX = graft(g, visited, parent, root, leaf, mate);

	} while (path_found);

	delete[] visited;
}

std::queue<Vertex> top_down_bfs(
	const Graph& g,
	std::set<Vertex>& unmatchedX,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate) {

	// this should be a thread safe queue once.
	std::queue<Vertex> queue;

	// here, we will use a parallel for in the parallel version
	for (Vertex v = 0; v < unmatchedX.size(); ++v) {
		if (is_active_tree(v, root, leaf)) {
			leaf[root[v]] = Graph::null_vertex();
		}
	}

	return queue;
}

void updatePointers(
	Vertex x,
	Vertex y,
	bool* visited,
	std::queue<Vertex>& queue,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate) {
	parent[y] = x;
	visited[y] = true;
	root[y] = root[x];
	if (is_matched(y, mate)) {
		queue.push(mate[y]);
		root[mate[y]] = root[y];
	} else {
		// an augmenting path is found, so we end the augmenting path
		leaf[root[x]] = y;
	}
}

