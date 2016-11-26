#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "tree_grafting.h"

// value taken from the paper
#define ALPHA 5.f

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const int NO_THREADS = numThreads; 
	const vertex_size_t n = boost::num_vertices(g);
	const int nt = std::min(static_cast<int>(n), NO_THREADS);

	volatile bool path_found;
	bool* visited = new bool[n];

	// initialize the root vector
	VertexVector root(n);
	std::fill(root.begin(), root.end(), Graph::null_vertex());
	// initialize the parent vector
	VertexVector parent(n);
	std::fill(parent.begin(), parent.end(), Graph::null_vertex());
	// initialize the leaf vector
	VertexVector leaf(n);
	std::fill(leaf.begin(), leaf.end(), Graph::null_vertex());
	// initialize the visited array
	memset(visited, 0, sizeof(bool) * n);
	
	// unmatchedX <- all unmatched X vertices
	// for all those unmatched vertices, set the root to itself
	std::deque<Vertex> F;
	for (Vertex x = 0; x < first_right; ++x) {
		if (is_unmatched(x, mate)) {
			F.push_back(x);
			root[x] = x;
		}
	}

	// init stack
	std::vector<FindPathElement> stack;

	do {
		path_found = false;
		
		// construct alternating BFS Forest
		while (!F.empty()) {
			// find numUnvisited
			size_t numUnvisitedY = 0; 
			for (Vertex y = first_right; y < n; ++y) {
				if (!visited[y]) {
					++numUnvisitedY;
				}
			}

			if (F.size() < numUnvisitedY / ALPHA) {
				F = top_down_bfs(g, F, visited, parent, root, leaf, mate);
			}
			else {
				// fill unvisited Y vertices
				std::deque<Vertex> R;
				for (Vertex y = first_right; y < n; ++y) {
					if (!visited[y]) {
						R.push_back(y);
					}
				}
				F = bottom_up_bfs(g, R, visited, parent, root, leaf, mate);
			}
		}

		// step 2: augment matching. This should be parallel.
		for (Vertex x = 0; x < first_right; ++x) {
			if (is_unmatched(x, mate)) {
				// if an augmenting path P from x is found then augment matching by P
				path_found = path_found || find_path_tg(x, g, first_right, mate, visited, stack);
			}
		}

		// step 3: construct frontier for the next phase
		F = graft(g, first_right, visited, parent, root, leaf, mate);

	} while (path_found);

	delete[] visited;
}

std::deque<Vertex> top_down_bfs(
	const Graph& g,
	std::deque<Vertex>& F,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate) {

	// this should be a thread safe queue once.
	std::deque<Vertex> queue;

	// TODO: in parallel
	for (Vertex x : F) {
		if (is_active_tree(x, root, leaf)) {
			leaf[root[x]] = Graph::null_vertex();
			AdjVertexIterator start, end;
			std::tie(start, end) = boost::adjacent_vertices(x, g);
			for (auto it = start; it != end; ++it) {
				if (!visited[*it]) {
					updatePointers(x, *it, visited, queue, parent, root, leaf, mate);
				}
			}
		}
	}

	return queue;
}

std::deque<Vertex> bottom_up_bfs(
	const Graph& g,
	std::deque<Vertex>& R,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate) {

	// this should be a thread safe queue once.
	std::deque<Vertex> queue;

	// TODO: in parallel
	for (Vertex y : R) {
		AdjVertexIterator start, end;
		std::tie(start, end) = boost::adjacent_vertices(y, g);
		for (auto it = start; it != end; ++it) {
			Vertex x = *it;
			if (is_active_tree(x, root, leaf)) {
				leaf[root[x]] = Graph::null_vertex();
				updatePointers(x, y, visited, queue, parent, root, leaf, mate);
				break;
			}
		}
	}

	return queue;
}

std::deque<Vertex> graft(
	const Graph& g,
	const Vertex first_right,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate) {

	const vertex_size_t n = boost::num_vertices(g);
	std::deque<Vertex> queue;

	std::deque<Vertex> activeX;
	for (Vertex x = 0; x < first_right; ++x) {
		if (is_active_tree(x, root, leaf)) {
			activeX.push_back(x);
		}
	}

	std::deque<Vertex> activeY;
	std::deque<Vertex> renewableY;
	for (Vertex y = first_right; y < n; ++y) {
		if (is_active_tree(y, root, leaf)) {
			activeY.push_back(y);
		}
		if (is_renewable_tree(y, root, leaf)) {
			renewableY.push_back(y);
		}
	}

	// TODO: in parallel
	for (auto y : renewableY) {
		visited[y] = false;
		root[y] = Graph::null_vertex();
	}

	if (activeX.size() > renewableY.size() / ALPHA) {
		queue = bottom_up_bfs(g, renewableY, visited, parent, root, leaf, mate);
	}
	else {
		// queue <- unmatched X vertices 
		for (Vertex x = first_right; x < n; ++x) {
			if (is_unmatched(x, mate)) {
				queue.push_back(x);
			}
		}

		// TODO: in parallel
		for (Vertex y : activeY) {
			visited[y] = false;
			root[y] = Graph::null_vertex();
		}
		
		// TODO: in parallel
		for (Vertex x : activeX) {
			// if x is in unmatchedX, exclude it
			if (is_unmatched(x, mate)) continue;
			root[x] = Graph::null_vertex();
		}
	}

	return queue;
}

void updatePointers(
	Vertex x,
	Vertex y,
	bool* visited,
	std::deque<Vertex>& queue,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate) {

	parent[y] = x;
	visited[y] = true;
	root[y] = root[x];
	if (is_matched(y, mate)) {
		queue.push_back(mate[y]);
		root[mate[y]] = root[y];
	} else {
		// an augmenting path is found, so we end the augmenting path
		leaf[root[x]] = y;
	}
}

bool find_path_tg(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::vector<FindPathElement>& stack) {

	stack.clear();

	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, g);
	stack.push_back(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.back();
		if (vars.found) {
			mate[vars.y] = vars.x0;
			mate[vars.x0] = vars.y;

			// handle return value
			stack.pop_back();
			if (stack.empty()) {
				return true;
			}

			stack.back().found = true;
			continue;
		}

		bool leaveWhile = false;

		for (; vars.start != vars.end; ++vars.start) {
			vars.y = *vars.start;

			if (visited[vars.y]) continue;
			visited[vars.y] = true;

			leaveWhile = true;

			if (is_unmatched(vars.y, mate)) { // y is unmatched
				mate[vars.y] = vars.x0;
				mate[vars.x0] = vars.y;

				// handle return value
				stack.pop_back();
				if (stack.empty()) {
					return true;
				}

				stack.back().found = true;
				break;
			}

			// y is matched with x1
			FindPathElement e2;
			e2.x0 = mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, g);
			stack.push_back(e2);

			break;
		}

		if (leaveWhile) {
			continue;
		}

		if (!stack.empty()) {
			stack.pop_back();
		}
	}

	return false;
}
