#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "tree_grafting.h"

// value taken from the paper
#define ALPHA 5.f

void ms_bfs_graft(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const vertex_size_t n = boost::num_vertices(g);
	const vertex_size_t n_right = n - first_right;
	const int nt = std::min(static_cast<int>(n), numThreads);

	volatile bool path_found;
	bool* visited = new bool[n];
	bool* augment_visited = new bool[n_right];

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
	
	// F <- all unmatched X vertices
	// for all those unmatched vertices, set the root to itself
	// parallel + lock vector?
	VertexVector F;
	for (Vertex x = 0; x < first_right; ++x) {
		if (is_unmatched(x, mate)) {
			F.push_back(x);
			root[x] = x;
		}
	}

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
				F = top_down_bfs(g, F, visited, parent, root, leaf, mate, nt);
			}
			else {
				// fill unvisited Y vertices
				VertexVector R;
				for (Vertex y = first_right; y < n; ++y) {
					if (!visited[y]) {
						R.push_back(y);
					}
				}
				F = bottom_up_bfs(g, R, visited, parent, root, leaf, mate, nt);
			}
		}

		memset(augment_visited, 0, sizeof(bool) * n_right);
		// step 2: augment matching. 
#pragma omp parallel num_threads(nt)
#pragma omp for
		for (int x = 0; x < first_right; ++x) {
			if (is_unmatched(x, mate)) {
				// if an augmenting path P from x is found then augment matching by P
				bool path_found_now = find_path_tg(x, g, first_right, mate, augment_visited);
				path_found = path_found_now || path_found;
			}
		}

		// step 3: construct frontier for the next phase
		if (path_found) {
			F = graft(g, first_right, visited, parent, root, leaf, mate, nt);
		}

	} while (path_found);

	delete[] visited;
	delete[] augment_visited;
}

VertexVector top_down_bfs(
	const Graph& g,
	VertexVector& F,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate,
	int numThreads) {

	// this should be a thread safe queue once.
	VertexVector queue;

	// TODO: in parallel
	for (int i = 0; i < F.size(); ++i) {
		Vertex x = F[i];
		if (is_active_tree(x, root, leaf)) {
			leaf[root[x]] = Graph::null_vertex();
			AdjVertexIterator start, end;
			for (std::tie(start, end) = boost::adjacent_vertices(x, g); start != end; ++start) {
				Vertex y = *start;
				if (!visited[y]) {
					updatePointers(x, y, visited, queue, parent, root, leaf, mate);
				}
			}
		}
	}

	return queue;
}

VertexVector bottom_up_bfs(
	const Graph& g,
	VertexVector& R,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate,
	int numThreads) {

	// this should be a thread safe queue once.
	VertexVector queue;

	// TODO: in parallel
	for (Vertex y : R) {
		AdjVertexIterator start, end;
		for (std::tie(start, end) = boost::adjacent_vertices(y, g); start != end; ++start) {
			Vertex x = *start;
			if (is_active_tree(x, root, leaf)) {
				leaf[root[x]] = Graph::null_vertex();
				updatePointers(x, y, visited, queue, parent, root, leaf, mate);
				break;
			}
		}
	}

	return queue;
}

VertexVector graft(
	const Graph& g,
	const Vertex first_right,
	bool* visited,
	VertexVector& parent,
	VertexVector& root,
	VertexVector& leaf,
	VertexVector& mate,
	int numThreads) {

	const vertex_size_t n = boost::num_vertices(g);
	VertexVector queue;

	VertexVector activeX;
	for (Vertex x = 0; x < first_right; ++x) {
		if (is_active_tree(x, root, leaf)) {
			activeX.push_back(x);
		}
	}

	VertexVector activeY;
	VertexVector renewableY;
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
		queue = bottom_up_bfs(g, renewableY, visited, parent, root, leaf, mate, numThreads);
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
	VertexVector& queue,
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

bool find_path_tg(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited) {

	// init stack
	std::vector<FindPathElement> stack;

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

			if (visited[vars.y - first_right]) continue;
			visited[vars.y - first_right] = true;

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
