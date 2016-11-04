//
// Created by suem on 10/27/16.
//
#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "graphtypes.h"
#include "pothen_fan.h"

void parallel_pothen_fan(const Graph& g, VertexVector& mate) {

	vertex_size_t  n = boost::num_vertices(g);
	VertexIterator start, end;

	// create initial greedy matching
	match_greedy(g, mate);

	Vertex null_vertex = g.null_vertex();

	std::atomic_flag* visited = new std::atomic_flag[n];

	//	std::vector<std::atomic_flag> visited(n);
	std::atomic_flag path_found = ATOMIC_FLAG_INIT;

	do {
		path_found.clear();
		for (vertex_size_t i = 0; i < n; i++) visited[i].clear();
		//        for (auto& flag : visited) flag.clear();

		std::tie(start, end) = boost::vertices(g);

		// we only do it in parallel if each thread has at least one vertex assigned.
		if (n > omp_get_num_threads()) {
#pragma omp parallel 
			{
				size_t nthreads = omp_get_num_threads();
				size_t ithread = omp_get_thread_num();

				int numberOfNodes = static_cast<int>(std::round(n / nthreads));
				VertexIterator startIt = start + ithread * numberOfNodes;
				VertexIterator endIt = ithread == nthreads - 1 ? end : startIt + numberOfNodes;

				for (; startIt != endIt; ++startIt) {
					const Vertex v = *startIt;
					if (is_right(v) || mate[v] != null_vertex) {
						// skip if vertex is on the right or if vertex is already matched
						continue;
					}
					// assert: v is on the left and unmatched
					bool path_found_v = find_path_atomic(v, g, mate, visited);
					if (path_found_v) path_found.test_and_set();
				}
			}
		}
		else {
			// do it sequential
			for (; start != end; ++start) {
				const Vertex v = *start;
				if (is_right(v) || mate[v] != null_vertex) {
					// skip if vertex is on the right or if vertex is already matched
					continue;
				}
				// assert: v is on the left and unmatched
				bool path_found_v = find_path_atomic(v, g, mate, visited);
				if (path_found_v) path_found.test_and_set();
			}
		}
	} while (path_found.test_and_set());

	delete[] visited;
}

void parallel_pothen_fan_opt(const Graph& g, VertexVector& mate) {

	size_t n = boost::num_vertices(g);
	VertexIterator start, end;

	// create initial greedy matching
	match_greedy(g, mate);

	const Vertex null_vertex = g.null_vertex();

	std::atomic_flag* visited = new std::atomic_flag[n];
	std::atomic_flag path_found = ATOMIC_FLAG_INIT;
	std::atomic<bool> path_found_check;

	do {
		path_found_check = false;
		path_found.clear();
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			visited[i].clear();
		}

		std::tie(start, end) = boost::vertices(g);

		const int nt = std::min(static_cast<int>(n), omp_get_num_threads());
#pragma omp parallel num_threads(nt)
		{
			size_t nthreads = omp_get_num_threads();
			size_t ithread = omp_get_thread_num();

			int numberOfNodes = static_cast<int>(std::round(n / nthreads));
			VertexIterator startIt = start + ithread * numberOfNodes;
			VertexIterator endIt = ithread == nthreads - 1 ? end : startIt + numberOfNodes;

			for (; startIt != endIt; ++startIt) {
				const Vertex v = *startIt;
				if (is_right(v) || mate[v] != null_vertex) {
					// skip if vertex is on the right or if vertex is already matched
					continue;
				}
				// assert: v is on the left and unmatched
				if (!path_found_check) {
					bool path_found_v = find_path_atomic(v, g, mate, visited);
					if (path_found_v) {
						path_found.test_and_set();
						path_found_check = true;
					}
				}
			}
		}

	} while (path_found.test_and_set());

	delete[] visited;
}

bool find_path_atomic(const Vertex x0, const Graph& g, VertexVector& mate, std::atomic_flag* visited) {

	std::stack<FindPathElement> stack;
	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, g);
	stack.push(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.top();
		if (vars.found) {
			mate[vars.y] = vars.x0;
			mate[vars.x0] = vars.y;

			// handle return value
			stack.pop();
			if (stack.empty()) {
				return true;
			}

			stack.top().found = true;
			continue;
		}

		bool leaveWhile = false;

		for (; vars.start != vars.end; ++vars.start) {
			vars.y = *vars.start;

			if (visited[vars.y].test_and_set()) {
				continue;
			}

			leaveWhile = true;

			if (mate[vars.y] == g.null_vertex()) { // y is unmatched
				mate[vars.y] = vars.x0;
				mate[vars.x0] = vars.y;

				// handle return value
				stack.pop();
				if (stack.empty()) {
					return true;
				}

				stack.top().found = true;
				break;
			}

			// y is matched with x1
			FindPathElement e2;
			e2.x0 = mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, g);
			stack.push(e2);

			break;
		}

		if (leaveWhile) {
			continue;
		}

		if (!stack.empty()) {
			stack.pop();
		}
	}

	return false;
}

void pothen_fan(const Graph& g, VertexVector& mate) {

	vertex_size_t  n = boost::num_vertices(g);

	// create initial greedy matching
	match_greedy(g, mate);

	Vertex null_vertex = g.null_vertex();

	bool path_found = true;
	bool* visited = new bool[n];

	while (path_found) {
		path_found = false;
		std::fill_n(visited, n, false);

		VertexIterator start, end;
		for (std::tie(start, end) = boost::vertices(g); start != end; ++start) {
			const Vertex v = *start;
			if (is_right(v) || mate[v] != null_vertex) {
				// skip if vertex is on the right or if vertex is already matched
				continue;
			}
			// assert: v is on the left and unmatched
			path_found = find_path(v, g, mate, visited) || path_found;
		}
	}

	delete[] visited;
}

bool find_path(const Vertex x0, const Graph& g, VertexVector& mate, bool* visited) {

	std::stack<FindPathElement> stack;
	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, g);
	stack.push(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.top();
		if (vars.found) {
			mate[vars.y] = vars.x0;
			mate[vars.x0] = vars.y;

			// handle return value
			stack.pop();
			if (stack.empty()) {
				return true;
			}

			stack.top().found = true;
			continue;
		}

		bool leaveWhile = false;

		for (; vars.start != vars.end; ++vars.start) {
			vars.y = *vars.start;

			if (visited[vars.y]) continue;
			visited[vars.y] = true;

			leaveWhile = true;

			if (mate[vars.y] == g.null_vertex()) { // y is unmatched
				mate[vars.y] = vars.x0;
				mate[vars.x0] = vars.y;

				// handle return value
				stack.pop();
				if (stack.empty()) {
					return true;
				}

				stack.top().found = true;
				break;
			}

			// y is matched with x1
			FindPathElement e2;
			e2.x0 = mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, g);
			stack.push(e2);

			break;
		}

		if (leaveWhile) {
			continue;
		}

		if (!stack.empty()) {
			stack.pop();
		}
	}

	return false;
}

bool find_path_recursive(const Vertex x0, const Graph& g, VertexVector& mate, bool* visited) {

	AdjVertexIterator start, end;

	for (std::tie(start, end) = boost::adjacent_vertices(x0, g); start != end; ++start) {
		Vertex y = *start;

		if (visited[y]) continue;
		visited[y] = true;

		if (mate[y] == g.null_vertex()) { // y is unmatched
			mate[y] = x0;
			mate[x0] = y;

			return true;
		}

		// y is matched with x1
		Vertex x1 = mate[y];

		bool path_found = find_path_recursive(x1, g, mate, visited);

		if (!path_found) {
			continue;
		}

		mate[y] = x0;
		mate[x0] = y;

		return true;
	}

	return false;
}

void match_greedy(const Graph& g, VertexVector& mate) {

	Vertex null_vertex = g.null_vertex();

	// set all mates to null vector
	std::fill(mate.begin(), mate.end(), null_vertex);

	// do greedy matching over all edges
	EdgeIterator start, end;
	for (std::tie(start, end) = boost::edges(g); start != end; start++) {
		Edge e = *start;
		Vertex u = boost::source(e, g);
		Vertex v = boost::target(e, g);
		if (mate[u] == null_vertex && mate[v] == null_vertex) {
			mate[u] = v;
			mate[v] = u;
		}
	}
}
