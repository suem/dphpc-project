#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "tree_grafting.h"

// value taken from the paper
const float TreeGrafting::ALPHA = 5.f;

// whether the Vertex v is in an active tree
inline bool is_active_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	return root[v] != Graph::null_vertex() && leaf[root[v]] == Graph::null_vertex();
}

// whether the Vertex v is in a renewable tree
inline bool is_renewable_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	return root[v] != Graph::null_vertex() && leaf[root[v]] != Graph::null_vertex();
}

TreeGrafting::TreeGrafting(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) :
	m_graph(g), m_mate(mate), m_firstRight(first_right) {
	n = boost::num_vertices(g);
	n_right = n - first_right;
	m_numThreads = std::min(static_cast<int>(n), numThreads);

	init();
	ms_bfs_graft();
}

TreeGrafting::~TreeGrafting() {
	delete[] m_visited;
	delete[] m_augmentVisited;
}

void TreeGrafting::init() {
	m_visited = new bool[n];
	m_augmentVisited = new std::atomic_flag[n_right];

	// initialize the root vector
	m_root.resize(n);
	std::fill(m_root.begin(), m_root.end(), Graph::null_vertex());
	// initialize the parent vector
	m_parent.resize(n);
	std::fill(m_parent.begin(), m_parent.end(), Graph::null_vertex());
	// initialize the leaf vector
	m_leaf.resize(n);
	std::fill(m_leaf.begin(), m_leaf.end(), Graph::null_vertex());
	// initialize the visited array
	memset(m_visited, 0, sizeof(bool) * n);
}

void TreeGrafting::ms_bfs_graft() {
	volatile bool path_found;

	// F <- all unmatched X vertices
	// for all those unmatched vertices, set the root to itself
	// parallel + lock vector?
	VertexVector F;
	for (Vertex x = 0; x < m_firstRight; ++x) {
		if (is_unmatched(x, m_mate)) {
			F.push_back(x);
			m_root[x] = x;
		}
	}

	do {
		path_found = false;

		// construct alternating BFS Forest
		while (!F.empty()) {
			// find numUnvisited
			size_t numUnvisitedY = 0;
			for (Vertex y = m_firstRight; y < n; ++y) {
				if (!m_visited[y]) {
					++numUnvisitedY;
				}
			}

			if (F.size() < numUnvisitedY / ALPHA) {
				F = top_down_bfs(F);
			}
			else {
				// fill unvisited Y vertices
				VertexVector R;
				for (Vertex y = m_firstRight; y < n; ++y) {
					if (!m_visited[y]) {
						R.push_back(y);
					}
				}
				F = bottom_up_bfs(R);
			}
		}

		memset(m_augmentVisited, 0, sizeof(std::atomic_flag) * n_right);
		// step 2: augment matching. 
#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
		for (int x = 0; x < m_firstRight; ++x) {
			if (is_unmatched(x, m_mate)) {
				// if an augmenting path P from x is found then augment matching by P
				path_found = find_path_tg(x) || path_found;
			}
		}

		// step 3: construct frontier for the next phase
		if (path_found) {
			F = graft();
		}

	} while (path_found);
}

VertexVector TreeGrafting::top_down_bfs(VertexVector& F) {
	// this should be a thread safe queue once.
	VertexVector queue;

	// TODO: in parallel
	for (int i = 0; i < F.size(); ++i) {
		Vertex x = F[i];
		if (is_active_tree(x, m_root, m_leaf)) {
			m_leaf[m_root[x]] = Graph::null_vertex();
			AdjVertexIterator start, end;
			for (std::tie(start, end) = boost::adjacent_vertices(x, m_graph); start != end; ++start) {
				Vertex y = *start;
				if (!m_visited[y]) {
					updatePointers(x, y, queue);
				}
			}
		}
	}

	return queue;
}

VertexVector TreeGrafting::bottom_up_bfs(VertexVector& R) {
	// this should be a thread safe queue once.
	VertexVector queue;

	// TODO: in parallel
	for (Vertex y : R) {
		AdjVertexIterator start, end;
		for (std::tie(start, end) = boost::adjacent_vertices(y, m_graph); start != end; ++start) {
			Vertex x = *start;
			if (is_active_tree(x, m_root, m_leaf)) {
				m_leaf[m_root[x]] = Graph::null_vertex();
				updatePointers(x, y, queue);
				break;
			}
		}
	}

	return queue;
}

VertexVector TreeGrafting::graft() {
	VertexVector queue;

	VertexVector activeX;
	for (Vertex x = 0; x < m_firstRight; ++x) {
		if (is_active_tree(x, m_root, m_leaf)) {
			activeX.push_back(x);
		}
	}

	VertexVector activeY;
	VertexVector renewableY;
	for (Vertex y = m_firstRight; y < n; ++y) {
		if (is_active_tree(y, m_root, m_leaf)) {
			activeY.push_back(y);
		}
		if (is_renewable_tree(y, m_root, m_leaf)) {
			renewableY.push_back(y);
		}
	}

	// TODO: in parallel
	for (auto y : renewableY) {
		m_visited[y] = false;
		m_root[y] = Graph::null_vertex();
	}

	if (activeX.size() > renewableY.size() / ALPHA) {
		queue = bottom_up_bfs(renewableY);
	}
	else {
		// queue <- unmatched X vertices 
		for (Vertex x = m_firstRight; x < n; ++x) {
			if (is_unmatched(x, m_mate)) {
				queue.push_back(x);
			}
		}

		// TODO: in parallel
		for (Vertex y : activeY) {
			m_visited[y] = false;
			m_root[y] = Graph::null_vertex();
		}

		// TODO: in parallel
		for (Vertex x : activeX) {
			// if x is in unmatchedX, exclude it
			if (is_unmatched(x, m_mate)) continue;
			m_root[x] = Graph::null_vertex();
		}
	}

	return queue;
}

void TreeGrafting::updatePointers(const Vertex x, const Vertex y, VertexVector& queue) {
	m_parent[y] = x;
	m_visited[y] = true;
	m_root[y] = m_root[x];
	if (is_matched(y, m_mate)) {
		queue.push_back(m_mate[y]);
		m_root[m_mate[y]] = m_root[y];
	}
	else {
		// an augmenting path is found, so we end the augmenting path
		m_leaf[m_root[x]] = y;
	}
}

bool TreeGrafting::find_path_tg(const Vertex x0) {
	// init stack
	std::vector<FindPathElement> stack;

	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, m_graph);
	stack.push_back(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.back();
		if (vars.found) {
			m_mate[vars.y] = vars.x0;
			m_mate[vars.x0] = vars.y;

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

			// skip vertex if this was alredy visited by another thread
			if (m_augmentVisited[vars.y - m_firstRight].test_and_set()) continue;

			leaveWhile = true;

			if (is_unmatched(vars.y, m_mate)) { // y is unmatched
				m_mate[vars.y] = vars.x0;
				m_mate[vars.x0] = vars.y;

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
			e2.x0 = m_mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, m_graph);
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
