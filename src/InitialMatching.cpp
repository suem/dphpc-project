#include "InitialMatching.h"

void InitialMatching::pks_matchAndUpdate(
        Vertex u,
        const Graph& g,
        VertexVector& matching,
        std::atomic_int* deg,
        std::atomic_flag* visited) {

    if (visited[u].test_and_set()) return;

    AdjVertexIterator start, end;
    std::tie(start, end) = boost::adjacent_vertices(u, g);
    Vertex v;
    while (start != end) {
        v = *start;
        if (!visited[v].test_and_set()) {
            break;
        }
        start++;
    }

    if (start == end) {
        return;
    }

    // add (u,v) to matching
    matching[u] = v;
    matching[v] = u;

    for (std::tie(start, end) = boost::adjacent_vertices(v, g); start != end; start++) {
        Vertex w = *start;
        if (std::atomic_fetch_add(&deg[w], -1) == 2) {
            pks_matchAndUpdate(w, g, matching, deg, visited);
        }
    }
}

void InitialMatching::parallel_karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching, const int numThreads) {
    vertex_size_t n = boost::num_vertices(g);
    matching.assign(n, Graph::null_vertex());

    std::atomic_flag* visited = new std::atomic_flag[n];
    memset(visited, 0, sizeof(std::atomic_flag) * n);

    std::atomic_int* deg = new std::atomic_int[first_right];

#pragma omp parallel num_threads(numThreads)
#pragma omp for
    for (int x = 0; x < first_right; x++)
        deg[x] = static_cast<int>(boost::degree(x, g));

#pragma omp parallel num_threads(numThreads)
#pragma omp for
    for (int x = 0; x < first_right; x++)
        if (deg[x] == 1) pks_matchAndUpdate(x, g, matching, deg, visited);

#pragma omp parallel num_threads(numThreads)
#pragma omp for
    for (int x = 0; x < first_right; x++)
        if (deg[x] > 1) pks_matchAndUpdate(x, g, matching, deg, visited);

}

void InitialMatching::ks_matchAndUpdate(
        const Vertex u,
        const Graph& g,
        VertexVector& matching,
        std::vector<size_t>& deg,
        std::vector<bool>& visited) {

    if (visited[u]) return;
    visited[u] = true;

    AdjVertexIterator start, end;
    std::tie(start, end) = boost::adjacent_vertices(u, g);
    Vertex v;
    while (start != end) {
        v = *start;

        if (!visited[v]) {
            visited[v] = true;
            break;
        }

        start++;
    }

    if (start == end) {
        return;
    }

    // add (u,v) to matching
    matching[u] = v;
    matching[v] = u;

    for (std::tie(start, end) = boost::adjacent_vertices(v, g); start != end; start++) {
        Vertex w = *start;
        deg[w] -= 1;
        if (deg[w] == 1) {
            ks_matchAndUpdate(w, g, matching, deg, visited);
        }
    }
}

void InitialMatching::karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching) {
    vertex_size_t n = boost::num_vertices(g);
    matching.assign(n, Graph::null_vertex());

    std::vector<bool> visited(n, false);
    std::vector<size_t> deg(first_right);

    for (int x = 0; x < first_right; x++)
        deg[x] = boost::degree(x, g);

    for (int x = 0; x < first_right; x++)
        if (deg[x] == 1) ks_matchAndUpdate(x, g, matching, deg, visited);

    for (int x = 0; x < first_right; x++)
        if (deg[x] > 1) ks_matchAndUpdate(x, g, matching, deg, visited);
}


void InitialMatching::greedy_matching(const Graph& g, VertexVector& matching) {
	// set all mates to null vector
	std::fill(matching.begin(), matching.end(), Graph::null_vertex());

	// do greedy matching over all edges
	EdgeIterator start, end;
	for (std::tie(start, end) = boost::edges(g); start != end; ++start) {
		Edge e = *start;
		Vertex u = boost::source(e, g);
		Vertex v = boost::target(e, g);
		if (matching[u] == Graph::null_vertex() && matching[v] == Graph::null_vertex()) {
			matching[u] = v;
			matching[v] = u;
		}
	}
}

void InitialMatching::markAdjacentEdges(Vertex v, const Graph& g, std::vector<size_t>& degree) {
	auto startOutEdges = boost::out_edges(v, g).first;
	auto endOutEdges = boost::out_edges(v, g).second;
	degree[v] = 0;
	for (auto o = startOutEdges; o != endOutEdges; ++o) {
		Vertex b = target(*o, g);
		--degree[b];
	}
}

void InitialMatching::karp_sipser_fast(const Graph& g, VertexVector& matching) {
	auto const n = boost::num_vertices(g);

	std::vector<Vertex> deg(n);
	for (Vertex v = 0; v < n; v++) {
		matching[v] = Graph::null_vertex();
		deg[v] = boost::degree(v, g);
	}

	bool found_edges;

	do {
		found_edges = false;

		for (auto e = boost::edges(g).first; e != boost::edges(g).second; ++e) {
			Vertex u = boost::source(*e, g);
			Vertex v = boost::target(*e, g);

			if (is_matched(u, matching) || is_matched(v, matching)) {
				continue;
			}

			if (deg[v] == 1 || deg[u] == 1) {
				matching[v] = u;
				matching[u] = v;
				markAdjacentEdges(v, g, deg);
				markAdjacentEdges(u, g, deg);
				found_edges = true;
				break;
			}
		}

		if (found_edges) continue;

		// If no edge found, select a random unmatched one
		EdgeIterator eStart, eEnd, eRandom, eIt;
		std::tie(eStart, eEnd) = boost::edges(g);
		int randEdge = rand() % num_edges(g);
		eRandom = eStart;
		for (int i = 0; i < randEdge; i++)
			++eRandom;

		eIt = eRandom;
		do {
			Edge e = *eIt;
			Vertex u = boost::source(e, g);
			Vertex v = boost::target(e, g);

			if (!is_matched(u, matching) && !is_matched(v, matching)) {
				matching[v] = u;
				matching[u] = v;
				markAdjacentEdges(v, g, deg);
				markAdjacentEdges(u, g, deg);
				found_edges = true;
				break;
			}

			++eIt;
			if (eIt == eEnd) {
				eIt = eStart;
			}

		} while (eIt != eRandom);

	} while (found_edges);
}

void InitialMatching::ks(const Graph& g, VertexVector& initialMatching) {
	auto const n = boost::num_vertices(g);

	VertexVector matching(n);
	std::vector<Vertex> deg(n);
	for (Vertex v = 0; v < n; v++) {
		matching[v] = Graph::null_vertex();
		deg[v] = boost::degree(v, g);
	}

	bool found_edges;
	do {
		// match all vertices with degree one
		bool found_deg_one;
		do {
			found_deg_one = false;

			for (Vertex v = 0; v < n; v++) {
				if (deg[v] != 1 || is_matched(v, matching)) continue;

				AdjVertexIterator s, e;
				for (std::tie(s, e) = boost::adjacent_vertices(v, g); s != e; s++) {
					Vertex u = *s;
					if (is_matched(u, matching)) continue;

					matching[v] = u;
					matching[u] = v;
					deg[u]--;
					deg[v]--;
					found_deg_one = true;
				}
			}
		} while (found_deg_one);

		found_edges = false;
		// take first unmatched edge
		EdgeIterator estart, eend;
		for (std::tie(estart, eend) = boost::edges(g); estart != eend; estart++) {
			Edge e = *estart;
			Vertex u = boost::source(e, g);
			Vertex v = boost::target(e, g);

			if (is_matched(u, matching) || is_matched(v, matching)) continue;

			matching[v] = u;
			matching[u] = v;
			deg[u]--;
			deg[v]--;
			found_edges = true;
		}

	} while (found_edges);
}

