#include "karpsipser.h"

void pks_matchAndUpdate(
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

void parallel_karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching, const int numThreads) {
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

void ks_matchAndUpdate(
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

void karp_sipser(const Graph& g, const Vertex first_right, std::vector<Vertex>& matching) {
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
