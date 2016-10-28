//
// Created by suem on 10/27/16.
//

#pragma once

#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

inline void verify_bipartite(const Graph& g) {
    bool bp = boost::is_bipartite(g);
    if (!bp) throw "Graph not bipartite";
}

inline void verify_matching(const Graph& g, const VertexVector& mate) {
    vertex_size_t n = num_vertices(g);
    VertexVector solution(n);
    bool success = boost::checked_edmonds_maximum_cardinality_matching(g, &solution[0]);
    if (!success) throw "Graph not bipartite";

    vertex_size_t matching_size = boost::matching_size(g, &mate[0]);
    vertex_size_t matching_size_solution = boost::matching_size(g, &solution[0]);

    VertexIterator start, end;
    for (std::tie(start, end) = boost::vertices(g); start != end; start++) {
        Vertex v = *start;
        if (mate[v] != g.null_vertex() && !boost::edge(v, mate[v], g).second) {
            std::cout << v << ", " << mate[v] << std::endl;
            throw "matching incorrect, invalid edge";
        }
    }

    if (matching_size < matching_size_solution) {
		std::cout << std::to_string(matching_size) + " < " + std::to_string(matching_size_solution) << std::endl;
		throw "matching not maximal";
    } else if (matching_size > matching_size_solution) {
        throw "matching incorrect, too large";
    }
}