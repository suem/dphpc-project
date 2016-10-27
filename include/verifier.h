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

inline void verify_matching(const Graph& g, const MateMap& mate) {
    vertex_size_t n = num_vertices(g);
    MateMap solution(n);
    bool success = boost::checked_edmonds_maximum_cardinality_matching(g, &solution[0]);
    if (!success) throw "Graph not bipartite";

    vertex_size_t matching_size = boost::matching_size(g, &mate[0]);
    vertex_size_t matching_size_solution = boost::matching_size(g, &solution[0]);

    if (matching_size < matching_size_solution) {
        throw "matching not maximal";
    } else if (matching_size > matching_size_solution) {
        throw "matching incorrect, too large";
    }
}