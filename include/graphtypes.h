//
// Created by suem on 10/27/16.
//

#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjVertexIterator;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef std::vector<Vertex> VertexVector;
