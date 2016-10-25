#include <iostream>
#include "../include/graphloader.h"

using namespace std;
using namespace boost;

void loadGraph(Graph& graph) {
    unsigned int u,v;
    while(cin >> u && cin >> v) {
        add_edge(u, v, graph);
    }
}

