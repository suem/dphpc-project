//
// Created by suem on 10/27/16.
//

#pragma once

#include "graphtypes.h"

void pothen_fan(const Graph& g, MateMap& mate) {

    for (auto& m : mate) {
        m = g.null_vertex();
    }


}
