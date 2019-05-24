#ifndef __BLOCK__GLOBAL_H
#define __BLOCK__GLOBAL_H

#include <random>
#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "graph.hpp"
#include "dyngraph.hpp"

static bool print_debug{false};

#define DBG(code)                                                              \
  if (print_debug) {                                                           \
    code                                                                       \
  }

typedef float wtype;
auto max_weight{std::numeric_limits<wtype>::max()};
auto min_weight{std::numeric_limits<wtype>::min()};

typedef block::graph<wtype> sgraph;
typedef block::weighted_edge<wtype> edge;

typedef block::dyngraph<wtype> dgraph;
typedef block::neighbor_info<wtype> dedge;

#endif //_GLOBAL_H

