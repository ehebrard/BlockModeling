#ifndef __BLOCK__ALGORITHM_HH
#define __BLOCK__ALGORITHM_HH

#include <vector>

#include "intstack.hpp"

namespace block
{

template <class graph_struct>
void BFS(graph_struct &g, intstack &order, int start = 0) {
  for (auto i{start}; i < order.size(); ++i)
    for (auto e : g[order[i]])
      order.add(e.endpoint());
}

// compute the error delta of a move w.r.t. blocks i -> j
// there are 'mij' edges from i to j and 'deltaij' are added
float get_error(const int mij, const int psize, const int nsize,
                const int deltaij) {

  float delta{0};

  float dij = static_cast<float>(mij) / static_cast<float>(psize);
  if (mij > 0 and mij < psize)
    delta -= static_cast<float>(psize) *
             ((dij - 1.0) * std::log2(1.0 - dij) - dij * std::log2(dij));

  if (nsize > 0) {
    float dpij = static_cast<float>(mij + deltaij) / static_cast<float>(nsize);
    if (mij + deltaij > 0 and mij + deltaij < nsize) {
      delta += static_cast<float>(nsize) *
               ((dpij - 1.0) * std::log2(1.0 - dpij) - dpij * std::log2(dpij));
    }
  }

  return delta;
}

// float get_prev_error(const int mij, const int psize) {
//
//   auto delta = 0;
//
//   float dij = static_cast<float>(mij) / static_cast<float>(psize);
//   if (mij > 0 and mij < psize)
//     delta -=
//         static_cast<float>(psize) * ((dij - 1.0) * std::log2(1.0 - dij) - dij
//         * std::log2(dij));
//
//   if (nsize > 0) {
//     float dpij = static_cast<float>(mij + deltaij) /
//     static_cast<float>(nsize);
//     if (mij + deltaij > 0 and mij + deltaij < nsize)
//       delta += static_cast<float>(nsize) *
//                ((dpij - 1.0) * std::log2(1.0 - dpij) - dpij *
//                std::log2(dpij));
//   }
//
//   return delta;
// }

template <class graph_struct> class block_model {

public:
  block_model(graph_struct &g)
      : data(g), model(g), N{g.capacity()}, fwd_neighbors(N), bwd_neighbors(N),
        origin_fwd_edgecount(N), target_fwd_edgecount(N),
        origin_bwd_edgecount(N), target_bwd_edgecount(N), prior_fwd_edge(N),
        prior_bwd_edge(N), color(N), block_size(N, 1) {

    for (auto i{0}; i < N; ++i)
      color[i] = i;
  }

  void compress(options& opt);

  float get_cost(const int v, const int t);
  void move(const int v, const int t);

  graph_struct &data;
  graph_struct model;

private:
  size_t N;

  std::vector<int> fwd_neighbors;
  std::vector<int> bwd_neighbors;

  std::vector<wtype> origin_fwd_edgecount;
  std::vector<wtype> target_fwd_edgecount;
  std::vector<wtype> origin_bwd_edgecount;
  std::vector<wtype> target_bwd_edgecount;

  std::vector<bool> prior_fwd_edge;
  std::vector<bool> prior_bwd_edge;

  std::vector<int> color;
  std::vector<int> block_size;

  std::vector<int> new_fwd_edges;
  std::vector<int> new_bwd_edges;
};

template <class graph_struct>
float block_model<graph_struct>::get_cost(const int v, const int t) {

  float delta{0}, ijdelta{0};

  int o{color[v]};

  std::fill(begin(prior_fwd_edge), end(prior_fwd_edge), false);
  std::fill(begin(prior_bwd_edge), end(prior_bwd_edge), false);

  std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
  std::fill(begin(bwd_neighbors), end(bwd_neighbors), 0);

  std::fill(begin(origin_fwd_edgecount), end(origin_fwd_edgecount), 0);
  std::fill(begin(target_fwd_edgecount), end(target_fwd_edgecount), 0);
  std::fill(begin(origin_bwd_edgecount), end(origin_bwd_edgecount), 0);
  std::fill(begin(target_bwd_edgecount), end(target_bwd_edgecount), 0);

  for (auto e : model.successors[o])
    origin_fwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.successors[t]) {
    target_fwd_edgecount[e.endpoint()] = e.weight();
    prior_fwd_edge[e.endpoint()] = true;
  }
  for (auto e : model.predecessors[o])
    origin_bwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.predecessors[t]) {
    target_bwd_edgecount[e.endpoint()] = e.weight();
    prior_bwd_edge[e.endpoint()] = true;
  }

  for (auto e : data.successors[v]) {
    ++fwd_neighbors[color[e.endpoint()]];
  }

  for (auto e : data.predecessors[v]) {
    ++bwd_neighbors[color[e.endpoint()]];
  }

  // std::cout << "     ";
  // for (auto j : model.nodes) {
  //   std::cout << std::setw(5) << j << " ";
  // }
  // std::cout << std::endl << "d_oj ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << origin_fwd_edgecount[j] << " ";
  // std::cout << std::endl << "d_tj ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << target_fwd_edgecount[j] << " ";
  // std::cout << std::endl << "d_jo ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << origin_bwd_edgecount[j] << " ";
  // std::cout << std::endl << "d_jt ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << target_bwd_edgecount[j] << " ";
  // std::cout << std::endl << "N+_j ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << fwd_neighbors[j] << " ";
  // std::cout << std::endl << "N-_j ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << bwd_neighbors[j] << " ";
  // std::cout << std::endl << "pszo ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << block_size[o] * block_size[j] << " ";
  // std::cout << std::endl << "nszo ";
  // for (auto j : model.nodes) {
  //   std::cout << std::setw(5);
  //   if (j != o and j != t) {
  //     std::cout << (block_size[o] - 1) * block_size[j];
  //   } else if (j == o) {
  //     std::cout << (block_size[o] - 1) * (block_size[o] - 1);
  //   } else {
  //     std::cout << (block_size[o] - 1) * (block_size[t] + 1);
  //   }
  //   std::cout << " ";
  // }
  // std::cout << std::endl << "pszt ";
  // for (auto j : model.nodes)
  //   std::cout << std::setw(5) << block_size[t] * block_size[j] << " ";
  // std::cout << std::endl << "nszt ";
  // for (auto j : model.nodes) {
  //   std::cout << std::setw(5);
  //   if (j != o and j != t) {
  //     std::cout << (block_size[t] + 1) * block_size[j];
  //   } else if (j == o) {
  //     std::cout << (block_size[t] + 1) * (block_size[o] - 1);
  //   } else {
  //     std::cout << (block_size[t] + 1) * (block_size[t] + 1);
  //   }
  //   std::cout << " ";
  // }
  // std::cout << std::endl << "delt ";
  // for (auto j : model.nodes) {
  //   auto psizeo{block_size[o] * block_size[j]};
  //   auto psizet{block_size[t] * block_size[j]};
  //   auto nsizeo{0};
  //   auto nsizet{0};
  //
  //   if (j != o and j != t) {
  //     nsizeo = (block_size[o] - 1) * block_size[j];
  //     nsizet = (block_size[t] + 1) * block_size[j];
  //   } else if (j == o) {
  //     nsizeo = (block_size[o] - 1) * (block_size[o] - 1);
  //     nsizet = (block_size[o] - 1) * (block_size[t] + 1);
  //   } else {
  //     nsizet = (block_size[t] + 1) * (block_size[t] + 1);
  //     nsizeo = (block_size[o] - 1) * (block_size[t] + 1);
  //   }
  //
  //   ijdelta = 0.0;
  //   if (nsizeo > 0)
  //     ijdelta += std::log2(static_cast<float>(nsizeo + 1));
  //   ijdelta -= std::log2(static_cast<float>(psizeo + 1));
  //
  //   if (nsizet > 0)
  //     ijdelta += std::log2(static_cast<float>(nsizet + 1));
  //   ijdelta -= std::log2(static_cast<float>(psizet + 1));
  //
  //   std::cout << std::setw(5) << std::setprecision(2) << ijdelta << " ";
  // }
  // std::cout << std::endl << "     ";

  // update the densities
  for (auto j : model.nodes) {

    auto psizeo{block_size[o] * block_size[j]};
    auto psizet{block_size[t] * block_size[j]};
    auto nsizeo{0};
    auto nsizet{0};

    if (j != o and j != t) {
      nsizeo = (block_size[o] - 1) * block_size[j];
      nsizet = (block_size[t] + 1) * block_size[j];
    } else if (j == o) {
      nsizeo = (block_size[o] - 1) * (block_size[o] - 1);
      nsizet = (block_size[o] - 1) * (block_size[t] + 1);
    } else {
      nsizet = (block_size[t] + 1) * (block_size[t] + 1);
      nsizeo = (block_size[o] - 1) * (block_size[t] + 1);
    }

    ijdelta = 0.0;
    if (nsizeo > 0)
      ijdelta += std::log2(static_cast<float>(nsizeo + 1));
    ijdelta -= std::log2(static_cast<float>(psizeo + 1));

    if (nsizet > 0)
      ijdelta += std::log2(static_cast<float>(nsizet + 1));
    ijdelta -= std::log2(static_cast<float>(psizet + 1));

    // error on j -> t
    ijdelta +=
        get_error(target_bwd_edgecount[j], psizet, nsizet, bwd_neighbors[j]);

    // error on t -> j
    ijdelta +=
        get_error(target_fwd_edgecount[j], psizet, nsizet, fwd_neighbors[j]);

    // error on j -> o
    ijdelta +=
        get_error(origin_bwd_edgecount[j], psizeo, nsizeo, -bwd_neighbors[j]);

    // error on o -> j
    ijdelta +=
        get_error(origin_fwd_edgecount[j], psizeo, nsizeo, -fwd_neighbors[j]);

    delta += ijdelta;

    // std::cout << std::setw(5) << std::setprecision(2) << ijdelta << " ";
  }
  // std::cout << " => " << delta << std::endl;

  return delta;
}

template <class graph_struct>
void block_model<graph_struct>::move(const int v, const int t) {
  assert(new_bwd_edges.size() == 0);
  assert(new_fwd_edges.size() == 0);

  auto o{color[v]};
  color[v] = t;

  ++block_size[t];
  --block_size[o];

  for (auto j : model.nodes) {
    if (!prior_bwd_edge[j] and bwd_neighbors[j] > 0) {
      new_bwd_edges.push_back(j);
      new_bwd_edges.push_back(t);
      //
      // std::cout << "there were " << target_bwd_edgecount[j] << " edges from "
      //           << j << " to " << t << " but " << bwd_neighbors[j] << " to "
      //           << v << std::endl;
    }
    if (!prior_fwd_edge[j] and fwd_neighbors[j] > 0) {
      new_fwd_edges.push_back(j);
      new_fwd_edges.push_back(t);
      //
      // std::cout << "there were " << target_fwd_edgecount[j] << " edges from "
      //           << t << " to " << j << " but " << fwd_neighbors[j] << " from
      //           "
      //           << v << std::endl;
    }
    // if (block_size[o] > 0) {
    //   if (origin_bwd_edgecount[j] == 0 and bwd_neighbors[j] > 0) {
    //     new_bwd_edges.push_back(j);
    //     new_bwd_edges.push_back(o);
    //
    // 	      std::cout << "there were " << origin_bwd_edgecount[j] << " edges
    // from "
    // 	                << j << " to " << o << " but " << bwd_neighbors[j] << "
    // to "
    // 	                << v << std::endl;
    //   }
    //   if (origin_fwd_edgecount[j] == 0 and fwd_neighbors[j] > 0) {
    //     new_fwd_edges.push_back(j);
    //     new_fwd_edges.push_back(o);
    //   }
    // }

    target_bwd_edgecount[j] += bwd_neighbors[j];
    target_fwd_edgecount[j] += fwd_neighbors[j];
    origin_bwd_edgecount[j] -= bwd_neighbors[j];
    origin_fwd_edgecount[j] -= fwd_neighbors[j];
  }

  while (new_bwd_edges.size() > 0) {
    auto i{new_bwd_edges.back()};
    new_bwd_edges.pop_back();
    auto j{new_bwd_edges.back()};
    new_bwd_edges.pop_back();

    // std::cout << "add " << j << "-(" << target_bwd_edgecount[j] << ")->" << i
    //           << std::endl;

    model.add_edge(j, i, target_bwd_edgecount[j]);
		
  }
  while (new_fwd_edges.size() > 0) {
    auto i{new_fwd_edges.back()};
    new_fwd_edges.pop_back();
    auto j{new_fwd_edges.back()};
    new_fwd_edges.pop_back();

    // std::cout << "add " << i << "-(" << target_fwd_edgecount[j] << ")->" << j
    //           << std::endl;

    model.add_edge(i, j, target_fwd_edgecount[j]);
		
  }

  if (block_size[o] == 0) {

    // std::cout << "  - remove " << o << std::endl;

    model.rem_node(o);
  } else {
    for (auto e : model.successors[o])
      if (origin_fwd_edgecount[e.endpoint()] != e.weight())
        model.set_weight(o, e.endpoint(), e, origin_fwd_edgecount[e.endpoint()]);
    for (auto e : model.predecessors[o])
      if (origin_bwd_edgecount[e.endpoint()] != e.weight())
        model.set_weight(e.endpoint(), o, e, origin_bwd_edgecount[e.endpoint()]);
  }

  for (auto e : model.successors[t])
    if (target_fwd_edgecount[e.endpoint()] != e.weight())
      model.set_weight(t, e.endpoint(), e, target_fwd_edgecount[e.endpoint()]);

  for (auto e : model.predecessors[t])
    if (target_bwd_edgecount[e.endpoint()] != e.weight())
      model.set_weight(e.endpoint(), t, e, target_bwd_edgecount[e.endpoint()]);
}

template <class graph_struct> void block_model<graph_struct>::compress(options& opt) {

  intstack movable(data.capacity());
  movable.fill();

  while (true) {
    bool moved;

    // try to find a move of negative cost
    for (auto v : movable) {

      moved = false;
      for (auto t : model.nodes)
        if (color[v] != t) {

          float delta{get_cost(v, t)};

          if (delta < 0) {
						
						if(opt.verbosity > 0)
		            std::cout << "move " << v << ": " << color[v] << " -> " << t << " ("
		                      << delta << ")" << std::endl;

            move(v, t);
            movable.remove(t);
            moved = true;
            break;
          }
        }

      if (moved)
        break;
    }

    if (!moved)
      break;
  }
}


} // namespace block

#endif
