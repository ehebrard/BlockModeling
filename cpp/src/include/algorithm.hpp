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
      : N{g.capacity()}, data(g), model(g), fwd_neighbors(N), bwd_neighbors(N),
        origin_fwd_edgecount(N), target_fwd_edgecount(N),
        origin_bwd_edgecount(N), target_bwd_edgecount(N), color(N),
        block_size(N, 1) {

    for (auto i{0}; i < N; ++i)
      color[i] = i;
  }

  void compress();

  float get_cost(const int v, const int t);
  void move(const int v, const int t);

private:
  size_t N;

  graph_struct &data;
  graph_struct model;

  std::vector<int> fwd_neighbors;
  std::vector<int> bwd_neighbors;

  std::vector<wtype> origin_fwd_edgecount;
  std::vector<wtype> target_fwd_edgecount;
  std::vector<wtype> origin_bwd_edgecount;
  std::vector<wtype> target_bwd_edgecount;

  std::vector<int> color;
  std::vector<int> block_size;

  std::vector<int> new_fwd_edges;
  std::vector<int> new_bwd_edges;
};

template <class graph_struct>
float block_model<graph_struct>::get_cost(const int v, const int t) {

  float delta{0}, ijdelta{0};

  int o{color[v]};

  std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
  std::fill(begin(bwd_neighbors), end(bwd_neighbors), 0);

  std::fill(begin(origin_fwd_edgecount), end(origin_fwd_edgecount), 0);
  std::fill(begin(target_fwd_edgecount), end(target_fwd_edgecount), 0);
  std::fill(begin(origin_bwd_edgecount), end(origin_bwd_edgecount), 0);
  std::fill(begin(target_bwd_edgecount), end(target_bwd_edgecount), 0);

  for (auto e : model.successors[o])
    origin_fwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.successors[t])
    target_fwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.predecessors[o])
    origin_bwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.predecessors[t])
    target_bwd_edgecount[e.endpoint()] = e.weight();

  for (auto e : data.successors[v]) {
    ++fwd_neighbors[color[e.endpoint()]];
  }

  for (auto e : data.predecessors[v]) {
    ++bwd_neighbors[color[e.endpoint()]];
  }

  std::cout << "     ";
  for (auto j : model.nodes) {
    std::cout << std::setw(5) << j << " ";
  }
  std::cout << std::endl << "d_oj ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << origin_fwd_edgecount[j] << " ";
  std::cout << std::endl << "d_tj ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << target_fwd_edgecount[j] << " ";
  std::cout << std::endl << "d_jo ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << origin_bwd_edgecount[j] << " ";
  std::cout << std::endl << "d_jt ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << target_bwd_edgecount[j] << " ";
  std::cout << std::endl << "N+_j ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << fwd_neighbors[j] << " ";
  std::cout << std::endl << "N-_j ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << bwd_neighbors[j] << " ";
  std::cout << std::endl << "pszo ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << block_size[o] * block_size[j] << " ";
  std::cout << std::endl << "nszo ";
  for (auto j : model.nodes) {
    std::cout << std::setw(5);
    if (j != o and j != t) {
      std::cout << (block_size[o] - 1) * block_size[j];
    } else if (j == o) {
      std::cout << (block_size[o] - 1) * (block_size[o] - 1);
    } else {
      std::cout << (block_size[o] - 1) * (block_size[t] + 1);
    }
    std::cout << " ";
  }
  std::cout << std::endl << "pszt ";
  for (auto j : model.nodes)
    std::cout << std::setw(5) << block_size[t] * block_size[j] << " ";
  std::cout << std::endl << "nszt ";
  for (auto j : model.nodes) {
    std::cout << std::setw(5);
    if (j != o and j != t) {
      std::cout << (block_size[t] + 1) * block_size[j];
    } else if (j == o) {
      std::cout << (block_size[t] + 1) * (block_size[o] - 1);
    } else {
      std::cout << (block_size[t] + 1) * (block_size[t] + 1);
    }
    std::cout << " ";
  }
  std::cout << std::endl << "delt ";
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

    std::cout << std::setw(5) << std::setprecision(2) << ijdelta << " ";
  }
  std::cout << std::endl << "     ";

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

    // // here we need to check for new
    // edges!!!!!
    // if(target_fwd_edgecount[j])
    // 	{
    //
    // 	}
    //
    // target_bwd_edgecount[j] +=
    // bwd_neighbors[j];
    // target_fwd_edgecount[j] +=
    // fwd_neighbors[j];
    // origin_bwd_edgecount[j] -=
    // bwd_neighbors[j];
    // origin_fwd_edgecount[j] -=
    // fwd_neighbors[j];

    delta += ijdelta;

    std::cout << std::setw(5) << std::setprecision(2) << ijdelta << " ";
  }
  std::cout << " => " << delta << std::endl;

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
    if (target_bwd_edgecount[j] == 0 and bwd_neighbors[j] > 0) {
      new_bwd_edges.push_back(j);
      new_bwd_edges.push_back(t);
    }
    if (target_fwd_edgecount[j] == 0 and fwd_neighbors[j] > 0) {
      new_fwd_edges.push_back(j);
      new_fwd_edges.push_back(t);
    }
    if (block_size[o] > 0) {
      if (origin_bwd_edgecount[j] == 0 and bwd_neighbors[j] > 0) {
        new_bwd_edges.push_back(j);
        new_bwd_edges.push_back(o);
      }
      if (origin_fwd_edgecount[j] == 0 and fwd_neighbors[j] > 0) {
        new_fwd_edges.push_back(j);
        new_fwd_edges.push_back(o);
      }
    }

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

		std::cout << "add " << j << "-(" << target_bwd_edgecount[j] << ")->" << i << std::endl;

    model.add_edge(j, i, target_bwd_edgecount[j]);
		
  }
  while (new_fwd_edges.size() > 0) {
    auto i{new_fwd_edges.back()};
    new_fwd_edges.pop_back();
    auto j{new_fwd_edges.back()};
    new_fwd_edges.pop_back();

		std::cout << "add " << i << "-(" << target_fwd_edgecount[j] << ")->" << j << std::endl;

    model.add_edge(i, j, target_fwd_edgecount[j]);
		
  }

  if (block_size[o] == 0) {
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

template <class graph_struct> void block_model<graph_struct>::compress() {

  while (true) {
    bool moved;

    // try to find a move of negative cost
    for (auto v : data.nodes) {

      moved = false;
      for (auto t : model.nodes)
        if (color[v] != t) {

          float delta{get_cost(v, t)};

          if (delta < 0) {
            std::cout << "move " << v << " -> " << t << std::endl;

            move(v, t);
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

// template <class graph_struct> void compress(graph_struct &forward) {
//
//   // reversed graph
//   dgraph backward(forward, true);
//
//   std::cout << backward << " " << backward.directed << std::endl;
//
//   // block model
//   dgraph fwd_block_model(forward, false);
//
//   dgraph bwd_block_model(backward, false);
//
//   //   std::cout << fwd_block_model << std::endl;
//   // std::cout << bwd_block_model << std::endl;
//
//   auto N{forward.size()};
//
//   float objective{static_cast<float>(N * N)};
//
//   // std::vector<int> origin_fwd_neighbors(N);
//   // std::vector<int> origin_bwd_neighbors(N);
//   // std::vector<int> target_fwd_neighbors(N);
//   // std::vector<int> target_bwd_neighbors(N);
//   std::vector<int> fwd_neighbors(N);
//   std::vector<int> bwd_neighbors(N);
//
//   std::vector<wtype> origin_fwd_edgecount(N);
//   std::vector<wtype> target_fwd_edgecount(N);
//   std::vector<wtype> origin_bwd_edgecount(N);
//   std::vector<wtype> target_bwd_edgecount(N);
//
//   std::vector<int> color(N);
//   for (auto i{0}; i < N; ++i)
//     color[i] = i;
//
//   std::vector<int> block_size(N, 1);
//
//   std::vector<int> new_edges;
//
//   while (true) {
//     // try to find a move of negative cost
//     for (auto v : forward.nodes) {
//       auto o{color[v]};
//
//       for (auto t : fwd_block_model.nodes)
//         if (o != t) {
//
//           float delta{0}, ijdelta{0};
//
//           std::cout << "check " << v << " -> b" << t << std::endl;
//
//           std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
//           std::fill(begin(bwd_neighbors), end(bwd_neighbors), 0);
//
//           std::fill(begin(origin_fwd_edgecount), end(origin_fwd_edgecount),
//           0);
//           std::fill(begin(target_fwd_edgecount), end(target_fwd_edgecount),
//           0);
//           std::fill(begin(origin_bwd_edgecount), end(origin_bwd_edgecount),
//           0);
//           std::fill(begin(target_bwd_edgecount), end(target_bwd_edgecount),
//           0);
//
//           for (auto e : fwd_block_model[o])
//             origin_fwd_edgecount[e.endpoint()] = e.weight();
//           for (auto e : fwd_block_model[t])
//             target_fwd_edgecount[e.endpoint()] = e.weight();
//           std::cout << o << ":";
//           for (auto e : bwd_block_model[o]) {
//             std::cout << " " << e.endpoint();
//             origin_bwd_edgecount[e.endpoint()] = e.weight();
//           }
//           std::cout << std::endl;
//           for (auto e : bwd_block_model[t])
//             target_bwd_edgecount[e.endpoint()] = e.weight();
//
//           for (auto e : forward[v]) {
//             ++fwd_neighbors[color[e.endpoint()]];
//           }
//
//           for (auto e : backward[v]) {
//             ++bwd_neighbors[color[e.endpoint()]];
//           }
//
//           std::cout << "     ";
//           for (auto j : fwd_block_model.nodes) {
//             std::cout << std::setw(5) << j << " ";
//           }
//           std::cout << std::endl << "d_oj ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << origin_fwd_edgecount[j] << " ";
//           std::cout << std::endl << "d_tj ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << target_fwd_edgecount[j] << " ";
//           std::cout << std::endl << "d_jo ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << origin_bwd_edgecount[j] << " ";
//           std::cout << std::endl << "d_jt ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << target_bwd_edgecount[j] << " ";
//           std::cout << std::endl << "N+_j ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << fwd_neighbors[j] << " ";
//           std::cout << std::endl << "N-_j ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << bwd_neighbors[j] << " ";
//           std::cout << std::endl << "pszo ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << block_size[o] * block_size[j] << "
//             ";
//           std::cout << std::endl << "nszo ";
//           for (auto j : fwd_block_model.nodes) {
//             std::cout << std::setw(5);
//             if (j != o and j != t) {
//               std::cout << (block_size[o] - 1) * block_size[j];
//             } else if (j == o) {
//               std::cout << (block_size[o] - 1) * (block_size[o] - 1);
//             } else {
//               std::cout << (block_size[o] - 1) * (block_size[t] + 1);
//             }
//             std::cout << " ";
//           }
//           std::cout << std::endl << "pszt ";
//           for (auto j : fwd_block_model.nodes)
//             std::cout << std::setw(5) << block_size[t] * block_size[j] << "
//             ";
//           std::cout << std::endl << "nszt ";
//           for (auto j : fwd_block_model.nodes) {
//             std::cout << std::setw(5);
//             if (j != o and j != t) {
//               std::cout << (block_size[t] + 1) * block_size[j];
//             } else if (j == o) {
//               std::cout << (block_size[t] + 1) * (block_size[o] - 1);
//             } else {
//               std::cout << (block_size[t] + 1) * (block_size[t] + 1);
//             }
//             std::cout << " ";
//           }
//           std::cout << std::endl << "delt ";
//           for (auto j : fwd_block_model.nodes) {
//             auto psizeo{block_size[o] * block_size[j]};
//             auto psizet{block_size[t] * block_size[j]};
//             auto nsizeo{0};
//             auto nsizet{0};
//
//             if (j != o and j != t) {
//               nsizeo = (block_size[o] - 1) * block_size[j];
//               nsizet = (block_size[t] + 1) * block_size[j];
//             } else if (j == o) {
//               nsizeo = (block_size[o] - 1) * (block_size[o] - 1);
//               nsizet = (block_size[o] - 1) * (block_size[t] + 1);
//             } else {
//               nsizet = (block_size[t] + 1) * (block_size[t] + 1);
//               nsizeo = (block_size[o] - 1) * (block_size[t] + 1);
//             }
//
//             ijdelta = 0.0;
//             if (nsizeo > 0)
//               ijdelta += std::log2(static_cast<float>(nsizeo + 1));
//             ijdelta -= std::log2(static_cast<float>(psizeo + 1));
//
//             if (nsizet > 0)
//               ijdelta += std::log2(static_cast<float>(nsizet + 1));
//             ijdelta -= std::log2(static_cast<float>(psizet + 1));
//
//             std::cout << std::setw(5) << std::setprecision(2) << ijdelta << "
//             ";
//           }
//           std::cout << std::endl << "     ";
//
//           // update the densities
//           for (auto j : fwd_block_model.nodes) {
//
//             auto psizeo{block_size[o] * block_size[j]};
//             auto psizet{block_size[t] * block_size[j]};
//             auto nsizeo{0};
//             auto nsizet{0};
//
//             if (j != o and j != t) {
//               nsizeo = (block_size[o] - 1) * block_size[j];
//               nsizet = (block_size[t] + 1) * block_size[j];
//             } else if (j == o) {
//               nsizeo = (block_size[o] - 1) * (block_size[o] - 1);
//               nsizet = (block_size[o] - 1) * (block_size[t] + 1);
//             } else {
//               nsizet = (block_size[t] + 1) * (block_size[t] + 1);
//               nsizeo = (block_size[o] - 1) * (block_size[t] + 1);
//             }
//
//             ijdelta = 0.0;
//             if (nsizeo > 0)
//               ijdelta += std::log2(static_cast<float>(nsizeo + 1));
//             ijdelta -= std::log2(static_cast<float>(psizeo + 1));
//
//             if (nsizet > 0)
//               ijdelta += std::log2(static_cast<float>(nsizet + 1));
//             ijdelta -= std::log2(static_cast<float>(psizet + 1));
//
//             // error on j -> t
//             ijdelta += get_error(target_bwd_edgecount[j], psizet, nsizet,
//                                  bwd_neighbors[j]);
//
//             // error on t -> j
//             ijdelta += get_error(target_fwd_edgecount[j], psizet, nsizet,
//                                  fwd_neighbors[j]);
//
//             // error on j -> o
//             ijdelta += get_error(origin_bwd_edgecount[j], psizeo, nsizeo,
//                                  -bwd_neighbors[j]);
//
//             // error on o -> j
//             ijdelta += get_error(origin_fwd_edgecount[j], psizeo, nsizeo,
//                                  -fwd_neighbors[j]);
//
//             // // here we need to check for new
//             // edges!!!!!
//             // if(target_fwd_edgecount[j])
//             // 	{
//             //
//             // 	}
//             //
//             // target_bwd_edgecount[j] +=
//             // bwd_neighbors[j];
//             // target_fwd_edgecount[j] +=
//             // fwd_neighbors[j];
//             // origin_bwd_edgecount[j] -=
//             // bwd_neighbors[j];
//             // origin_fwd_edgecount[j] -=
//             // fwd_neighbors[j];
//
//             delta += ijdelta;
//
//             std::cout << std::setw(5) << std::setprecision(2) << ijdelta << "
//             ";
//           }
//           std::cout << " => " << delta << std::endl;
//
// 					if(delta < 0)
// 						exit(1);
//
//           //           for (auto e :
//           //           fwd_block_model[o])
//           // e.setweight(origin_fwd_edgecount[e.endpoint()]);
//           //           for (auto e :
//           //           fwd_block_model[t])
//           // e.setweight(target_fwd_edgecount[e.endpoint()]);
//           //           for (auto e :
//           //           bwd_block_model[o])
//           // e.setweight(origin_bwd_edgecount[e.endpoint()]);
//           //           for (auto e :
//           //           bwd_block_model[t])
//           // e.setweight(target_bwd_edgecount[e.endpoint()]);
//         }
//     }
//
//     break;
//   }
// }

} // namespace block

#endif
