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
float get_error_delta(const int mij, const int psize, const int nsize,
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
        prior_bwd_edge(N), color(N), block_size(N, 1), bitsize(N) {

    for (auto i{0}; i < N; ++i)
      color[i] = i;
  }

  void compress(options& opt);

  float get_cost(options &opt, const int v, const int t);
  void move(const int v, const int t);

  float epsilon{1e-3};

  graph_struct &data;
  graph_struct model;

  float get_objective();

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

  // for the brute-force thing
  std::vector<float> bitsize;

  std::vector<int> util;

  void check(float nbits, float incr_nbits);
  std::vector<std::vector<float>> block_bits;
  std::vector<std::vector<float>> error_bits;
  std::vector<std::vector<float>> bldlt_bits;
  std::vector<std::vector<float>> erdlt_bits;
  std::vector<std::vector<float>> iblck_bits;
  std::vector<std::vector<float>> ierrr_bits;
};

template <class graph_struct>
void block_model<graph_struct>::check(float nbits, float incr_nbits) {

  // std::cout << "CHECK STRUCT\n";

  auto n{data.capacity()};
  auto m{model.size()};

  // check block_size
  util.clear();
  util.resize(n, 0);
  for (auto v : data.nodes) {
    ++util[color[v]];
  }
  for (auto u : model.nodes) {
    assert(block_size[u] == util[u]);
  }

  // check edge count
  auto total{0};
  util.clear();
  util.resize(m * m, 0);
  for (auto u : data.nodes) {
    auto i{model.nodes.index(color[u])};

    for (auto e : data[u]) {
      auto j{model.nodes.index(color[e.endpoint()])};
      ++util[i * m + j];
    }
  }
  for (auto u : model.nodes) {
    auto i{model.nodes.index(color[u])};
    for (auto e : model[u]) {
      auto j{model.nodes.index(color[e.endpoint()])};

      if (e.weight() != util[i * m + j]) {
        std::cout << "problem with " << model.nodes[i] << "," << model.nodes[j]
                  << "(" << e.weight() << "/" << util[i * m + j] << "):\n"
                  << model << std::endl
                  << "   ";
        for (int ii = 0; ii < m; ++ii)
          std::cout << std::setw(3) << model.nodes[ii];
        std::cout << std::endl;
        for (int ii = 0; ii < m; ++ii) {
          std::cout << std::setw(3) << model.nodes[ii];
          for (int jj = 0; jj < m; ++jj) {
            std::cout << std::setw(3) << util[ii * m + jj];
          }
          std::cout << std::endl;
        }
        exit(1);
      }
      assert(e.weight() == util[i * m + j]);
      total += e.weight();
    }
  }

  // std::cout << "CHECK BF/INCR\n";
	
	
	// auto N{model.capacity()};

  // int ok{(std::abs(nbits - incr_nbits) <= 1e-3 ? 0 : N*N)};

	auto ok{std::abs(nbits - incr_nbits) <= 1e-3};

  if (!ok) {
    std::cout << "global problem: " << nbits << " / " << incr_nbits
              << std::endl;
    // exit(1);
  }

	int ierr{-1};
	int jerr{-1};
	int terr{-1};

  float cibits{0}, cnbits{0};
  for (int i = 0; i < data.capacity(); ++i) {
    for (int j = 0; j < data.capacity(); ++j) {
      iblck_bits[i][j] += bldlt_bits[i][j];
      ierrr_bits[i][j] += erdlt_bits[i][j];
      cibits += iblck_bits[i][j];
      cibits += ierrr_bits[i][j];
      cnbits += block_bits[i][j];
      cnbits += error_bits[i][j];
      if (std::abs(block_bits[i][j] - iblck_bits[i][j]) > 1e-3) {
				ierr = i; jerr = j; terr = 0;
        ok = false;
      }
      if (std::abs(error_bits[i][j] - ierrr_bits[i][j]) > 1e-3) {
				ierr = i; jerr = j; terr = 1;
        ok = false;
      }
    }
  }

  if (std::abs(nbits - cnbits) > 1e-3) {
    std::cout << "BF count: " << cnbits << " / " << nbits << std::endl;
  }

  if (std::abs(incr_nbits - cibits) > 1e-3) {
    std::cout << "INCR count: " << cibits << " / " << incr_nbits << std::endl;
  }

  if (!ok) {
    std::cout << "\nDELTA MODEL SIZE:\n";
    for (int i = 0; i < data.capacity(); ++i) {
      for (int j = 0; j < data.capacity(); ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << bldlt_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nINCR MODEL SIZE:\n";
    for (int i = 0; i < data.capacity(); ++i) {
      for (int j = 0; j < data.capacity(); ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << iblck_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nREAL MODEL SIZE:\n";
    for (int i = 0; i < data.capacity(); ++i) {
      for (int j = 0; j < data.capacity(); ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << block_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }

    std::cout << "\nDELTA ERROR SIZE:\n";
    for (int i = 0; i < data.capacity(); ++i) {
      for (int j = 0; j < data.capacity(); ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << erdlt_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nINCR ERROR SIZE:\n";
    for (int i = 0; i < data.capacity(); ++i) {
      for (int j = 0; j < data.capacity(); ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << ierrr_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nREAL ERROR SIZE:\n";
    for (int i = 0; i < data.capacity(); ++i) {
      for (int j = 0; j < data.capacity(); ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << error_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }

	// if (!ok) {
		
		std::cout << (terr ? "error" : "block") << " bits discrepancy on " << ierr << "," << jerr << std::endl;
		
    exit(1);
  }

  assert(total == data.num_edges);
}

template <class graph_struct> float block_model<graph_struct>::get_objective() {
  float nbits{0};

  std::cout << "BRUTE FORCE\n";
  // std::cout << "    ";
  // for (auto u : model.nodes) {
  //   std::cout << std::setw(5) << u << " ";
  //   // std::fill(begin(_bits))
  // }

  for (int i = 0; i < model.capacity(); ++i) {
    std::fill(begin(block_bits[i]), end(block_bits[i]), 0);
    std::fill(begin(error_bits[i]), end(error_bits[i]), 0);
  }

  for (auto u : model.nodes) {
    // std::fill(begin(bitsize), end(bitsize), 0);
    // std::cout << std::setw(3) << u << " ";

    // model size
    for (auto v : model.nodes) {

      assert((block_size[u] * block_size[v]) > 0);

      auto block_sz{
          std::log2(static_cast<float>(block_size[u] * block_size[v] + 1))};

      nbits += block_sz;

      block_bits[u][v] = block_sz;
    }

    // error size
    for (auto e : model[u]) {
      auto v{e.vertex};
      auto nedge{e.weight()};
      auto size{block_size[u] * block_size[v]};

      if (nedge == 0 or nedge == size)
        continue;

      float fnedge{static_cast<float>(nedge)};
      float fsize{static_cast<float>(size)};
      float density{fnedge / fsize};

      auto error_sz{size * ((density - 1) * std::log2(1 - density) -
                            density * std::log2(density))};
      nbits += error_sz;

      error_bits[u][v] += error_sz;
    }

    // for (auto v : model.nodes) {
    //   std::cout << std::setw(5) << std::setprecision(2)
    //             << (error_bits[u][v] - block_bits[u][v] - bldlt_bits[u][v])
    //             << " ";
    // }
    // std::cout << std::endl;
  }

  return nbits;
}

template <class graph_struct>
float block_model<graph_struct>::get_cost(options &opt, const int v,
                                          const int t) {

  // std::cout << "COMPUTE DELTA\n";

  if (opt.checked) {
    for (int i = 0; i < data.capacity(); ++i) {
      std::fill(begin(bldlt_bits[i]), end(bldlt_bits[i]), 0);
      std::fill(begin(erdlt_bits[i]), end(erdlt_bits[i]), 0);
    }
  }

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

	assert(target_bwd_edgecount[o] == origin_fwd_edgecount[t]);
	assert(target_fwd_edgecount[o] == origin_bwd_edgecount[t]);


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

    auto ojdelta = 0.0;
    if (nsizeo > 0)
      ojdelta += std::log2(static_cast<float>(nsizeo + 1));
    ojdelta -= std::log2(static_cast<float>(psizeo + 1));
    if (j != o and j != t)
      ojdelta *= 2;

    auto tjdelta = 0.0;
    if (nsizet > 0)
      tjdelta += std::log2(static_cast<float>(nsizet + 1));
    tjdelta -= std::log2(static_cast<float>(psizet + 1));
    if (j != o and j != t)
      tjdelta *= 2;

    ijdelta = (ojdelta + tjdelta);

    float err;

    if (opt.checked) {
      err = std::log2(static_cast<float>(nsizeo + 1)) -
            std::log2(static_cast<float>(psizeo + 1));
      bldlt_bits[j][o] += err;
      if (j != o and j != t)
        bldlt_bits[o][j] += err;

      err = std::log2(static_cast<float>(nsizet + 1)) -
            std::log2(static_cast<float>(psizet + 1));
      bldlt_bits[j][t] += err;
      if (j != o and j != t)
        bldlt_bits[t][j] += err;
    }

    // if (j == o) {
    //   std::cout << "delta of removing " << v << " from bag " << o
    //             << " that has " << fwd_neighbors[j] << " fwd neighbors, and "
    //             << bwd_neighbors[j] << " bwd neighbors out of "
    //             << origin_bwd_edgecount[j] << ". old size of the bag was "
    //             << psizeo << ", new size is " << nsizeo << std::endl;
    //
    //   std::cout << "delta of adding " << v << " to bag " << t
    //             << " that has " << fwd_neighbors[j] << " fwd neighbors, and "
    //             << bwd_neighbors[j] << " bwd neighbors out of "
    //             << origin_bwd_edgecount[j] << " in " << o << ". old size of the bag was "
    //             << psizet << ", new size is " << nsizet << std::endl;
    // }
		

    if (j == t) {

      assert(target_bwd_edgecount[j] == target_fwd_edgecount[j]);
			

      // error on t -> t
      err =
          get_error_delta(target_bwd_edgecount[t] // + target_fwd_edgecount[j]
                          ,
                          psizet, nsizet, bwd_neighbors[t] + fwd_neighbors[t]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[t][t] += err;
			
			// error on t -> o
			
			err = get_error_delta(origin_bwd_edgecount[t], psizeo, nsizeo, fwd_neighbors[o] - bwd_neighbors[t]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[t][o] += err;
			
			
		} else if (j == o) {

      assert(origin_bwd_edgecount[o] == origin_fwd_edgecount[o]);

      // error_delta on o -> o
      err =
          get_error_delta(origin_bwd_edgecount[o] //+ origin_fwd_edgecount[j]
                          ,
                          psizeo, nsizeo, -bwd_neighbors[o] - fwd_neighbors[o]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[o][o] += err;

      // if (j == o) {
      // 	std::cout << "delta[" << o << "][" << j << "]" << err << std::endl;
      // }
			
			
			// error on o -> t
			err = get_error_delta(target_bwd_edgecount[o], psizet, nsizet, bwd_neighbors[o] - fwd_neighbors[t]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[o][t] += err;
			

    } else {

      // error on j -> t
      err = get_error_delta(target_bwd_edgecount[j], psizet, nsizet,
                            bwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[j][t] += err;
			
      // if (j == o) {
      // 	std::cout << "(t) delta[" << o << "][" << t << "]" << err << std::endl;
      // }

      // error_delta on t -> j
      err = get_error_delta(target_fwd_edgecount[j], psizet, nsizet,
                            fwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[t][j] += err;
			
      // if (j == o) {
      // 	std::cout << "(t) delta[" << t << "][" << o << "]" << err << std::endl;
      // }


      // error_delta on j -> o
      err = get_error_delta(origin_bwd_edgecount[j], psizeo, nsizeo,
                            -bwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[j][o] += err;

      // if (j == t) {
      // 	std::cout << "(o) delta[" << t << "][" << o << "]" << err << std::endl;
      // }

      // error_delta on o -> j
      err = get_error_delta(origin_fwd_edgecount[j], psizeo, nsizeo,
                            -fwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked)
        erdlt_bits[o][j] += err;

      // if (j == t) {
      // 	std::cout << "(o) delta[" << o << "][" << t << "]" << err << std::endl;
      // }
    }

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

    // if (j == 5 or j == 3) {
    //   std::cout << j << " -- "
    //             << "bwd: " << target_bwd_edgecount[j] << "+" <<
    //             bwd_neighbors[j]
    //             << ", fwd:" << target_fwd_edgecount[j] << "+"
    //             << fwd_neighbors[j] << std::endl;
    //   std::cout << j << " -- "
    //             << "bwd: " << origin_bwd_edgecount[j] << "+" <<
    //             bwd_neighbors[j]
    //             << ", fwd:" << origin_fwd_edgecount[j] << "+"
    //             << fwd_neighbors[j] << std::endl;
    // }

    target_bwd_edgecount[j] += bwd_neighbors[j];
    target_fwd_edgecount[j] += fwd_neighbors[j];
    if (j == t) {
      // when adding a node to a bag t, add both its forward and backward edges
      // in the mix
      target_bwd_edgecount[j] += fwd_neighbors[j];
      target_fwd_edgecount[j] += bwd_neighbors[j];
      origin_bwd_edgecount[j] += fwd_neighbors[o];
      origin_fwd_edgecount[j] += bwd_neighbors[o];
    }

    origin_bwd_edgecount[j] -= bwd_neighbors[j];
    origin_fwd_edgecount[j] -= fwd_neighbors[j];
    if (j == o) {
      // when removing a node from a bag o, remove both its forward and backward
      // edges from the mix
      origin_bwd_edgecount[j] -= fwd_neighbors[j];
      origin_fwd_edgecount[j] -= bwd_neighbors[j];
      target_bwd_edgecount[o] -= fwd_neighbors[t];
      target_fwd_edgecount[o] -= bwd_neighbors[t];
    }

    // std::cout << "bwd: " << target_bwd_edgecount[3] << ", fwd:" <<
    // target_fwd_edgecount[3] << std::endl;
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

  for (auto e : model.successors[t])
    if (target_fwd_edgecount[e.endpoint()] != e.weight())
      model.set_weight(t, e.endpoint(), e, target_fwd_edgecount[e.endpoint()]);

  for (auto e : model.predecessors[t])
    if (target_bwd_edgecount[e.endpoint()] != e.weight())
      model.set_weight(e.endpoint(), t, e, target_bwd_edgecount[e.endpoint()]);

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
}

template <class graph_struct> void block_model<graph_struct>::compress(options& opt) {
  if (opt.checked) {
    block_bits.resize(data.capacity());
    error_bits.resize(data.capacity());
    iblck_bits.resize(data.capacity());
    ierrr_bits.resize(data.capacity());
    bldlt_bits.resize(data.capacity());
    erdlt_bits.resize(data.capacity());
    for (auto i = 0; i < data.capacity(); ++i) {
      block_bits[i].resize(data.capacity(), 1);
      error_bits[i].resize(data.capacity(), 0);
      bldlt_bits[i].resize(data.capacity(), 0);
      erdlt_bits[i].resize(data.capacity(), 0);
      iblck_bits[i].resize(data.capacity(), 1);
      ierrr_bits[i].resize(data.capacity(), 0);
    }
  }

  intstack movable(data.capacity());
  movable.fill();

  float nbits{0};
	if (opt.checked) {
			nbits = get_objective();
			std::cout << nbits << " / " << (model.size() * model.size()) << std::endl;
		}
  auto incr_nbits{nbits};

  //
  // exit(1);

  auto prev{-1};

  while (true) {
    bool moved;

    // try to find a move of negative cost
    for (auto v : movable) {

      moved = false;
      for (auto t : model.nodes) {
        auto o{color[v]};
        if (o != t) {

          float delta{get_cost(opt, v, t)};

          // std::cout << delta << std::endl;

          if (delta < -epsilon) {

            incr_nbits += delta;
            // if (opt.checked) {
            //   check();
            //   nbits = get_objective();
            // }

            // std::cout << model << std::endl;

            move(v, t);
            movable.remove(t);
            moved = true;

            float pnbits{nbits};
            if (opt.checked) {
              nbits = get_objective();
              // check();
            }

            if (opt.verbosity > 0)
              std::cout << std::setw(4) << model.size() << " " << std::setw(4)
                        << model.num_edges << " "
                        << "move " << v << ": " << o << " -> " << t << " ("
                        << delta;

            if (opt.checked)
              std::cout << "/" << (nbits - pnbits);

            std::cout << ")" << std::endl;

            // std::cout << nbits << " / " << incr_nbits
            //           << std::endl;

            if (opt.checked) {
              // nbits = get_objective();
              check(nbits, incr_nbits);
            }

            // if (prev == o * t)
            //   exit(1);
            // else
            //   prev = o * t;

            break;
          }
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
