#ifndef __BLOCK__ALGORITHM_HH
#define __BLOCK__ALGORITHM_HH

#include <vector>

#include "global.hpp"
#include "intstack.hpp"
#include "options.hpp"

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


template <class graph_struct> class block_model {

public:
  block_model(graph_struct &g)
      : data(g), model(g), N{g.capacity()}, fwd_neighbors(N), bwd_neighbors(N),
        origin_fwd_edgecount(N), target_fwd_edgecount(N),
        origin_bwd_edgecount(N), target_bwd_edgecount(N), prior_fwd_edge(N),
        prior_bwd_edge(N), block(N), block_size(N, 1) {

    for (auto u : g.nodes)
      block[u] = u;
  }

  template <class block_struct>
  block_model(graph_struct &g, block_struct &blocks)
      : data(g), N{blocks.size()}, fwd_neighbors(N), bwd_neighbors(N),
        origin_fwd_edgecount(N), target_fwd_edgecount(N),
        origin_bwd_edgecount(N), target_bwd_edgecount(N), prior_fwd_edge(N),
        prior_bwd_edge(N), block(g.capacity()), block_size(N, 0) {

    model.initialise(N);
		

    for (auto b{begin(blocks)}; b!=end(blocks); ++b) {
      // std::cout << (b-begin(blocks)) << ": ";
      block_size[b-begin(blocks)] = b->size();
      for (auto u{b->begin()}; u!=b->end(); ++u) {
        // std::cout << " " << *u;
        block[*u] = (b-begin(blocks));
      }
      // std::cout << std::endl;
    }

		for (auto b{begin(blocks)}; b!=end(blocks); ++b) {
      std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
			for (auto u{b->begin()}; u!=b->end(); ++u) {
        // std::cout << " " << *u << "( ";
        for (auto e : g.successors[*u]) {
          // std::cout << e.endpoint() << ":" << block[e.endpoint()] << " ";
          ++fwd_neighbors[block[e.endpoint()]];
        }
        // std::cout << ")" << std::endl;
      }
      // std::cout << std::endl;

			auto bi{(b-begin(blocks))};
      for (int bj{0}; bj < N; ++bj) {
        // std::cout << " " << fwd_neighbors[bj];
        if (fwd_neighbors[bj] > 0) {
          // std::cout << ":(" << bi << "," << bj << ")";
          model.add_edge(bi, bj, fwd_neighbors[bj]);
        }
      }
      // std::cout << std::endl;
			
			// std::cout << model << std::endl;
    }

    // std::cout << model << std::endl;
  }

  void get_blocks(std::vector<std::vector<int>> &blocks);

  void compress(options &opt, std::mt19937 &rng);

  float get_cost(options &opt, const int v, const int t);
  void move(const int v, const int t);

  float epsilon{1e-3};

  graph_struct &data;
  graph_struct model;

  float get_objective(); // const float A, const float B);

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

  std::vector<int> block;
  std::vector<int> block_size;

  std::vector<int> new_fwd_edges;
  std::vector<int> new_bwd_edges;

  /********* VERIFICATION STUFF ***********/
  // for the brute-force thing
  // std::vector<float> bitsize;
  std::vector<int> util;
  // float err_epsilon{.5};
  void check(float nbits, float incr_nbits, float eps);
  std::vector<std::vector<float>> block_bits;
  std::vector<std::vector<float>> error_bits;
  std::vector<std::vector<float>> bldlt_bits;
  std::vector<std::vector<float>> erdlt_bits;
  std::vector<std::vector<float>> iblck_bits;
  std::vector<std::vector<float>> ierrr_bits;
};

template <class graph_struct>
void block_model<graph_struct>::get_blocks(
    std::vector<std::vector<int>> &blocks) {
  blocks.clear();
  blocks.resize(model.size());
  for (auto x : data.nodes)
    blocks[model.nodes.index(block[x])].push_back(x);
}

template <class graph_struct>
float block_model<graph_struct>::get_cost(options &opt, const int v,
                                          const int t) {

  // float A{opt.alpha};
  // float B{opt.beta};

  if (opt.checked()) {
    for (int i = 0; i < N; ++i) {
      std::fill(begin(bldlt_bits[i]), end(bldlt_bits[i]), 0);
      std::fill(begin(erdlt_bits[i]), end(erdlt_bits[i]), 0);
    }
  }

  float delta{0}, ijdelta{0};

  int o{block[v]};

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
    ++fwd_neighbors[block[e.endpoint()]];
  }

  for (auto e : data.predecessors[v]) {
    ++bwd_neighbors[block[e.endpoint()]];
  }

  assert(target_bwd_edgecount[o] == origin_fwd_edgecount[t]);
  assert(target_fwd_edgecount[o] == origin_bwd_edgecount[t]);

  float err;

  auto ojdelta = 0.0;
  auto tjdelta = 0.0;

  ijdelta = 0.0;

  // o.o
  if (block_size[o] > 1)
    ojdelta += std::log2(
        static_cast<float>((block_size[o] - 1) * (block_size[o] - 1) + 1));
  ojdelta -= std::log2(static_cast<float>(block_size[o] * block_size[o] + 1));
  // ojdelta *= A;
  if (opt.checked())
    bldlt_bits[o][o] += ojdelta;

  // t.t
  // if (nsizet > 0)
  tjdelta += std::log2(
      static_cast<float>((block_size[t] + 1) * (block_size[t] + 1) + 1));
  tjdelta -= std::log2(static_cast<float>(block_size[t] * block_size[t] + 1));
  // tjdelta *= A;
  if (opt.checked())
    bldlt_bits[t][t] += tjdelta;

  // o.t and t.o
  if (block_size[o] > 1)
    ijdelta += std::log2(
        static_cast<float>((block_size[o] - 1) * (block_size[t] + 1) + 1));
  ijdelta -= std::log2(static_cast<float>(block_size[o] * block_size[t] + 1));
  // ijdelta *= A;
  if (opt.checked()) {
    bldlt_bits[t][o] += ijdelta;
    bldlt_bits[o][t] += ijdelta;
  }
  ijdelta *= 2;

  ijdelta += (ojdelta + tjdelta);

  // error on t -> t
  err = // B *
      get_error_delta(target_bwd_edgecount[t], block_size[t] * block_size[t],
                      (block_size[t] + 1) * (block_size[t] + 1),
                      bwd_neighbors[t] + fwd_neighbors[t]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[t][t] += err;

  // error on t -> o
  err = // B *
      get_error_delta(origin_bwd_edgecount[t], block_size[t] * block_size[o],
                      (block_size[t] + 1) * (block_size[o] - 1),
                      fwd_neighbors[o] - bwd_neighbors[t]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[t][o] += err;

  assert(origin_bwd_edgecount[o] == origin_fwd_edgecount[o]);

  // error_delta on o -> o
  err = // B *
      get_error_delta(origin_bwd_edgecount[o], block_size[o] * block_size[o],
                      (block_size[o] - 1) * (block_size[o] - 1),
                      -bwd_neighbors[o] - fwd_neighbors[o]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[o][o] += err;

  // error on o -> t
  err = // B *
      get_error_delta(target_bwd_edgecount[o], block_size[t] * block_size[o],
                      (block_size[t] + 1) * (block_size[o] - 1),
                      bwd_neighbors[o] - fwd_neighbors[t]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[o][t] += err;

  delta += ijdelta;

  // update the densities
  for (auto j : model.nodes)
    if (j != o and j != t) {

      auto psizeo{block_size[o] * block_size[j]};
      auto psizet{block_size[t] * block_size[j]};
      auto nsizeo{(block_size[o] - 1) * block_size[j]};
      auto nsizet{(block_size[t] + 1) * block_size[j]};

      ojdelta = 0.0;
      if (nsizeo > 0)
        ojdelta += std::log2(static_cast<float>(nsizeo + 1));
      ojdelta -= std::log2(static_cast<float>(psizeo + 1));
      ojdelta *= 2;

      tjdelta = 0.0;
      if (nsizet > 0)
        tjdelta += std::log2(static_cast<float>(nsizet + 1));
      tjdelta -= std::log2(static_cast<float>(psizet + 1));
      tjdelta *= 2;

      ijdelta = (ojdelta + tjdelta);

      if (opt.checked()) {
        err = std::log2(static_cast<float>(nsizeo + 1)) -
              std::log2(static_cast<float>(psizeo + 1));
        bldlt_bits[j][o] += err;
        bldlt_bits[o][j] += err;

        err = std::log2(static_cast<float>(nsizet + 1)) -
              std::log2(static_cast<float>(psizet + 1));
        bldlt_bits[j][t] += err;
        bldlt_bits[t][j] += err;
      }

      // error on j -> t
      err = get_error_delta(target_bwd_edgecount[j], psizet, nsizet,
                            bwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked())
        erdlt_bits[j][t] += err;

      // error_delta on t -> j
      err = get_error_delta(target_fwd_edgecount[j], psizet, nsizet,
                            fwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked())
        erdlt_bits[t][j] += err;

      // error_delta on j -> o
      err = get_error_delta(origin_bwd_edgecount[j], psizeo, nsizeo,
                            -bwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked())
        erdlt_bits[j][o] += err;

      // error_delta on o -> j
      err = get_error_delta(origin_fwd_edgecount[j], psizeo, nsizeo,
                            -fwd_neighbors[j]);
      ijdelta += err;
      if (opt.checked())
        erdlt_bits[o][j] += err;

      delta += ijdelta;
    }

  return delta;
}

template <class graph_struct>
void block_model<graph_struct>::move(const int v, const int t) {

  assert(new_bwd_edges.size() == 0);
  assert(new_fwd_edges.size() == 0);

  auto o{block[v]};
  block[v] = t;

  ++block_size[t];
  --block_size[o];

  for (auto j : model.nodes) {
    if (j == t) {
      if ((!prior_bwd_edge[j] and bwd_neighbors[j] > 0) or
          (!prior_fwd_edge[j] and fwd_neighbors[j] > 0)) {
        new_fwd_edges.push_back(t);
        new_fwd_edges.push_back(t);
      }
    } else {
      if (!prior_bwd_edge[j] and bwd_neighbors[j] > 0) {
        new_bwd_edges.push_back(j);
        new_bwd_edges.push_back(t);
      }
      if (!prior_fwd_edge[j] and fwd_neighbors[j] > 0) {
        new_fwd_edges.push_back(j);
        new_fwd_edges.push_back(t);
      }
    }

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
    // model.rename(o, model.nodes[model.size()-1]);
    model.rem_node(o);
  } else {
    for (auto e : model.successors[o])
      if (origin_fwd_edgecount[e.endpoint()] != e.weight())
        model.set_weight(o, e.endpoint(), e,
                         origin_fwd_edgecount[e.endpoint()]);
    for (auto e : model.predecessors[o])
      if (origin_bwd_edgecount[e.endpoint()] != e.weight())
        model.set_weight(e.endpoint(), o, e,
                         origin_bwd_edgecount[e.endpoint()]);
  }
}

template <class graph_struct>
void block_model<graph_struct>::compress(options &opt, std::mt19937 &rng) {

  epsilon = opt.epsilon;

  if (opt.checked()) {
    block_bits.resize(N);
    error_bits.resize(N);
    iblck_bits.resize(N);
    ierrr_bits.resize(N);
    bldlt_bits.resize(N);
    erdlt_bits.resize(N);
    for (auto i = 0; i < N; ++i) {			
      block_bits[i].resize(N);
      error_bits[i].resize(N);
      bldlt_bits[i].resize(N,0);
      erdlt_bits[i].resize(N,0);
      iblck_bits[i].resize(N);
      ierrr_bits[i].resize(N);
    }
  }

  intstack movable(data.capacity());
  movable.fill();

  if (opt.randomized()) {
    movable.shuffle(rng);
  }

  float nbits{0};
  if (opt.checked()) {
    nbits = get_objective(); // opt.alpha, opt.beta);
		for (auto i = 0; i < N; ++i) {	
		for (auto j = 0; j < N; ++j) {
			iblck_bits[i][j] = block_bits[i][j];
				ierrr_bits[i][j] = error_bits[i][j];
			}	
		}
    if (opt.verbosity >= block::options::NORMAL)
      std::cout << nbits << " / " << (model.size() * model.size()) << std::endl;
  }
  auto incr_nbits{nbits};

  auto prev{-1};

  float delta, update;

  while (true) {
    bool moved;
    auto move_vertex{-1};
    auto move_target{-1};
		update = -epsilon;

    // try to find a move of negative cost
    for (auto v : movable) {

      moved = false;
      for (auto t : model.nodes) {
        auto o{block[v]};
        if (o != t) {
          delta = get_cost(opt, v, t);

          if (delta < update) {
            move_vertex = v;
            move_target = t;
            update = delta;
            moved = true;
            if (opt.policy == block::options::FIRST)
              break;
          }
        }
      }

      if (moved)
        break;
    }

    if (move_vertex < 0)
      break;
    else {
			if (opt.policy == block::options::BEST)
				update = get_cost(opt, move_vertex, move_target);
			
      incr_nbits += update;
      if (opt.verbosity >= block::options::NORMAL)
        std::cout << std::setw(4) << model.size() << " " << std::setw(4)
                  << model.num_edges << " "
                  << "move " << move_vertex << ": " << block[move_vertex]
                  << " -> " << move_target << " (" << update;

      move(move_vertex, move_target);
      if (opt.stable)
        movable.remove(move_target);

      if (opt.checked()) {
        float pnbits{nbits};
        nbits = get_objective(); // opt.alpha, opt.beta);
        if (opt.verbosity >= block::options::NORMAL)
          std::cout << "/" << (nbits - pnbits);
      }
      if (opt.verbosity >= block::options::NORMAL)
        std::cout << ")" << std::endl;
      if (opt.checked())
        check(nbits, incr_nbits, opt.check_epsilon);
    }
  }
}

/*************** VERIFICATION STUFF ******************/

template <class graph_struct>
void block_model<graph_struct>::check(float nbits, float incr_nbits,
                                      float err_epsilon) {

  auto n{data.capacity()};
  auto m{model.size()};

  // check block_size
  util.clear();
  util.resize(n, 0);
  for (auto v : data.nodes) {
    ++util[block[v]];
  }
  for (auto u : model.nodes) {
    assert(block_size[u] == util[u]);
  }

  // check edge count
  auto total{0};
  util.clear();
  util.resize(m * m, 0);
  for (auto u : data.nodes) {
    auto i{model.nodes.index(block[u])};

    for (auto e : data[u]) {
      auto j{model.nodes.index(block[e.endpoint()])};
      ++util[i * m + j];
    }
  }
  for (auto u : model.nodes) {
    auto i{model.nodes.index(u)};
    for (auto e : model[u]) {
      auto j{model.nodes.index(e.endpoint())};

      if (e.weight() != util[i * m + j]) {
        std::cout << "problem with " << u << "," << e.endpoint() << " = "
                  << model.nodes[i] << "," << model.nodes[j] << " ("
                  << e.weight() << "/" << util[i * m + j] << "):\n"
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

  auto ok{std::abs(nbits - incr_nbits) <= err_epsilon};

  if (!ok) {
    std::cout << "global problem: " << std::setprecision(10) << nbits << " / "
              << incr_nbits << std::endl;
  }

  int ierr{-1};
  int jerr{-1};
  int terr{-1};

  float cibits{0}, cnbits{0};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      iblck_bits[i][j] += bldlt_bits[i][j];
      ierrr_bits[i][j] += erdlt_bits[i][j];
      cibits += iblck_bits[i][j];
      cibits += ierrr_bits[i][j];
      cnbits += block_bits[i][j];
      cnbits += error_bits[i][j];
      if (std::abs(block_bits[i][j] - iblck_bits[i][j]) > err_epsilon) {
        ierr = i;
        jerr = j;
        terr = 0;
        ok = false;
      }
      if (std::abs(error_bits[i][j] - ierrr_bits[i][j]) > err_epsilon) {
        ierr = i;
        jerr = j;
        terr = 1;
        ok = false;
      }
    }
  }

  if (std::abs(nbits - cnbits) > err_epsilon) {
    std::cout << "BF count: " << cnbits << " / " << nbits << std::endl;
  }

  if (std::abs(incr_nbits - cibits) > err_epsilon) {
    std::cout << "INCR count: " << cibits << " / " << incr_nbits << std::endl;
  }

  if (!ok) {
    std::cout << "\nDELTA MODEL SIZE:\n";
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << bldlt_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nINCR MODEL SIZE:\n";
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << iblck_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nREAL MODEL SIZE:\n";
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << block_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }

    std::cout << "\nDELTA ERROR SIZE:\n";
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << erdlt_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nINCR ERROR SIZE:\n";
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << ierrr_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nREAL ERROR SIZE:\n";
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << std::setw(5) << std::setprecision(2) << error_bits[i][j]
                  << " ";
      }
      std::cout << std::endl;
    }

    std::cout << (terr ? "error" : "block") << " bits discrepancy on " << ierr
              << "," << jerr << std::endl;

    exit(1);
  }

  assert(total == data.num_edges);
}

template <class graph_struct>
float block_model<graph_struct>::get_objective() { // const float A, const float
                                                   // B) {
  float nbits{0};

  for (int i = 0; i < N; ++i) {
    std::fill(begin(block_bits[i]), end(block_bits[i]), 0);
    std::fill(begin(error_bits[i]), end(error_bits[i]), 0);
  }

  for (auto u : model.nodes) {

    for (auto v : model.nodes) {

      assert((block_size[u] * block_size[v]) > 0);

      auto block_sz{
          // A *
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

      auto error_sz{// B *
                    size * ((density - 1) * std::log2(1 - density) -
                            density * std::log2(density))};
      nbits += error_sz;

      error_bits[u][v] += error_sz;
    }
  }

  return nbits;
}

} // namespace block

#endif
