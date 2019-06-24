#ifndef __BLOCK__ALGORITHM_HH
#define __BLOCK__ALGORITHM_HH

#include <vector>

#include "global.hpp"
#include "intstack.hpp"
#include "options.hpp"
#include "partition.hpp"
#include "timing.hpp"

namespace block
{

template <class graph_struct>
void BFS(graph_struct &g, intstack &order, int start = 0) {
  for (auto i{start}; i < order.size(); ++i)
    for (auto e : g[order[i]])
      order.add(e.endpoint());
}

double get_error(const int m, const int size) {

  double delta{0};

  double dij = static_cast<double>(m) / static_cast<double>(size);
  if (m > 0 and m < size)
    delta = static_cast<double>(size) *
            ((dij - 1.0) * std::log2(1.0 - dij) - dij * std::log2(dij));

  return delta;
}

// compute the error delta of a move w.r.t. blocks i -> j
// there are 'mij' edges from i to j and 'deltaij' are added
double get_error_delta(const int mij, const int psize, const int nsize,
                       const int deltaij) {

  auto delta{get_error(mij + deltaij, nsize)};
  delta -= get_error(mij, psize);

  return delta;
}

template <class graph_struct> class block_model {

public:
  block_model(graph_struct &g, options &opt)
      : opt(opt), data(g), model(g), N{g.capacity()}, fwd_neighbors(N),
        bwd_neighbors(N), origin_fwd_edgecount(N), target_fwd_edgecount(N),
        origin_bwd_edgecount(N), target_bwd_edgecount(N), prior_fwd_edges(N),
        prior_bwd_edges(N), best_move_target(N, -1), best_move_delta(N, 0), sample(N),
        forbidden_move(N, -1), num_cost{0} {

    bag.resize(N, N);
		
		sample.fill();

    int i{0};
    for (auto u : g.nodes) {
      bag.add_elt(u, i++);
    }
  }

  template <class block_struct>
  block_model(graph_struct &g, block_struct &blocks, options &opt)
      : opt(opt), data(g), N{blocks.size()}, fwd_neighbors(N), bwd_neighbors(N),
        origin_fwd_edgecount(N), target_fwd_edgecount(N),
        origin_bwd_edgecount(N), target_bwd_edgecount(N), prior_fwd_edges(N),
        prior_bwd_edges(N), best_move_target(g.capacity(), -1),
        best_move_delta(g.capacity(), 0), sample(N), forbidden_move(g.capacity(), -1), num_cost{0} {

    model.initialise(N);

    bag.resize(g.capacity(), N);
		
		sample.fill();

    int i{0};
    for (auto b{begin(blocks)}; b != end(blocks); ++b) {
      for (auto u{b->begin()}; u != b->end(); ++u) {
        bag.add_elt(*u, i);
      }
      ++i;
    }

    for (auto b{begin(blocks)}; b != end(blocks); ++b) {
      std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
      for (auto u{b->begin()}; u != b->end(); ++u) {
        for (auto e : g.successors[*u]) {
          ++fwd_neighbors[get_block(e.endpoint())];
        }
      }

      auto bi{(b - begin(blocks))};
      for (int bj{0}; bj < N; ++bj) {
        if (fwd_neighbors[bj] > 0) {
          model.add_edge(bi, bj, fwd_neighbors[bj]);
        }
      }
    }
  }

  void get_blocks(std::vector<std::vector<int>> &blocks);

  void compress(std::mt19937 &rng);

  double
  get_move_cost(const int v,
                const int t); // cost of moving vertex v from its block to t
  double get_merge_cost(const int o,
                        const int t); // cost of merging blocks o and t
  void move(const int v, const int t);

  options &opt;

  double epsilon{1e-3};

  graph_struct &data;
  graph_struct model;

  double get_objective(); // const double A, const double B);

private:
  size_t N;

  std::vector<int> fwd_neighbors;
  std::vector<int> bwd_neighbors;

  std::vector<wtype> origin_fwd_edgecount;
  std::vector<wtype> target_fwd_edgecount;
  std::vector<wtype> origin_bwd_edgecount;
  std::vector<wtype> target_bwd_edgecount;

  intstack prior_fwd_edges;
  intstack prior_bwd_edges;

  partition bag;
  inline int blocksize(const int b) { return bag[model.nodes.index(b)].size(); }
  inline int get_block(const int u) { return model.nodes[bag.bag_of(u)]; }

  void assign_block(const int v, const int b);

  std::vector<int> new_fwd_edges;
  std::vector<int> new_bwd_edges;

  /********* BEST MOVE STUFF ***********/
  std::vector<int> best_move_origin;
  std::vector<int> best_move_target;
  std::vector<int> forbidden_move;
  std::vector<double> best_move_delta;

  void store_best_move(const int v);
	void sample_best_move(const int v, const int size, std::mt19937 &rng);
	intstack sample;
	

  /********* VERIFICATION STUFF ***********/
  // for the brute-force thing
  // std::vector<double> bitsize;
  std::vector<int> util;
  // double err_epsilon{.5};
  void check(double nbits, double incr_nbits, double eps);
  std::vector<std::vector<double>> block_bits;
  std::vector<std::vector<double>> error_bits;
  std::vector<std::vector<double>> bldlt_bits;
  std::vector<std::vector<double>> erdlt_bits;
  std::vector<std::vector<double>> iblck_bits;
  std::vector<std::vector<double>> ierrr_bits;

  double init_check_structs();

  long long int num_cost;
};

template <class graph_struct>
void block_model<graph_struct>::get_blocks(
    std::vector<std::vector<int>> &blocks) {
  blocks.clear();
  blocks.resize(model.size());
  for (auto x : data.nodes)
    blocks[bag.bag_of(x)].push_back(x);
  // blocks[model.nodes.index(block[x])].push_back(x);
}

template <class graph_struct>
void block_model<graph_struct>::assign_block(const int v, const int t) {

  auto o{get_block(v)};
  auto oi{model.nodes.index(o)};
  auto ti{model.nodes.index(t)};

  // assert(bag.bag_of(v)==oi);
  bag.move(v, ti);

  if (blocksize(o) == 0) {
    model.rem_node(o);
    bag.remove_bag(oi);
		sample.remove(o);
  }

  forbidden_move[v] = o;
}

template <class graph_struct>
double block_model<graph_struct>::get_merge_cost(const int o, const int t) {

  double delta{0};

  // + loop (ot) - loop o - loop t - arcs ot - arcs to
  double ijdelta{static_cast<double>(
      (blocksize(o) + blocksize(t)) * (blocksize(o) + blocksize(t)) + 1)};
  ijdelta -= static_cast<double>((blocksize(o) * blocksize(o)) + 1);
  delta += ijdelta;
  if (opt.checked())
    bldlt_bits[o][o] = ijdelta;

  ijdelta = -static_cast<double>((blocksize(t) * blocksize(t)) + 1);
  delta += ijdelta;
  if (opt.checked())
    bldlt_bits[t][t] = ijdelta;

  ijdelta = -static_cast<double>((blocksize(o) * blocksize(t)) + 1);
  delta += 2 * ijdelta;
  if (opt.checked())
    bldlt_bits[o][t] = bldlt_bits[t][o] = ijdelta;

  // we need to compute error delta only with respect to blocks with at least
  // one edge with o or t
  for (auto e : model.successors[o])
    origin_fwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.predecessors[o])
    origin_bwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.successors[t])
    target_fwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.successors[o])
    target_bwd_edgecount[e.endpoint()] = e.weight();

  for (auto j : model.nodes) {
    if (j == o) {

      // error on oo
      ijdelta = get_error(origin_fwd_edgecount[o] + origin_fwd_edgecount[t],
                          (blocksize(o) + blocksize(t)) *
                              (blocksize(o) + blocksize(t)));
      -get_error(origin_fwd_edgecount[o], blocksize(o) * blocksize(o));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[o][o] = ijdelta;

      // error on ot
      ijdelta - get_error(origin_fwd_edgecount[t], blocksize(o) * blocksize(t));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[o][t] = ijdelta;

      // error on ot
      ijdelta - get_error(origin_bwd_edgecount[t], blocksize(o) * blocksize(t));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[t][o] = ijdelta;

    } else if (j == t) {

      // error on tt
      ijdelta =
          -get_error(target_fwd_edgecount[t], (blocksize(t) * blocksize(t)));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[t][t] = ijdelta;

      // // error on ot [already counted]
      // ijdelta = -get_error(target_fwd_edgecount[t], (blocksize(t) *
      // blocksize(t)));
      // 		 	delta += ijdelta;
      // 		 	if(opt.checked())
      // 		 		erdlt_bits[o][t] = ijdelta;

      // // error on to [already counted]
      // ijdelta = -get_error(target_bwd_edgecount[t], (blocksize(t) *
      // blocksize(t)));
      // 		 	delta += ijdelta;
      // 		 	if(opt.checked())
      // 		 		erdlt_bits[t][o] = ijdelta;

    } else {
      // remove the edges jt and tj
      ijdelta = static_cast<double>((blocksize(j) * blocksize(t)) + 1);
      delta += 2 * ijdelta;
      if (opt.checked())
        bldlt_bits[j][t] = bldlt_bits[t][j] = ijdelta;

      // add the edges j(ot) and (ot)j remove the edges jo and oj
      ijdelta = (static_cast<double>(
                     (blocksize(j) * (blocksize(o) + blocksize(t))) + 1) -
                 static_cast<double>((blocksize(j) * blocksize(o)) + 1));
      delta += 2 * ijdelta;
      if (opt.checked())
        bldlt_bits[j][o] = bldlt_bits[o][j] = ijdelta;

      // error on oj
      ijdelta = get_error(origin_fwd_edgecount[j] + target_fwd_edgecount[j],
                          (blocksize(o) + blocksize(t)) * blocksize(j));
      -get_error(origin_fwd_edgecount[j], blocksize(o) * blocksize(j));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[o][j] = ijdelta;

      // error on tj
      ijdelta = get_error(target_fwd_edgecount[j] + target_fwd_edgecount[j],
                          (blocksize(o) + blocksize(t)) * blocksize(j));
      -get_error(target_fwd_edgecount[j], blocksize(t) * blocksize(j));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[t][j] = ijdelta;

      // error on jo
      ijdelta = get_error(origin_bwd_edgecount[j] + target_bwd_edgecount[j],
                          (blocksize(o) + blocksize(t)) * blocksize(j));
      -get_error(origin_bwd_edgecount[j], blocksize(o) * blocksize(j));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[j][o] = ijdelta;

      // error on tj
      ijdelta = get_error(target_bwd_edgecount[j] + target_bwd_edgecount[j],
                          (blocksize(o) + blocksize(t)) * blocksize(j));
      -get_error(target_bwd_edgecount[j], blocksize(t) * blocksize(j));
      delta += ijdelta;
      if (opt.checked())
        erdlt_bits[j][t] = ijdelta;
    }
  }

  // // error on t -> t
  // err = // B *
  //     get_error_delta(target_bwd_edgecount[t], blocksize(t) * blocksize(t),
  //                     (blocksize(t) + 1) * (blocksize(t) + 1),
  //                     bwd_neighbors[t] + fwd_neighbors[t]);
  // ijdelta += err;
  // if (opt.checked())
  //   erdlt_bits[t][t] += err;
  //
  // // error on t -> o
  // err = // B *
  //     get_error_delta(origin_bwd_edgecount[t], blocksize(t) * blocksize(o),
  //                     (blocksize(t) + 1) * (blocksize(o) - 1),
  //                     fwd_neighbors[o] - bwd_neighbors[t]);
  // ijdelta += err;
  // if (opt.checked())
  //   erdlt_bits[t][o] += err;
  //
  // assert(origin_bwd_edgecount[o] == origin_fwd_edgecount[o]);
  //
  //
  //
  // // error on o -> t
  // err = // B *
  //     get_error_delta(target_bwd_edgecount[o], blocksize(t) * blocksize(o),
  //                     (blocksize(t) + 1) * (blocksize(o) - 1),
  //                     bwd_neighbors[o] - fwd_neighbors[t]);

  //
  // //
  // for(auto j : model.nodes)
  // {
  // 	// update the
  // }

  return delta;
}

template <class graph_struct>
double block_model<graph_struct>::get_move_cost(const int v, const int t) {

	// std::cout << " compute for " << v << " -> " << t << std::endl;

  ++num_cost;

  // double A{opt.alpha};
  // double B{opt.beta};

  if (opt.checked()) {
    for (int i = 0; i < N; ++i) {
      std::fill(begin(bldlt_bits[i]), end(bldlt_bits[i]), 0);
      std::fill(begin(erdlt_bits[i]), end(erdlt_bits[i]), 0);
    }
  }

  double delta{0}, ijdelta{0};

  int o{get_block(v)};

  prior_fwd_edges.clear();
  prior_bwd_edges.clear();

  for (auto u : model.nodes) {
    fwd_neighbors[u] = bwd_neighbors[u] = origin_fwd_edgecount[u] =
        origin_bwd_edgecount[u] = target_fwd_edgecount[u] =
            target_bwd_edgecount[u] = 0;
	}

  // std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
  // std::fill(begin(bwd_neighbors), end(bwd_neighbors), 0);
  //
  //
  // std::fill(begin(origin_fwd_edgecount), end(origin_fwd_edgecount), 0);
  // std::fill(begin(target_fwd_edgecount), end(target_fwd_edgecount), 0);
  // std::fill(begin(origin_bwd_edgecount), end(origin_bwd_edgecount), 0);
  // std::fill(begin(target_bwd_edgecount), end(target_bwd_edgecount), 0);

	
	auto u{0};
  for (auto e : model.successors[o])
    origin_fwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.successors[t]) {
		u = e.endpoint();
    target_fwd_edgecount[u] = e.weight();
		prior_fwd_edges.add(u);
  }
  for (auto e : model.predecessors[o])
    origin_bwd_edgecount[e.endpoint()] = e.weight();
  for (auto e : model.predecessors[t]) {
		u = e.endpoint();
    target_bwd_edgecount[u] = e.weight();
		prior_bwd_edges.add(u);
  }
	
	new_fwd_edges.clear();
	new_bwd_edges.clear();
	
	
  for (auto e : data.successors[v]) {
    u = get_block(e.endpoint());
    if (!prior_fwd_edges.contain(u) and fwd_neighbors[u] == 0) {
      new_fwd_edges.push_back(u);
		}
    ++fwd_neighbors[u];
  }

  for (auto e : data.predecessors[v]) {
    u = get_block(e.endpoint());
    if (!prior_bwd_edges.contain(u) and bwd_neighbors[u] == 0 and
        (u != t or fwd_neighbors[u] == 0)) {
      new_bwd_edges.push_back(u);
		}
    ++bwd_neighbors[u];
  }


  assert(target_bwd_edgecount[o] == origin_fwd_edgecount[t]);
  assert(target_fwd_edgecount[o] == origin_bwd_edgecount[t]);

  double err;

  auto ojdelta = 0.0;
  auto tjdelta = 0.0;

  ijdelta = 0.0;

  // o.o
  if (blocksize(o) > 1)
    ojdelta += std::log2(
        static_cast<double>((blocksize(o) - 1) * (blocksize(o) - 1) + 1));
  ojdelta -= std::log2(static_cast<double>(blocksize(o) * blocksize(o) + 1));
  // ojdelta *= A;
  if (opt.checked())
    bldlt_bits[o][o] += ojdelta;

  // t.t
  // if (nsizet > 0)
  tjdelta += std::log2(
      static_cast<double>((blocksize(t) + 1) * (blocksize(t) + 1) + 1));
  tjdelta -= std::log2(static_cast<double>(blocksize(t) * blocksize(t) + 1));
  // tjdelta *= A;
  if (opt.checked())
    bldlt_bits[t][t] += tjdelta;

  // o.t and t.o
  if (blocksize(o) > 1)
    ijdelta += std::log2(
        static_cast<double>((blocksize(o) - 1) * (blocksize(t) + 1) + 1));
  ijdelta -= std::log2(static_cast<double>(blocksize(o) * blocksize(t) + 1));
  // ijdelta *= A;
  if (opt.checked()) {
    bldlt_bits[t][o] += ijdelta;
    bldlt_bits[o][t] += ijdelta;
  }
  ijdelta *= 2;

  ijdelta += (ojdelta + tjdelta);

  // error on t -> t
  err = // B *
      get_error_delta(target_bwd_edgecount[t], blocksize(t) * blocksize(t),
                      (blocksize(t) + 1) * (blocksize(t) + 1),
                      bwd_neighbors[t] + fwd_neighbors[t]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[t][t] += err;

  // error on t -> o
  err = // B *
      get_error_delta(origin_bwd_edgecount[t], blocksize(t) * blocksize(o),
                      (blocksize(t) + 1) * (blocksize(o) - 1),
                      fwd_neighbors[o] - bwd_neighbors[t]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[t][o] += err;

  assert(origin_bwd_edgecount[o] == origin_fwd_edgecount[o]);

  // error_delta on o -> o
  err = // B *
      get_error_delta(origin_bwd_edgecount[o], blocksize(o) * blocksize(o),
                      (blocksize(o) - 1) * (blocksize(o) - 1),
                      -bwd_neighbors[o] - fwd_neighbors[o]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[o][o] += err;

  // error on o -> t
  err = // B *
      get_error_delta(target_bwd_edgecount[o], blocksize(t) * blocksize(o),
                      (blocksize(t) + 1) * (blocksize(o) - 1),
                      bwd_neighbors[o] - fwd_neighbors[t]);
  ijdelta += err;
  if (opt.checked())
    erdlt_bits[o][t] += err;

  delta += ijdelta;

  // update the densities
  for (auto j : model.nodes)
    if (j != o and j != t) {

      auto psizeo{blocksize(o) * blocksize(j)};
      auto psizet{blocksize(t) * blocksize(j)};
      auto nsizeo{(blocksize(o) - 1) * blocksize(j)};
      auto nsizet{(blocksize(t) + 1) * blocksize(j)};

      ojdelta = 0.0;
      if (nsizeo > 0)
        ojdelta += std::log2(static_cast<double>(nsizeo + 1));
      ojdelta -= std::log2(static_cast<double>(psizeo + 1));
      ojdelta *= 2;

      tjdelta = 0.0;
      if (nsizet > 0)
        tjdelta += std::log2(static_cast<double>(nsizet + 1));
      tjdelta -= std::log2(static_cast<double>(psizet + 1));
      tjdelta *= 2;

      ijdelta = (ojdelta + tjdelta);

      if (opt.checked()) {
        err = std::log2(static_cast<double>(nsizeo + 1)) -
              std::log2(static_cast<double>(psizeo + 1));
        bldlt_bits[j][o] += err;
        bldlt_bits[o][j] += err;

        err = std::log2(static_cast<double>(nsizet + 1)) -
              std::log2(static_cast<double>(psizet + 1));
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

  auto o{get_block(v)};

  assign_block(v, t);

  for (auto j : model.nodes) {
		
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
    auto j{new_bwd_edges.back()};
    new_bwd_edges.pop_back();
    model.add_edge(j, t, target_bwd_edgecount[j]);
  }
  while (new_fwd_edges.size() > 0) {
    auto j{new_fwd_edges.back()};
    new_fwd_edges.pop_back();
    model.add_edge(t, j, target_fwd_edgecount[j]);
  }

  for (auto e : model.successors[t])
    if (target_fwd_edgecount[e.endpoint()] != e.weight())
      model.set_weight(t, e.endpoint(), e, target_fwd_edgecount[e.endpoint()]);

  for (auto e : model.predecessors[t])
    if (target_bwd_edgecount[e.endpoint()] != e.weight())
      model.set_weight(e.endpoint(), t, e, target_bwd_edgecount[e.endpoint()]);

  if (blocksize(o) == 0) {
    // model.rename(o, model.nodes[model.size()-1]);
    // model.rem_node(o);
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


// compute (brute force) the best move with u and store it
template <class graph_struct>
void block_model<graph_struct>::store_best_move(const int v) {

  best_move_delta[v] = static_cast<double>(N*N);
  double delta{0};

  auto o{get_block(v)};

  for (auto t : model.nodes) {

    if (o != t and t != forbidden_move[v]) {
      delta = get_move_cost(v, t);

      if (opt.verbosity >= options::YACKING)
        std::cout << "  probe move " << v << ": " << get_block(v) << " -> " << t
                  << " (" << delta;

      if (best_move_delta[v] > delta) {
        best_move_delta[v] = delta;
        best_move_target[v] = t;

        if (opt.verbosity >= options::YACKING)
          std::cout << "*";
      }

      if (opt.verbosity >= options::YACKING)
        std::cout << ")\n";
    }
  }

  if (opt.verbosity >= options::YACKING)
    std::cout << " ==> " << best_move_target[v] << " (" << best_move_delta[v]
              << ")\n";
}

// compute (brute force) the best move with u and store it
template <class graph_struct>
void block_model<graph_struct>::sample_best_move(const int v, const int size, std::mt19937 &rng) {

  best_move_delta[v] = static_cast<double>(N*N);
  double delta{0};

  auto o{get_block(v)};
	
	for (auto i{0}; i<size; ++i)
	{
		auto j{rng() % (sample.size() - i)};
		sample.move(sample[i+j], i);
	}
	
	for (auto i{0}; i<size; ++i)
	{
		auto t{sample[i]};
		if (o != t and t != forbidden_move[v]) {
      delta = get_move_cost(v, t);

      if (opt.verbosity >= options::YACKING)
        std::cout << "  probe move " << v << ": " << get_block(v) << " -> " << t
                  << " (" << delta;

      if (best_move_delta[v] > delta) {
        best_move_delta[v] = delta;
        best_move_target[v] = t;

        if (opt.verbosity >= options::YACKING)
          std::cout << "*";
      }

      if (opt.verbosity >= options::YACKING)
        std::cout << ")\n";
		}
	}

  if (opt.verbosity >= options::YACKING)
    std::cout << " ==> " << best_move_target[v] << " (" << best_move_delta[v]
              << ")\n";
}


template <class graph_struct>
double block_model<graph_struct>::init_check_structs() {
block_bits.resize(N);
error_bits.resize(N);
iblck_bits.resize(N);
ierrr_bits.resize(N);
bldlt_bits.resize(N);
erdlt_bits.resize(N);
for (auto i = 0; i < N; ++i) {
  block_bits[i].resize(N);
  error_bits[i].resize(N);
  bldlt_bits[i].resize(N, 0);
  erdlt_bits[i].resize(N, 0);
  iblck_bits[i].resize(N);
  ierrr_bits[i].resize(N);
}

double nbits = get_objective(); // opt.alpha, opt.beta);
for (auto i = 0; i < N; ++i) {	
for (auto j = 0; j < N; ++j) {
  iblck_bits[i][j] = block_bits[i][j];
  ierrr_bits[i][j] = error_bits[i][j];
}
}
if (opt.verbosity >= block::options::NORMAL)
  std::cout << nbits << " / " << (model.size() * model.size()) << std::endl;

return nbits;
}

template <class graph_struct>
void block_model<graph_struct>::compress(std::mt19937 &rng) {

  epsilon = opt.epsilon;

  double mean_delta{0};

  double nbits{0}, incr_nbits{0}, update, delta;
  if (opt.checked())
    incr_nbits = nbits = init_check_structs();
	else
		incr_nbits = nbits = get_objective();

	std::vector<int> node;
	for(auto v : data.nodes)
		node.push_back(v);
	std::shuffle(begin(node), end(node), rng);
	
  best_move_origin.clear();
  // auto current{0};

  int iter{0}, last_improvement{-1};
  double slack{10};

  while (iter - last_improvement < 2 * node.size()) {
    auto v{node[iter % node.size()]};

		if(model.size() > 100)
			sample_best_move(v, 100, rng);
		else
			store_best_move(v);

    mean_delta += best_move_delta[v];
    // std::cout << best_move_delta[v] << " " << (mean_delta / (iter+1)) <<
    // std::endl;

    if ((best_move_delta[v] < -epsilon) or
        (best_move_delta[v] < slack and
         (iter - last_improvement) >= node.size())) {
      auto delta{get_move_cost(v, best_move_target[v])};

      incr_nbits += delta;
      if (opt.verbosity >= block::options::NORMAL)
        std::cout << std::left << std::setw(7)
                  << (2 * node.size() - iter + last_improvement) << " "
                  << std::setw(9) << std::setprecision(4)
                  << (mean_delta / (iter + 1)) << " " << std::setw(8)
                  << std::setprecision(3) << cpuTime() << " " << std::setw(10)
                  << num_cost << " | " << std::setw(5) << std::right
                  << model.size() << "/" << std::left << std::setw(4)
                  << model.num_edges << std::right << " " << std::setw(8)
                  << (int)incr_nbits << "  -- "
                  << "move " << v << ": " << get_block(v) << " -> "
                  << best_move_target[v] << " (" << delta;

      if (best_move_delta[v] < -epsilon)
        last_improvement = iter;

      move(v, best_move_target[v]);

      if (opt.checked()) {
        double pnbits{nbits};
        nbits = get_objective(); // opt.alpha, opt.beta);
        if (opt.verbosity >= block::options::NORMAL)
          std::cout << "/" << (nbits - pnbits);
        check(nbits, incr_nbits, opt.check_epsilon);
      }

      if (opt.verbosity >= block::options::NORMAL)
        std::cout << ")" << std::endl;
    }

    ++iter;
  }

  // std::fill(begin(forbidden_move), end(forbidden_move), -1)
  // {
  //
  // }
}

/*************** VERIFICATION STUFF ******************/

template <class graph_struct>
void block_model<graph_struct>::check(double nbits, double incr_nbits,
                                      double err_epsilon) {

  auto n{data.capacity()};
  auto m{model.size()};

  // check block_size
  util.clear();
  util.resize(n, 0);
  for (auto v : data.nodes) {
    ++util[get_block(v)];
  }
  for (auto u : model.nodes) {
    assert(blocksize(u) == util[u]);
  }

  // check edge count
  auto total{0};
  util.clear();
  util.resize(m * m, 0);
  for (auto u : data.nodes) {
    auto i{bag.bag_of(u)}; // model.nodes.index(block[u])};
    // assert(i == bag.bag_of(u));

    for (auto e : data[u]) {
      auto j{
          bag.bag_of(e.endpoint())}; // model.nodes.index(block[e.endpoint()])};
      // assert(j==bag.bag_of(e.endpoint()));
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

  double cibits{0}, cnbits{0};
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
double block_model<graph_struct>::get_objective() { // const double A, const double
                                                    // B) {
  double nbits{0};

  if (!opt.checked() and model.size() == data.size())
    return static_cast<double>(model.size() * model.size());

  if (opt.checked())
    for (int i = 0; i < N; ++i) {
      std::fill(begin(block_bits[i]), end(block_bits[i]), 0);
      std::fill(begin(error_bits[i]), end(error_bits[i]), 0);
    }

  for (auto u : model.nodes) {

    for (auto v : model.nodes) {

      assert((blocksize(u) * blocksize(v)) > 0);

      auto block_sz{
          // A *
          std::log2(static_cast<double>(blocksize(u) * blocksize(v) + 1))};

      nbits += block_sz;

      if (opt.checked())
        block_bits[u][v] = block_sz;

      // std::cout << " + " << block_sz << " = " << nbits << std::endl;
    }

    // error size
    for (auto e : model.successors[u]) {
      auto v{e.vertex};
      auto nedge{e.weight()};
      auto size{blocksize(u) * blocksize(v)};

      if (nedge == 0 or nedge == size)
        continue;

      double fnedge{static_cast<double>(nedge)};
      double fsize{static_cast<double>(size)};
      double density{fnedge / fsize};
			
			

      auto error_sz{// B *
                    fsize * ((density - 1) * std::log2(1 - density) -
                            density * std::log2(density))};
										
										
										// if(u==1 and v==13)
										// {
										// 	std::cout << "HERE: " << density << " " << size << " => " << error_sz << std::endl;
										// }
										
      nbits += error_sz;

      if (opt.checked())
        error_bits[u][v] = error_sz;
    }
  }
  // return 0;

  return nbits;
}

} // namespace block

#endif
