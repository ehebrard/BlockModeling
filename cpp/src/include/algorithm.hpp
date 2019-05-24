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

template <class graph_struct> void compress(graph_struct &forward) {

  // reversed graph
  dgraph backward(forward, true);

  std::cout << backward << std::endl;

  // block model
  dgraph fwd_block_model(forward);

  dgraph bwd_block_model(backward);

  //   std::cout << fwd_block_model << std::endl;
  // std::cout << bwd_block_model << std::endl;

  auto N{forward.size()};

  double objective{static_cast<double>(N * N)};

  // std::vector<int> origin_fwd_neighbors(N);
  // std::vector<int> origin_bwd_neighbors(N);
  // std::vector<int> target_fwd_neighbors(N);
  // std::vector<int> target_bwd_neighbors(N);
  std::vector<int> fwd_neighbors(N);
  std::vector<int> bwd_neighbors(N);

  std::vector<wtype> origin_fwd_density(N);
  std::vector<wtype> target_fwd_density(N);
  std::vector<wtype> origin_bwd_density(N);
  std::vector<wtype> target_bwd_density(N);

  std::vector<int> color(N);
  for (auto i{0}; i < N; ++i)
    color[i] = i;

  while (true) {
    // try to find a move of negative cost
    for (auto v : forward.nodes) {
      for (auto i : fwd_block_model.nodes)
        if (color[v] != i) {
          std::cout << "check " << v << " -> b" << i << std::endl;

          std::fill(begin(fwd_neighbors), end(fwd_neighbors), 0);
          std::fill(begin(bwd_neighbors), end(bwd_neighbors), 0);

          std::fill(begin(origin_fwd_density), end(origin_fwd_density), 0);
          std::fill(begin(target_fwd_density), end(target_fwd_density), 0);
          std::fill(begin(origin_bwd_density), end(origin_bwd_density), 0);
          std::fill(begin(target_bwd_density), end(target_bwd_density), 0);

          for (auto e : fwd_block_model[color[v]])
            origin_fwd_density[e.endpoint()] = e.weight();
          for (auto e : fwd_block_model[i])
            target_fwd_density[e.endpoint()] = e.weight();
          for (auto e : bwd_block_model[color[v]])
            origin_bwd_density[e.endpoint()] = e.weight();
          for (auto e : bwd_block_model[i])
            target_bwd_density[e.endpoint()] = e.weight();

          for (auto e : forward[v]) {
            ++fwd_neighbors[color[e.endpoint()]];
          }

          for (auto e : backward[v]) {
            ++bwd_neighbors[color[e.endpoint()]];
          }
					
					std::cout << "     ";
					for(auto j : fwd_block_model.nodes)
					{
						std::cout << std::setw(3) << j << " ";
					}
					std::cout << std::endl << "d_oj ";
					for(auto j : fwd_block_model.nodes)
						std::cout << std::setw(3) << std::setprecision(1) << origin_fwd_density[j] << " ";
					std::cout << std::endl << "d_jo ";
					for(auto j : fwd_block_model.nodes)
						std::cout << std::setw(3) << std::setprecision(1) << origin_bwd_density[j] << " ";
					std::cout << std::endl << "d_tj ";
					for(auto j : fwd_block_model.nodes)
						std::cout << std::setw(3) << std::setprecision(1) << target_fwd_density[j] << " ";
					std::cout << std::endl << "d_jt ";
					for(auto j : fwd_block_model.nodes)
						std::cout << std::setw(3) << std::setprecision(1) << target_bwd_density[j] << " ";
					std::cout << std::endl << "N+";
					for(auto j : fwd_block_model.nodes)
						std::cout << std::setw(3) << std::setprecision(1) << target_bwd_density[j] << " ";
					std::cout << std::endl;
					
					
					
        }
    }

    break;
  }
}



} // namespace block

#endif
