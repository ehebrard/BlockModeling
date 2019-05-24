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



template <class graph_struct>
void compress(graph_struct &in)
{

	// reversed graph
	dgraph out(in, true);

	// block model
	dgraph block_model(in);
	
	
	std::cout << block_model << std::endl;
	
	
}



} // namespace block

#endif
