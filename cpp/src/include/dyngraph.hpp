#ifndef __BLOCK__DYNGRAPH_HH
#define __BLOCK__DYNGRAPH_HH

#include <vector>
#include <assert.h>

#include "intstack.hpp"

// #define _VERIFY_MCGRAPH 1

namespace block
{

template <class WEIGHT> struct edge_info {
  // ranks of the vertices in eachother neighbors list
  int rank[2];
  // weight of the edge
  WEIGHT weight;

  edge_info(const int a, const int b, const WEIGHT w) {
    // rank of origin vertex in in_neighbor list of the target
    rank[0] = a;
    // rank of target vertex in out_neighbor list of the origin
    rank[1] = b;
    weight = w;
  }

  int operator[](const int i) const { return rank[i]; }
  int &operator[](const int i) { return rank[i]; }
};

template <class W, class T>
std::ostream &operator<<(std::ostream &os, const edge_info<W> &x) {
  return x.describe(os);
}

template <class WEIGHT> struct neighbor_info {

  WEIGHT weight() const { return edge->weight; }
  void set_weight(WEIGHT w) const { edge->weight = w; }
  int endpoint() const { return vertex; }

  int vertex;

  edge_info<WEIGHT> *edge;

  neighbor_info(const int v, edge_info<WEIGHT> *e) : vertex{v}, edge(e) {}
  neighbor_info(const neighbor_info<WEIGHT> &ni, edge_info<WEIGHT> *e)
      : vertex{ni.vertex}, edge(e) {}

  int get_origin_rank() const { return edge->rank[0]; }
  int get_target_rank() const { return edge->rank[1]; }

  // set the rank of u in vertex's neighbor list
  void set_origin_rank(const int r) const { edge->rank[0] = r; }
  // set the rank of vertex in u's neighbor list
  void set_target_rank(const int r) const { edge->rank[1] = r; }
};

template <class WEIGHT> class dyngraph {

public:
  /** nodes */
  // the list of nodes of the graph, in no particular order
  intstack nodes;

  /** edges */
  // the current number of edges (the list of edges is not updated)
  size_t num_edges;

  // ranks stores the indices of the nodes in eachother's neighbor list so
  // that we can remove the edge efficiently
  std::vector<edge_info<WEIGHT>> edges;

  /** neighborhood */
  // the list of successors of each node
  std::vector<std::vector<neighbor_info<WEIGHT>>> successors;
  // the list of predecessors of each node
  std::vector<std::vector<neighbor_info<WEIGHT>>> predecessors;

  // whether the vertex has a loop
  std::vector<int> loop;

  // the weights on the nodes
  std::vector<WEIGHT> weight;
  // the total weight of the successors
  std::vector<WEIGHT> successors_weight;
  // the total weight of the predecessors
  std::vector<WEIGHT> predecessors_weight;

  // total node weight
  WEIGHT node_weight;
  // total edge weight
  WEIGHT edge_weight;


  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::iterator begin() {
    return begin(successors);
  }
  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::reverse_iterator
  rbegin() {
    return rbegin(predecessors);
  }

  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::iterator end() {
    return end(successors);
  }
  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::reverse_iterator
  rend() {
    return rend(predecessors);
  }

  std::vector<neighbor_info<WEIGHT>> &operator[](const int v) {
    return successors[v];
  }

  dyngraph() {}
  dyngraph(const int n) { initialise(n); }

  template <class graph_struct> dyngraph(const graph_struct &g) {
    deep_copy(g);
  }

  void initialise(const int n) {
    weight.resize(n, 0);

    nodes.reserve(n);
    nodes.fill();
		
		edges.reserve(std::min(n*n,100*n));

    successors.resize(n);
    predecessors.resize(n);

    // loop.resize(n,0);

    successors_weight.resize(n, 0);
    predecessors_weight.resize(n, 0);

    num_edges = 0;
    node_weight = 0;
    edge_weight = 0;
  }

  void reinit() {
    for (auto i{nodes.size()}; i < nodes.capacity(); ++i)
      add_node(nodes[i]);
  }

  void undo_add() { rem_node(nodes.back()); }

  void undo_rem() { add_node(nodes[nodes.size()]); }

  template <class graph_struct> void deep_copy(const graph_struct &g) {
    weight = g.weight;
    nodes.reserve(g.capacity());
    if (g.size() == g.capacity())
      nodes.fill();
    else
      for (auto u : g.nodes)
        nodes.add(u);

    edges.reserve(g.edges.capacity());
    successors.resize(g.capacity());
    predecessors.resize(g.capacity());

    successors_weight.resize(g.capacity(), 0);
    predecessors_weight.resize(g.capacity(), 0);

    num_edges = 0;
    node_weight = g.node_weight;
    edge_weight = 0;

    for (auto u : g.nodes)
      for (auto n : g.neighbors[u])
        add_edge(n.vertex, u, n.weight());

#ifdef _VERIFY_MCGRAPH
    verify("after copy");
#endif
  }

  dyngraph<WEIGHT> &operator=(const dyngraph<WEIGHT> &g) {
    // copy operator doesn't work, need to make sure that the neighbors point to
    // the right graph's edges

    assert(false);
    //
    // deep_copy(g);
    // verify("after operator=");
    // return *this;
  }

  // helpers
  size_t size() const { return nodes.size(); }
  size_t capacity() const { return nodes.capacity(); }
  bool null() const { return size() == 0; }
  bool empty() const { return num_edges == 0; }
  bool full() const { return size() * size() == num_edges; }
  size_t outdegree(const int x) const { return successors[x].size(); }
  size_t indegree(const int x) const { return predecessors[x].size(); }
  float get_density() const {
    return (float)num_edges / (float)(size() * size());
  }

  void add_node(const int x, const WEIGHT w) {
    weight[x] = w;
    node_weight += w;
  }

  void add_node(const int x) {

#ifdef _VERIFY_MCGRAPH
    verify("before add node");
#endif

    assert(!nodes.contain(x));

    nodes.add(x);
    // num_edges += indegree(x);
    num_edges += outdegree(x);

    node_weight += weight[x];
    edge_weight += successors_weight[x];
    edge_weight += predecessors_weight[x];

    for (auto &y : successors[x]) {
      y.set_origin_rank(predecessors[y.vertex].size());
      predecessors[y.vertex].push_back(neighbor_info<WEIGHT>(x, y.edge));
      successors_weight[y.vertex] -= y.weight();
    }

    for (auto &y : predecessors[x]) {
      y.set_target_rank(successors[y.vertex].size());
      successors[y.vertex].push_back(neighbor_info<WEIGHT>(x, y.edge));
      predecessors_weight[y.vertex] -= y.weight();
      if (y.vertex != x)
        ++num_edges;
    }

#ifdef _VERIFY_MCGRAPH
    verify("after add node");
#endif
  }

  void rem_node(const int x) {
#ifdef _VERIFY_MCGRAPH
    verify("before rem node");
#endif

    assert(nodes.contain(x));

    nodes.remove(x);

    num_edges -= outdegree(x);
    // num_edges -= indegree(x);

    node_weight -= weight[x];
    edge_weight -= successors_weight[x];
    edge_weight -= predecessors_weight[x];

    for (auto &y : successors[x]) {
      predecessors_weight[y.vertex] -= y.weight();

      auto rx = y.get_origin_rank();

      // swap x and the last neighbor of y
      auto z{predecessors[y.vertex].back()};
      predecessors[y.vertex].pop_back();

      if (z.vertex != x) {
        predecessors[y.vertex][rx] = z;
        z.set_origin_rank(rx);
      }
    }

    for (auto &y : predecessors[x]) {
      successors_weight[y.vertex] -= y.weight();

      auto rx = y.get_target_rank();

      // swap x and the last neighbor of y
      auto z{successors[y.vertex].back()};
      successors[y.vertex].pop_back();

      if (z.vertex != x) {
        successors[y.vertex][rx] = z;
        z.set_target_rank(rx);
      }

      if (y.vertex != x)
        --num_edges;
    }

// std::cout << " -> " <<  node_weight << std::endl;

#ifdef _VERIFY_MCGRAPH
    verify("after rem node");
#endif
  }

  // edge addition/removal
  void add_edge(const int a, const int b, const WEIGHT w) {

#ifdef _VERIFY_MCGRAPH
    verify("before add edge");
#endif

    edge_info<WEIGHT> e(predecessors[b].size(), successors[a].size(), w);
		
		edges.push_back(e);
		
    neighbor_info<WEIGHT> ni_b(b, &(edges[0]) + edges.size() - 1);
    neighbor_info<WEIGHT> ni_a(a, &(edges[0]) + edges.size() - 1);

    successors[a].push_back(ni_b);
    predecessors[b].push_back(ni_a);

    successors_weight[a] += w;
    predecessors_weight[b] += w;

    edge_weight += w;

    
    num_edges += 1;

#ifdef _VERIFY_MCGRAPH
    verify("after add edge");
#endif
  }
	
  // change the weight on an edge [which is in the successor list of o]
  void set_weight(const int o, const int t, const neighbor_info<WEIGHT> n, const WEIGHT w) {

#ifdef _VERIFY_MCGRAPH
    verify("before set succ weight");
#endif
		
    auto p{n.weight()};
    n.set_weight(w);

		// std::cout << o << "." << t << ": " << p << "->" << w << std::endl;

    successors_weight[o] -= p;
    predecessors_weight[t] -= p;
    successors_weight[o] += w;
    predecessors_weight[t] += w;

    edge_weight += w;
    edge_weight -= p;

#ifdef _VERIFY_MCGRAPH
    verify("after set succ weight");
#endif
  }

//   // change the weight on an edge [which is in the successor list of o]
//   void set_successor_weight(const int o, const neighbor_info<WEIGHT> n, const WEIGHT w) {
//
// #ifdef _VERIFY_MCGRAPH
//     verify("before set succ weight");
// #endif
//
// 		auto t{n.vertex};
//     auto p{n.weight()};
//     n.set_weight(w);
//
// 		// std::cout << o << "." << t << ": " << p << "->" << w << std::endl;
//
//     successors_weight[o] -= p;
//     predecessors_weight[t] -= p;
//     successors_weight[o] += w;
//     predecessors_weight[t] += w;
//
//     edge_weight += w;
//     edge_weight -= p;
//
// #ifdef _VERIFY_MCGRAPH
//     verify("after set succ weight");
// #endif
//   }
//
//   // change the weight on an edge [which is in the predecessor list of t]
//   void set_predecessor_weight(const int t, const neighbor_info<WEIGHT> n, const WEIGHT w) {
//
// #ifdef _VERIFY_MCGRAPH
//     verify("before set succ weight");
// #endif
//
// 		auto o{n.vertex};
//     auto p{n.weight()};
//     n.set_weight(w);
//
// 		// std::cout << o << "." << t << ": " << p << "->" << w << std::endl;
//
//     successors_weight[o] -= p;
//     predecessors_weight[t] -= p;
//     successors_weight[o] += w;
//     predecessors_weight[t] += w;
//
//     edge_weight += w;
//     edge_weight -= p;
//
// #ifdef _VERIFY_MCGRAPH
//     verify("after set succ weight");
// #endif
//   }

  // printing
  std::ostream &describe(std::ostream &os, const int verbosity) const;

  // debug
  void verify(const char *msg);
};

template<class WEIGHT>
std::ostream& dyngraph<WEIGHT>::describe(std::ostream& os, const int verbosity) const
{

  switch (verbosity) {
  case 0:
    os << "n=" << nodes.size() << " m=" << num_edges << std::endl;
    break;

  case 1:
    os << "V=" << nodes << " m=" << num_edges;
    break;

  case 2:
    for (auto v : nodes) {
      os << std::setw(2) << v << " (" << weight[v] << "):";
      for (auto n : successors[v]) {
        os << " " << std::setw(2) << n.vertex;
      }
      os << std::endl;
    }
    break;
  case 3:
    for (auto v : nodes) {
      os << std::setw(2) << v << " (" << std::setw(3) << weight[v] << "):";
      for (auto u : successors[v]) {
        os << " " << std::setw(2) << u.vertex << "/" << std::setw(3)
           << std::left << u.weight() << std::right;
      }
      os << "[" << successors_weight[v] << "]";
      for (auto u : predecessors[v]) {
        os << " " << std::setw(2) << u.vertex << "/" << std::setw(3)
           << std::left << u.weight() << std::right;
      }
      os << "[" << predecessors_weight[v] << "]";
      os << std::endl;
    }
  }

  return os;
}

template<class WEIGHT>
void dyngraph<WEIGHT>::verify(const char* msg)
{

  // std::cout << *this << std::endl;

  auto count{num_edges};
  std::vector<int> f(nodes.capacity(), 0);
  std::vector<int> b(nodes.capacity(), 0);
  std::vector<int> Nghb(nodes.capacity());

  WEIGHT node_weight_count{0};
  WEIGHT edge_weight_count{0};
  std::vector<WEIGHT> fw(nodes.capacity(), 0);
  std::vector<WEIGHT> bw(nodes.capacity(), 0);

  for (auto u : nodes) {
    std::fill(Nghb.begin(), Nghb.end(), 0);

    node_weight_count += weight[u];
    for (auto n : successors[u]) {
      if (!nodes.contain(n.vertex)) {
        std::cout << msg << ": " << n.vertex << " (succ of " << u
                  << ") is not in the graph!\n";
        std::cout << *this << std::endl;
        exit(1);
      }

      if (predecessors[n.vertex][n.get_origin_rank()].vertex != u) {
        std::cout << msg << ": " << u << " is not at its rank ("
                  << n.get_origin_rank() << ") in " << n.vertex
                  << "'s successors list\n";
        std::cout << *this << std::endl;
        exit(1);
      }

      edge_weight_count += n.weight();
      fw[u] += n.weight();
      bw[n.vertex] += n.weight();

      if (++Nghb[n.vertex] > 1) {
        std::cout << msg << ": " << u << " has " << n.vertex
                  << " twice as successor!\n";
        exit(1);
      }

      ++f[u];
      ++b[n.vertex];
      --count;
    }
  }

  if (count) {
    std::cout << msg << " (1):" << num_edges << " / " << (num_edges - count)
              << std::endl;
    exit(1);
  }
  if (edge_weight != edge_weight_count) {
    std::cout << msg << " inconsistency in total edge weight (1) ("
              << edge_weight / 2 << "/" << edge_weight_count << ")\n";
    exit(1);
  }

  edge_weight_count = 0;
  count = num_edges;
  for (auto u : nodes) {
    std::fill(Nghb.begin(), Nghb.end(), 0);

    for (auto n : predecessors[u]) {
      if (!nodes.contain(n.vertex)) {
        std::cout << msg << ": " << n.vertex << " (pred of " << u
                  << ") is not in the graph!\n";
        std::cout << *this << std::endl;
        exit(1);
      }

      if (successors[n.vertex][n.get_target_rank()].vertex != u) {
        std::cout << msg << ": " << u << " is not at its rank ("
                  << n.get_target_rank() << ") in " << n.vertex
                  << "'s predecessors list\n";
        std::cout << *this << std::endl;
        exit(1);
      }

      if (++Nghb[n.vertex] > 1) {
        std::cout << msg << ": " << u << " has " << n.vertex
                  << " twice as predecessor!\n";
        exit(1);
      }

      edge_weight_count += n.weight();
      --count;
    }
  }

  if (count) {
    std::cout << msg << " (2):" << num_edges << " / " << (num_edges / 2 - count)
              << std::endl;
    exit(1);
  }
  if (edge_weight != edge_weight_count) {
    std::cout << msg << " inconsistency in total edge weight (2) ("
              << edge_weight / 2 << "/" << edge_weight_count << ")\n";
    exit(1);
  }

  if (node_weight != node_weight_count) {
    std::cout << msg << " inconsistency in total node weight (" << node_weight
              << "/" << node_weight_count << ")\n";
    exit(1);
  }

  for (auto u : nodes) {
    assert(b[u] == indegree(u));
    assert(f[u] == outdegree(u));
    // if (fw[u] != bw[u]) {
    //   std::cout << msg << " weight symmetry problem for " << u << ": " <<
    //   fw[u]
    //             << "/" << bw[u] << "\n";
    //   exit(1);
    // }
    if (fw[u] != successors_weight[u]) {
      std::cout << msg << " succ weight consistency problem for " << u << ": "
                << fw[u] << "/" << successors_weight[u] << "\n";

      for (auto n : successors[u]) {
        std::cout << " " << n.weight();
      }
      std::cout << std::endl;

      exit(1);
    }

    if (bw[u] != predecessors_weight[u]) {
      std::cout << msg << " pred weight consistency problem for " << u << ": "
                << bw[u] << "/" << predecessors_weight[u] << "\n";

      for (auto n : predecessors[u]) {
        std::cout << " " << n.weight();
      }
      std::cout << std::endl;

      exit(1);
    }

    // assert(fw[u] == bw[u] and bw[u] == neighbor_weight[u]);
  }
}

template<class WEIGHT>
std::ostream& operator<<(std::ostream& os, const dyngraph<WEIGHT>& g)
{
  return g.describe(os, 3);
}

} // namespace block

#endif
