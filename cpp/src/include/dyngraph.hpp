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
    // rank of smallest vertex in neighbor list of the highest
    rank[0] = a;
    // rank of highest vertex in neighbor list of the smallest
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
  int endpoint() const { return vertex; }

  int vertex;

  edge_info<WEIGHT> *edge;

  neighbor_info(const int v, edge_info<WEIGHT> *e) : vertex{v}, edge(e) {}
  neighbor_info(const neighbor_info<WEIGHT> &ni, edge_info<WEIGHT> *e)
      : vertex{ni.vertex}, edge(e) {}

  int get_rank(const int u) const { return edge->rank[u > vertex]; }

  // set the rank of u in vertex's neighbor list
  void set_other_rank(const int u, const int r) const {
    edge->rank[u > vertex] = r;
  }
  // set the rank of vertex in u's neighbor list
  void set_self_rank(const int u, const int r) const {
    edge->rank[u < vertex] = r;
  }
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
  // the list of neighbors of each node
  std::vector<std::vector<neighbor_info<WEIGHT>>> neighbors;

  // the weights on the nodes
  std::vector<WEIGHT> weight;
  // the weights on the edges
  std::vector<WEIGHT> neighbor_weight;

  // total node weight
  WEIGHT node_weight;
  // total edge weight
  WEIGHT edge_weight;

  bool directed;

  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::iterator begin() {
    return begin(neighbors);
  }
  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::reverse_iterator
  rbegin() {
    return rbegin(neighbors);
  }

  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::iterator end() {
    return end(neighbors);
  }
  typename std::vector<std::vector<neighbor_info<WEIGHT>>>::reverse_iterator
  rend() {
    return rend(neighbors);
  }

  std::vector<neighbor_info<WEIGHT>> &operator[](const int v) {
    return neighbors[v];
  }

  dyngraph() {}
  dyngraph(const int n, const bool d = true)
  // : weight(n,0), neighbor_weight(n,0)
  {
    initialise(n, d);
  }

  template <class graph_struct>
  dyngraph(const graph_struct &g, const bool r = false)
  // : weight(n,0), neighbor_weight(n,0)
  {
    deep_copy(g, r);
  }

  void initialise(const int n, const bool d = true) {
    directed = d;
    weight.resize(n, 0);
    neighbor_weight.resize(n, 0);

    nodes.reserve(n);
    nodes.fill();

    neighbors.resize(n);

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

  template <class graph_struct>
  void deep_copy(const graph_struct &g, const bool reverse) {
    weight = g.weight;
    neighbor_weight.resize(g.capacity());
    nodes.reserve(g.capacity());
    if (g.size() == g.capacity())
      nodes.fill();
    else
      for (auto u : g.nodes)
        nodes.add(u);

    edges.reserve(g.edges.size());
    neighbors.resize(g.capacity());

    num_edges = 0;
    node_weight = g.node_weight;
    edge_weight = 0;

    for (auto u : g.nodes)
      for (auto n : g.neighbors[u])
        if (directed) {
          if (reverse)
            add_directed_edge(n.vertex, u, n.weight());
          else
            add_directed_edge(u, n.vertex, n.weight());
        } else if (n.vertex > u) {
          add_undirected_edge(u, n.vertex, n.weight());
        }

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
  bool full() const { return size() * (size() - 1) == 2 * num_edges; }
  size_t degree(const int x) const { return neighbors[x].size(); }
  double get_density() const {
    return (double)num_edges / (double)(size() * (size() - 1) / 2);
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
    num_edges += degree(x);

    node_weight += weight[x];
    edge_weight += neighbor_weight[x];

    for (auto &y : neighbors[x]) {
      y.set_other_rank(x, neighbors[y.vertex].size());
      neighbors[y.vertex].push_back(neighbor_info<WEIGHT>(x, y.edge));
      neighbor_weight[y.vertex] += y.weight();
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
    // std::cout << "rem " << x << " (" << weight[x] << ") from " << node_weight
    // << "\n" << *this ;

    nodes.remove(x);
    num_edges -= degree(x);

    node_weight -= weight[x];
    edge_weight -= neighbor_weight[x];

    for (auto &y : neighbors[x]) {
      neighbor_weight[y.vertex] -= y.weight();

      auto rx = y.get_rank(x);

      // swap x and the last neighbor of y
      auto z{neighbors[y.vertex].back()};
      neighbors[y.vertex].pop_back();

      if (z.vertex != x) {
        neighbors[y.vertex][rx] = z;
        z.set_self_rank(y.vertex, rx);
      }
    }

// std::cout << " -> " <<  node_weight << std::endl;

#ifdef _VERIFY_MCGRAPH
    verify("after rem node");
#endif
  }

  void add_edge(const int a, const int b, const WEIGHT w) {
    if (directed)
      add_directed_edge(a, b, w);
    else if (a < b)
      add_undirected_edge(a, b, w);
    else
      add_undirected_edge(b, a, w);
  }

  // edge addition/removal
  void add_directed_edge(const int a, const int b, const WEIGHT w) {

#ifdef _VERIFY_MCGRAPH
    verify("before add edge");
#endif

    edge_info<WEIGHT> e(neighbors[a].size(), neighbors[b].size(), w);
    neighbor_info<WEIGHT> ni_b(b, &(edges[0]) + edges.size());

    neighbors[a].push_back(ni_b);

    neighbor_weight[a] += w;

    edge_weight += w;

    edges.push_back(e);
    ++num_edges;

#ifdef _VERIFY_MCGRAPH
    verify("after add edge");
#endif
  }

  // edge addition/removal
  void add_undirected_edge(const int a, const int b, const WEIGHT w) {

#ifdef _VERIFY_MCGRAPH
    verify("before add edge");
#endif

    edge_info<WEIGHT> e(neighbors[a].size(), neighbors[b].size(), w);
    neighbor_info<WEIGHT> ni_b(b, &(edges[0]) + edges.size());
    neighbor_info<WEIGHT> ni_a(a, &(edges[0]) + edges.size());

    neighbors[a].push_back(ni_b);
    neighbors[b].push_back(ni_a);

    neighbor_weight[a] += w;
    neighbor_weight[b] += w;

    edge_weight += w;

    edges.push_back(e);
    ++num_edges;

#ifdef _VERIFY_MCGRAPH
    verify("after add edge");
#endif
  }

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
      for (auto n : neighbors[v]) {
        os << " " << std::setw(2) << n.vertex;
      }
      os << std::endl;
    }
    break;
  case 3:
    for (auto v : nodes) {
      os << std::setw(2) << v << " (" << std::setw(3) << weight[v] << "):";
      for (auto u : neighbors[v]) {
        os << " " << std::setw(2) << u.vertex << "/" << std::setw(3)
           << std::left << u.weight() << std::right;
      }
      os << std::endl;
    }
  }

  return os;
}

template<class WEIGHT>
void dyngraph<WEIGHT>::verify(const char* msg)
{

  auto count{2 * num_edges};
  std::vector<int> f(nodes.capacity(), 0);
  std::vector<int> b(nodes.capacity(), 0);

  WEIGHT node_weight_count{0};
  WEIGHT edge_weight_count{0};
  std::vector<WEIGHT> fw(nodes.capacity(), 0);
  std::vector<WEIGHT> bw(nodes.capacity(), 0);

  for (auto u : nodes) {

    node_weight_count += weight[u];
    for (auto n : neighbors[u]) {

      if (!nodes.contain(n.vertex)) {
        std::cout << msg << ": " << n.vertex << " is not in the graph!\n";
        exit(1);
      }

      if (neighbors[n.vertex][n.get_rank(u)].vertex != u) {
        std::cout << msg << ": " << u << " is not at its rank ("
                  << n.get_rank(u) << ") in " << n.vertex << "'s list\n";
        std::cout << *this << std::endl;
        exit(1);
      }

      edge_weight_count += n.weight();
      fw[u] += n.weight();
      bw[n.vertex] += n.weight();

      ++f[u];
      ++b[n.vertex];
      --count;
    }
  }

  if (count)
    std::cout << count << std::endl;

  assert(count == 0);

  if (node_weight != node_weight_count) {
    std::cout << msg << " inconsistency in total node weight (" << node_weight
              << "/" << node_weight_count << ")\n";
    exit(1);
  }

  if (2 * edge_weight != edge_weight_count) {
    std::cout << msg << " inconsistency in total edge weight (" << edge_weight
              << "/" << edge_weight_count / 2 << ")\n";
    exit(1);
  }

  // assert(node_weight == node_weight_count);

  for (auto u : nodes) {
    assert(f[u] == b[u] and b[u] == degree(u));
    if (fw[u] != bw[u]) {
      std::cout << msg << " weight symmetry problem for " << u << ": " << fw[u]
                << "/" << bw[u] << "\n";
      exit(1);
    }
    if (fw[u] != neighbor_weight[u]) {
      std::cout << msg << " weight consistency problem for " << u << ": "
                << fw[u] << "/" << neighbor_weight[u] << "\n";

      for (auto n : neighbors[u]) {
        std::cout << " " << n.weight();
      }
      std::cout << std::endl;

      exit(1);
    }

    assert(fw[u] == bw[u] and bw[u] == neighbor_weight[u]);
  }
}

template<class WEIGHT>
std::ostream& operator<<(std::ostream& os, const dyngraph<WEIGHT>& g)
{
  return g.describe(os, 3);
}

} // namespace block

#endif
