

#include <random>
// #include <sys/resource.h>
// #include <unistd.h>

// #include "intstack.hpp"
#include "algorithm.hpp"
#include "global.hpp"
#include "timing.hpp"
#include "options.hpp"
#include "reader.hpp"

using namespace std;
using namespace block;

std::mt19937 random_generator;

// double cpuTime(void) {
//   struct rusage ru;
//   getrusage(RUSAGE_SELF, &ru);
//   return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
// }

template <class BM>
void compress(BM &m, std::vector<std::vector<int>> &blocks) {
  m.compress(random_generator);
  // std::vector<std::vector<int>> blocks;
  m.get_blocks(blocks);
}

int main(int argc, char *argv[]) {

  auto options = parse(argc, argv);

  if (options.verbosity >= options::QUIET)
    options.describe(std::cout);

  random_generator.seed(options.seed);

  dgraph g;

  double start_time{cpuTime()};

  dimacs::read_graph(options.instance_file,
                     [&](int nv, int ne) {
                       g.initialise(nv);
                       g.edges.reserve(ne);
                     },
                     [&](int u, wtype w) { g.add_node(u - 1, w); },
                     [&](int u, int v, wtype w = 1) {
                       if (u != v)
                         g.add_edge(u - 1, v - 1, w);
                     });

  double read_time{cpuTime()};

  // if (options.printgraph)
    std::cout << g << std::endl;

  std::vector<std::vector<int>> blocks;

  // block_model<dyngraph<wtype>>& m{NULL}; //(g, blocks);

  if (options.size > 0) {
    intstack nodes(g.capacity());
    nodes.fill();
    auto k{options.size};
    auto bi{0};
    blocks.resize(k);
    while (!nodes.empty()) {
      auto v{nodes[random_generator() % nodes.size()]};
      blocks[bi++ % k].push_back(v);
      nodes.remove(v);
    }
  } else {
    blocks.resize(g.capacity());
    for (auto v{0}; v < g.capacity(); ++v)
      blocks[v].push_back(v);
  }

  block_model<dyngraph<wtype>> m(g, blocks, options);
  // m = warm;
  // block_model m(g, blocks);
  compress(m, blocks);

  if (options.printgraph) {

    auto k = m.model.size();
    std::vector<std::vector<float>> densities(k);
    for (auto d{begin(densities)}; d != end(densities); ++d)
      d->resize(k, 0);

    for (auto u : m.model.nodes) {
      auto i{m.model.nodes.index(u)};
      for (auto e : m.model.successors[u]) {
        auto j{m.model.nodes.index(e.endpoint())};
        densities[i][j] =
            static_cast<float>(e.weight()) /
            static_cast<float>(blocks[i].size() * blocks[j].size());
      }
    }

    for (auto i{0}; i < k; ++i) {
      std::cout << std::setw(4) << blocks[i].size() << ":";
      for (auto j{0}; j < k; ++j) {
        std::cout << " " << std::setw(7) << std::setprecision(3)
                  << densities[i][j];
      }
      std::cout << std::endl;
    }
  }
  // std::cout << m.model << std::endl;

  // } else {
  //   block_model<dyngraph<wtype>> m(g, options);
  //   // m = cold;
  //   compress(m, blocks);
  // }

  // compress(m, blocks, options);

  std::ostream &os(options.outfile());
  for (auto &b : blocks) {
    for (auto x : b)
      os << " " << x;
    os << std::endl;
    }
}
