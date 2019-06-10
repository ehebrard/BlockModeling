

#include <random>
#include <sys/resource.h>
#include <unistd.h>

// #include "intstack.hpp"
#include "algorithm.hpp"
#include "global.hpp"
#include "options.hpp"
#include "reader.hpp"

using namespace std;
using namespace block;

std::mt19937 random_generator;

double cpuTime(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
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

  if (options.printgraph)
    std::cout << g << std::endl;

  block_model m(g);
  m.compress(options, random_generator);

  if (options.printsolution) {
    std::vector<std::vector<int>> blocks;
    m.get_blocks(blocks);
    for (auto &b : blocks) {
      for (auto x : b)
        std::cout << " " << x;
      std::cout << std::endl;
    }
  }
}
