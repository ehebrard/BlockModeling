

#include <random>
#include <sys/resource.h>
#include <unistd.h>

// #include "intstack.hpp"
#include "global.hpp"
#include "reader.hpp"
#include "options.hpp"
#include "algorithm.hpp"


using namespace std;
using namespace block;

std::mt19937 random_generator;

double cpuTime(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// typedef options<WEIGHT> params;

int main(int argc, char *argv[]) {

  auto options = parse(argc, argv);

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

  intstack q(g.size());
  q.add(0);

  BFS(g, q);

  std::cout << q << std::endl;

  block_model m(g);
  m.compress(options);
	
  if (options.printsolution)
    std::cout << m.model << std::endl;

  // g.rem_node(13);
  // g.rem_node(3);
  //
  // std::cout << g << std::endl;

  // compress(g);
}
