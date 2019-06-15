#ifndef __BLOCK__OPTIONS_HPP
#define __BLOCK__OPTIONS_HPP


#include <tclap/CmdLine.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

namespace block
{

struct options {
  // outputs a nice description of all options
  void describe(std::ostream &);

  // the actual options
  std::string cmdline; // for reference
  std::string instance_file;

  enum verbosity { SILENT = 0, QUIET, NORMAL, YACKING };
  int verbosity;

  enum policy { BEST = 0, FIRST };
  int policy;

  bool printgraph;
  bool stable;

  int seed;
  int size;

  float check_epsilon;
  float epsilon;
  // float alpha;
  // float beta;

  std::string output;

  std::ofstream *os;
  std::ostream &outfile() {

    if (output == "stdout")
      return std::cout;
    os = new std::ofstream(output.c_str(), std::ios_base::out);

    return *os;
  }

  bool randomized() { return seed > 0; }
  bool checked() { return check_epsilon > 0; }
};

struct argbase {
  virtual ~argbase() {}
  virtual void assign() = 0;
};

template <typename Opt, typename ClapArg, typename E = void>
struct arg : public argbase {
  ClapArg carg;
  Opt &opt;

  template <typename... T>
  arg(TCLAP::CmdLine &cmd, Opt &opt, T &&... args)
      : carg(std::forward<T>(args)...), opt(opt) {
    cmd.add(carg);
  }

  virtual void assign() override { opt = carg.getValue(); }
};

template <typename Opt, typename ClapArg>
struct arg<Opt, ClapArg, typename std::enable_if<std::is_enum<Opt>{}>::type>
    : public argbase {
  ClapArg carg;
  Opt &opt;

  template <typename... T>
  arg(TCLAP::CmdLine &cmd, Opt &opt, T &&... args)
      : carg(std::forward<T>(args)...), opt(opt) {
    cmd.add(carg);
  }

  virtual void assign() override {
    opt =
        static_cast<typename std::remove_reference<Opt>::type>(carg.getValue());
  }
};

struct cmdline {
  TCLAP::CmdLine cmd;
  std::vector<std::unique_ptr<argbase>> args;

  cmdline(const std::string &message, const char delimiter = ' ',
          const std::string &version = "none", bool helpAndVersion = true)
      : cmd(message, delimiter, version, helpAndVersion) {}

  template <typename ClapArg, typename Opt, typename... T>
  void add(Opt &opt, T &&... clapargs) {
    args.emplace_back(std::move(std::make_unique<arg<Opt, ClapArg>>(
        cmd, opt, std::forward<T>(clapargs)...)));
  }

  void parse(int argc, char *argv[]) {
    cmd.parse(argc, argv);
    for (auto &arg : args)
      arg->assign();
  }
};

options parse(int argc, char* argv[])
{
  using namespace TCLAP;
  using namespace std::string_literals;
  cmdline cmd("block modelling", ' ');

  options opt;
  opt.cmdline = std::accumulate(
      argv, argv + argc, ""s,
      [&](std::string acc, const char *arg) { return acc + " " + arg; });

  cmd.add<UnlabeledValueArg<std::string>>(opt.instance_file, "file",
                                          "Instance file in dimacs format", true, "", "string");

  cmd.add<ValueArg<int>>(
      opt.verbosity, "", "verbosity",
      "Verbosity level (0:silent,1:quiet,2:normal,3:verbose", false,
      2, "int");

  cmd.add<ValueArg<int>>(opt.policy, "", "policy",
                         "Move selection policy (0:best,1:first improving", false, 0, "int");

  cmd.add<ValueArg<int>>(opt.seed, "", "seed", "Set the random seed (default: not randomized)", false, 0, "int");

  cmd.add<ValueArg<int>>(opt.size, "", "size",
                         "Maximum number of blocks (default: no limit)", false,
                         0, "int");

  cmd.add<SwitchArg>(opt.printgraph, "", "printgraph", "Display the graph",
                     false);

  cmd.add<ValueArg<std::string>>(opt.output, "o", "output",
                                 "Output file (default stdout)", false,
                                 "stdout", "string");

  cmd.add<ValueArg<float>>(opt.check_epsilon, "", "checked",
                           "Check w.r.t. brute force (argument = acceptable error)", false,
                           0.0, "float");

  cmd.add<SwitchArg>(opt.stable, "", "stable",
                     "Do not move the representative of a non-singleton bag",
                     false);

  cmd.add<ValueArg<float>>(opt.epsilon, "", "epsilon", "Minimum gain to accept a move", false,
                           1e-3, "float");

  // cmd.add<ValueArg<float>>(opt.alpha, "", "alpha", "model weight", false,
  // 1.0,
  //                          "float");
  //
  // cmd.add<ValueArg<float>>(opt.beta, "", "beta", "error weight", false, 1.0,
  //                          "float");

  cmd.parse(argc, argv);
  return opt;
}

void options::describe(std::ostream& os)
{
  os << "[options] block modelling\n";
  os << "[options] cmdline = " << cmdline << "\n";
  os << "[options] instance file = " << instance_file << "\n";
  os << "[options] verbosity = "
     << (verbosity == SILENT
             ? "silent"
             : (verbosity == QUIET
                    ? "quiet"
                    : (verbosity == NORMAL ? "normal" : "yacking")))
     << "\n";
  os << "[options] policy = "
     << (policy == policy::FIRST ? "first improving move" : "best move")
     << (stable ? " (stable bags)" : "") << "\n";
  os << "[options] minimum gain = " << epsilon << "\n";
  os << "[options] ";
  if (randomized())
    os << "randomized, seed = " << seed << "\n";
  else
    os << "not randomized\n";
  os << "[options] checked = ";
  if (checked())
    os << "yes (" << check_epsilon << ")\n";
  else
    os << "no \n";
  os << std::endl;
}

} // namespace block

#endif
