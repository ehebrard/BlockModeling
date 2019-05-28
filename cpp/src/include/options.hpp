#ifndef __BLOCK__OPTIONS_HPP
#define __BLOCK__OPTIONS_HPP


#include <tclap/CmdLine.h>

#include <string>
#include <iostream>
#include <numeric>

namespace block
{
	
struct options {
    // outputs a nice description of all options
    void describe(std::ostream&);

    // the actual options
    std::string cmdline; // for reference
    std::string instance_file;
		
    enum verbosity { SILENT = 0, QUIET, NORMAL, YACKING };
    int verbosity;
				
		bool printgraph;
		bool printsolution;
                bool checked;

                int seed;
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
                                          "instance file", true, "", "string");

  cmd.add<ValueArg<int>>(
      opt.verbosity, "", "verbosity",
      "verbosity level (0:silent,1:quiet,2:improvements only,3:verbose", false,
      2, "int");

  cmd.add<ValueArg<int>>(opt.seed, "", "seed", "random seed", false, 12345,
                         "int");

  cmd.add<SwitchArg>(opt.printgraph, "", "printgraph", "display the graph",
                     false);

  cmd.add<SwitchArg>(opt.printsolution, "", "printsolution",
                     "display the solution", false);

  cmd.add<SwitchArg>(opt.checked, "", "checked",
                     "check the real objective at each move", false);

  cmd.parse(argc, argv);
  return opt;
}

void options::describe(std::ostream& os)
{
  os << "[options] block modelling\n";
  os << "[options] cmdline = " << cmdline << "\n";
  os << "[options] instance file = " << instance_file << "\n";
  os << "[options] verbosity = " << verbosity << "\n";
  os << "[options] random seed = " << seed << "\n";
  os << "[options] checked = " << (checked ? "yes" : "no") << "\n";
  os << std::endl;
}

} // namespace block

#endif
