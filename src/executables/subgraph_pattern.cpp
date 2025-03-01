#include <fstream>
#include <gurobi/F1plus.hpp>
#include <scip/scip_affine_F2.hpp>

#include "auxiliary/cxxopts.hpp"

#include "auxiliary/GXLGraphReader.hpp"

#include "auxiliary/options.hpp"
#include "auxiliary/io.hpp"
#include "scip/utils.hpp"
#include "gurobi/F2.hpp"
#include "gurobi/FORI.hpp"
#include "gurobi/F2plus.hpp"

int main(int argc, char **argv) {
        try {
                cxxopts::Options opts("experiment", "calc GED with given formulation");
                opts.add_options()
                        ("f, formulation", "which LP Formulation to use", cxxopts::value<std::string>())(
                        "g, graph1", "id of first graph", cxxopts::value<std::string>())(
                        "h, graph2", "id of second graph", cxxopts::value<std::string>())(
                        "t, threads", "number of threads", cxxopts::value<int>()->default_value("1"))(
                        "s, setting", "which branching priority setting to use", cxxopts::value<int>()->default_value("0"))(
                        "l, timelimit", "number of seconds the solver is allowed to run for", cxxopts::value<double>()->default_value("900"))(
                        "r, seed", "random seed", cxxopts::value<int>()->default_value("1"))(
                        "disablePresolving", "disables SCIPs presolving", cxxopts::value<bool>()->default_value("false"))(
                        "relax", "set all variables to real instead of binary", cxxopts::value<bool>()->default_value("false"))(
                        "b, branchingRule", "select branching rule, 0 = scip default, 1 = 1-hop cosine similarity", cxxopts::value<int>()->default_value("0"))(
                        "solutionFile", "absolute path to a .gedsol file for the corresponding instance. throws runtime error if file doesn't correspond to graph1id_graph2id", cxxopts::value<std::string>()->default_value(""))(
                        "d, branchingDirection", "select branching direction, 0 = scip default, 1 = x_ik->0 and y_ijkl->1, 2 = x_ik->1 and y_ijkl->1", cxxopts::value<int>()->default_value("0"));


                auto arguments = opts.parse(argc, argv);
                const std::string formulation = arguments["formulation"].as<std::string>();
                const std::string g1 = arguments["graph1"].as<std::string>();
                const std::string g2 = arguments["graph2"].as<std::string>();
                const std::string solutionFile = arguments["solutionFile"].as<std::string>();
                const int threads = arguments["threads"].as<int>();
                const int branchPrio = arguments["branchingRule"].as<int>();
                const int branchDirection = arguments["branchingDirection"].as<int>();
                const int setting = arguments["setting"].as<int>();
                const int seed = arguments["seed"].as<int>();
                const double timelimit = arguments["timelimit"].as<double>();
                const bool disablePresolving = arguments["disablePresolving"].as<bool>();
                const bool relax = arguments["relax"].as<bool>();


                std::string gr1 = "../data/SubgraphPattern/largeGraphs/" + g1 + ".edgelist";
                std::string gr2 = "../data/SubgraphPattern/Patterns/" + g2 + ".edgelist";

                std::cout << gr1 << "\n" << gr2 << std::endl;
                graph<int, int> graph1 = IO::read_edgelist(gr1, "subgraph");
                graph<int, int> graph2 = IO::read_edgelist(gr2, "subgraph");
                graph1.set_graph_id(g1);
                graph2.set_graph_id(g2);


                if (graph2.number_of_edges() > graph1.number_of_edges()) {
                        std::swap(gr1, gr2);
                        std::swap(graph1, graph2);
                }


                options opt;
                opt.dataset_name_ = "subgraph-pattern";
                opt.seed_ = seed;
                opt.setting_ = setting;
                opt.timelimit_ = timelimit;
                opt.threads_ = threads;
                opt.formulation_name_ = formulation;
                opt.solutionFilePath_ = solutionFile;
                opt.G_id_ = graph1.get_graph_id();
                opt.G_filename_ = graph1.get_graph_id();
                opt.G_num_nodes_ = std::to_string(graph1.number_of_nodes());
                opt.G_num_edges_ = std::to_string(graph1.number_of_edges());
                opt.H_filename_ = graph2.get_graph_id();
                opt.H_id_ = graph2.get_graph_id();
                opt.H_num_nodes_ = std::to_string(graph2.number_of_nodes());
                opt.H_num_edges_ = std::to_string(graph2.number_of_edges());
                opt.disablePresolving = disablePresolving;
                opt.branchPrio_ = branchPrio;
                opt.branchDirection_ = branchDirection;
                opt.relaxed_ = relax;
                opt.output_lpfile_ = opt.G_id_ + "_" + opt.H_id_ + "_" + opt.formulation_name_ +
                                                                  + ".mps";

                std::cout << "G=" << g1 << ", n=" << opt.G_num_nodes_ << ", m=" << opt.G_num_edges_ << std::endl;
                std::cout << "H=" << g2 << ", n=" << opt.H_num_nodes_ << ", m=" << opt.H_num_edges_ << std::endl;


                utils<std::string, int> utility;
                std::string file;
                if (formulation == "F2+") {
                        F2plus<int, int> ilp(opt);
                        IO::setFilename("subgraph", "F2+", opt);
                        ilp.ged(graph1, graph2);
                        utility.writeJsonToFile(opt);
                } else if(formulation == "FORI"){
                        FORI<int, int> ilp(opt);
                        IO::setFilename("subgraph", "FORI", opt);
                        ilp.ged(graph1, graph2);
                        utility.writeJsonToFile(opt);
                } else if(formulation == "F1+"){
                        F1plus<int, int> ilp(opt);
                        IO::setFilename("subgraph", "F1+", opt);
                        ilp.ged(graph1, graph2);
                        utility.writeJsonToFile(opt);
                }else if(formulation == "F2") {
                        F2<int, int> ilp(opt);
                        IO::setFilename("subgraph", "F2", opt);
                        ilp.ged(graph1, graph2);
                        utility.writeJsonToFile(opt);
                } else {
                        std::cout << "Unknown Formulation name: " << formulation << std::endl;
                        return 1;
                }
        }
        catch (std::exception &e) {
                std::cout << "exception " << e.what() << std::endl;
                return 1;
        }


}

