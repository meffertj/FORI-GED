#define GXL_GEDLIB_SHARED

#include "src/env/ged_env.hpp"

#include <fstream>

#include "auxiliary/cxxopts.hpp"
#include <auxiliary/gedlib_costs.hpp>
#include "auxiliary/GXLGraphReader.hpp"
#include "auxiliary/options.hpp"
#include "gurobi/F2.hpp"
#include "gurobi/F1plus.hpp"
#include "gurobi/F2plus.hpp"
#include "gurobi/FORI.hpp"


int main(int argc, char **argv) {
    try {
        cxxopts::Options opts("experiment", "calc GED with given formulation, different branching priorities possible");
        opts.add_options()
            ("f, formulation", "which LP Formulation to use", cxxopts::value<std::string>())(
            "g, graph1", "id of first graph", cxxopts::value<std::string>())(
            "h, graph2", "id of second graph", cxxopts::value<std::string>())(
            "t, threads", "number of threads", cxxopts::value<int>()->default_value("1"))(
            "s, setting", "which branching priority setting to use", cxxopts::value<int>()->default_value("0"))(
            "l, timelimit", "number of seconds the solver is allowed to run for", cxxopts::value<double>()->default_value("900"))(
            "r, seed", "random seed", cxxopts::value<int>()->default_value("1"))(
            "C2F1", "add class 2 inequalities for F1", cxxopts::value<bool>()->default_value("false"))(
            "C2F2", "add class 2 inequalities for F2", cxxopts::value<bool>()->default_value("false"))(
            "C3", "add class 3 inequalities", cxxopts::value<bool>()->default_value("false"))(
            "C3onlySub", "adds C3 inequalities but only for edge substitution", cxxopts::value<bool>()->default_value("false"))(
            "disablePresolving", "disables SCIPs presolving", cxxopts::value<bool>()->default_value("false"))(
            "relax", "set all variables to real instead of binary", cxxopts::value<bool>()->default_value("false"))(
            "b, branchingRule", "select branching rule, 0 = scip default, 1 = 1-hop cosine similarity", cxxopts::value<int>()->default_value("0"))(
            "solutionFile", "absolute path to a .gedsol file for the corresponding instance. throws runtime error if file doesn't correspond to graph1id_graph2id", cxxopts::value<std::string>()->default_value(""))(
            "d, branchingDirection", "select branching direction, 0 = scip default, 1 = x_ik->0 and y_ijkl->1, 2 = x_ik->1 and y_ijkl->1", cxxopts::value<int>()->default_value("0"))(
            "w, writeInFolder", "path in which to write solution file", cxxopts::value<std::string>()->default_value(""));


        auto arguments = opts.parse(argc, argv);
        const std::string formulation = arguments["formulation"].as<std::string>();
        const std::string g1 = arguments["graph1"].as<std::string>();
        const std::string g2 = arguments["graph2"].as<std::string>();
        const std::string solutionFile = arguments["solutionFile"].as<std::string>();
        const std::string writePath = arguments["writeInFolder"].as<std::string>();
        const int threads = arguments["threads"].as<int>();
        const int branchPrio = arguments["branchingRule"].as<int>();
        const int branchDirection = arguments["branchingDirection"].as<int>();
        const int setting = arguments["setting"].as<int>();
        const int seed = arguments["seed"].as<int>();
        const double timelimit = arguments["timelimit"].as<double>();
        const bool disablePresolving = arguments["disablePresolving"].as<bool>();
        const bool C2F1 = arguments["C2F1"].as<bool>();
        const bool C2F2 = arguments["C2F2"].as<bool>();
        const bool C3 = arguments["C3"].as<bool>();
        const bool C3onlySub = arguments["C3onlySub"].as<bool>();
        const bool relax = arguments["relax"].as<bool>();


        if (C2F1 and C2F2) {
            throw std::runtime_error("C2F1 and C2F2 must not be selected simultaneously!");
        }


        std::string gr1 = "../data/AIDS/" + g1 + ".gxl";
        std::string gr2 = "../data/AIDS/" + g2 + ".gxl";

        std::cout << gr1 << "\n" << gr2 << std::endl;
        graph<std::string, int> graph1 = GXLGraphReader::read_AIDS(gr1);
        graph<std::string, int> graph2 = GXLGraphReader::read_AIDS(gr2);

        if (graph1.number_of_edges() < graph2.number_of_edges()) {
            std::swap(graph1, graph2);
        }

        /// edit costs are calculated using GEDLIB
        getGEDLIBcosts<std::string, int> getAIDSEditCosts(&graph1, &graph2);
        getAIDSEditCosts.getEditCosts();

        options opt;
        opt.dataset_name_ = "aids";
        opt.seed_ = seed;
        opt.setting_ = setting;
        opt.timelimit_ = timelimit;
        opt.threads_ = threads;
        opt.formulation_name_ = formulation;
        opt.solutionFilePath_ = solutionFile;
        opt.G_id_ = g1;
        opt.G_filename_ = graph1.get_graph_id();
        opt.G_num_nodes_ = std::to_string(graph1.number_of_nodes());
        opt.G_num_edges_ = std::to_string(graph1.number_of_edges());
        opt.H_filename_ = graph2.get_graph_id();
        opt.H_id_ = g2;
        opt.H_num_nodes_ = std::to_string(graph2.number_of_nodes());
        opt.H_num_edges_ = std::to_string(graph2.number_of_edges());
        opt.addC2F1_ = C2F1;
        opt.addC2F2_ = C2F2;
        opt.addC3_ = C3;
        opt.addC3onlySub_ = C3onlySub;
        opt.disablePresolving = disablePresolving;
        opt.branchPrio_ = branchPrio;
        opt.branchDirection_ = branchDirection;
        opt.relaxed_ = relax;

        opt.output_fname_ = writePath + "/" + opt.G_id_ + "_" + opt.H_id_ + "_" + opt.formulation_name_ +
                            +"_" + std::to_string(opt.seed_) + ".json";
        opt.output_lpfile_ = writePath + "/" + opt.G_id_ + "_" + opt.H_id_ + "_" + opt.formulation_name_ +
                             +".mps";

        std::cout << "G=" << g1 << ", n=" << opt.G_num_nodes_ << ", m=" << opt.G_num_edges_ << std::endl;
        std::cout << "H=" << g2 << ", n=" << opt.H_num_nodes_ << ", m=" << opt.H_num_edges_ << std::endl;


        std::string file;
        if (formulation == "F2+") {
            F2plus<std::string, int> ilp(opt);
            ilp.ged(graph1, graph2);
            opt.output_ = IO::createGurobiOutput(opt);
            IO::writeJsonToFile(opt);
        } else if (formulation == "FORI") {
            FORI<std::string, int> ilp(opt);
            ilp.ged(graph1, graph2);
            opt.output_ = IO::createGurobiOutput(opt);
            IO::writeJsonToFile(opt);
        } else if (formulation == "F1+") {
            F1plus<std::string, int> ilp(opt);
            ilp.ged(graph1, graph2);
            opt.output_ = IO::createGurobiOutput(opt);
            IO::writeJsonToFile(opt);
        } else if (formulation == "F2") {
            F2<std::string, int> ilp(opt);
            ilp.ged(graph1, graph2);
            opt.output_ = IO::createGurobiOutput(opt);
            IO::writeJsonToFile(opt);
        } else {
            std::cout << "Unknown Formulation name: " << formulation << std::endl;
            return 1;
        }
    }
    catch (std::exception &e) {
        std::cout << "exception " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
