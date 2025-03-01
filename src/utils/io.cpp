#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "auxiliary/json.hpp"
//#include "gurobi_c++.h"

#include "auxiliary/io.hpp"
#include "auxiliary/graph.hpp"
#include "auxiliary/options.hpp"


namespace IO {

    void setFilename(const std::string &datasetName, const std::string &formulation, options &opt) {
        opt.log_fname_ =
            datasetName + "_" + formulation + "_" + opt.G_id_ + "_" + opt.H_id_ + "_" + std::to_string(opt.threads_) +
            "_" + std::to_string(opt.seed_);
        if (opt.setting_ > 0) {
            opt.log_fname_ = opt.log_fname_ + "_setting" + std::to_string(opt.setting_);
        }
        if (opt.addC2F1_) {
            opt.log_fname_ = opt.log_fname_ + "_C2F1";
        } else if (opt.addC2F2_) {
            opt.log_fname_ = opt.log_fname_ + "_C2F2";
        } else if (opt.addC3_) {
            opt.log_fname_ = opt.log_fname_ + "_C3";
        } else if (opt.addC3onlySub_) {
            opt.log_fname_ = opt.log_fname_ + "_C3onlySub";
        }
        if (opt.disablePresolving) {
            opt.log_fname_ = opt.log_fname_ + "_disablePresolve";
        }
        if (opt.branchPrio_) {
            opt.log_fname_ = opt.log_fname_ + "_branch" + std::to_string(opt.branchPrio_);
        }
        if (opt.branchDirection_) {
            opt.log_fname_ = opt.log_fname_ + "_dir" + std::to_string(opt.branchDirection_);
        }
        if (opt.includeObjCoefIntoFracDegreeBranch) {
            opt.log_fname_ = opt.log_fname_ + "_InclObjCoef";
        }
        if (!opt.solutionFilePath_.empty()) {
            opt.log_fname_ = opt.log_fname_ + "_optsol";
        }
        opt.output_fname_ = opt.log_fname_ + ".json";
    }

    std::string create_output(options &opt) {
        nlohmann::json j;
        j["formulation"] = opt.formulation_name_;
        j["constant"] = opt.constant_;
        j["setting"] = opt.setting_;
        j["seed"] = opt.seed_;

        j["G_filename"] = opt.G_filename_;
        j["G_id"] = opt.G_id_;
        j["G_Nodes"] = opt.G_num_nodes_;
        j["G_Edges"] = opt.G_num_edges_;
        j["H_filename"] = opt.H_filename_;
        j["H_id"] = opt.H_id_;
        j["H_Nodes"] = opt.H_num_nodes_;
        j["H_Edges"] = opt.H_num_edges_;

        j["time"] = opt.time_;
        j["bbNodeCount"] = opt.bbnodecount_;
        j["solved"] = opt.status_ == 2;

        j["vars-full"] = opt.n_vars_full_;
        j["constraints-full"] = opt.n_cons_full_;
        j["non-Zeros"] = opt.numNZ_;
        j["simplex-iterations"] = opt.iterCount_;

        j["lp-relaxation"] = opt.lprelval;
        j["lp-rel-time"] = opt.lpreltime;
        j["MIP-lp-relaxation"] = opt.ILP_lprelval;
        j["MIP-lp-rel-time"] = opt.ILP_lpreltime;
        j["finalPrimalBound"] = opt.objval_;
        j["finalDualBound"] = opt.final_dualbound_;
        j["rootPrimalBound"] = opt.root_primal_;
        j["rootDualBound"] = opt.root_dual_;
        j["rootGap"] = opt.root_gap_;
        if (opt.status_ == 2) {
            j["lpGap"] = (opt.objval_ - opt.lprelval) / opt.objval_;
        } else {
            j["lpGap"] = -1;
        }

        j["rootTime"] = opt.root_time_;
        j["finalGap"] = opt.mipgap_;

        j["C2F1"] = opt.addC2F1_;
        j["C2F2"] = opt.addC2F2_;
        j["C3"] = opt.addC3_;
        j["C3onlySub"] = opt.addC3onlySub_;

        j["branchPrio"] = opt.branchPrio_;
        j["branchDir"] = opt.branchDirection_;
        j["InclObjCoefInBranch"] = opt.includeObjCoefIntoFracDegreeBranch;

        return j.dump(4);
    }

    std::string createGurobiOutput(options &opt) {
        nlohmann::json j;
        j["formulation"] = opt.formulation_name_;
        j["constant"] = opt.constant_;
        j["setting"] = opt.setting_;
        j["seed"] = opt.seed_;

        j["G_filename"] = opt.G_filename_;
        j["G_id"] = opt.G_id_;
        j["G_Nodes"] = opt.G_num_nodes_;
        j["G_Edges"] = opt.G_num_edges_;
        j["H_filename"] = opt.H_filename_;
        j["H_id"] = opt.H_id_;
        j["H_Nodes"] = opt.H_num_nodes_;
        j["H_Edges"] = opt.H_num_edges_;

        j["time"] = opt.time_;
        j["objval"] = opt.objval_;
        j["mip-gap"] = opt.mipgap_;
        j["bbNodeCount"] = opt.bbnodecount_;
        j["solved"] = opt.mipgap_ < 0.0004;


        j["lp-relaxation-value"] = opt.lprelval;
        j["lp-relaxation-time"] = opt.lpreltime;
        j["C2F1"] = opt.addC2F1_;
        j["C2F2"] = opt.addC2F2_;
        j["C3"] = opt.addC3_;
        j["C3onlySub"] = opt.addC3onlySub_;

        j["branchPrio"] = opt.branchPrio_;
        j["branchDir"] = opt.branchDirection_;
        j["InclObjCoefInBranch"] = opt.includeObjCoefIntoFracDegreeBranch;

        return j.dump(4);
    }

    graph<int, int> read_edgelist(const std::string &filename, const std::string &dataset) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("In read_edgelist: could not open file " + filename);
        }

        graph<int, int> edgelist_graph;
        edgelist_graph.set_dataset(dataset);
        int cur_idx = 0;


        std::string node1, node2;
        std::string line;
        while (file >> node1 >> node2) {
            if (!edgelist_graph.find_file_label(node1)) {
                edgelist_graph.add_node(cur_idx, 0);
                edgelist_graph.set_file_to_zero_indexed(node1, cur_idx);
                cur_idx++;
            }
            if (!edgelist_graph.find_file_label(node2)) {
                edgelist_graph.add_node(cur_idx, 0);
                edgelist_graph.set_file_to_zero_indexed(node2, cur_idx);
                cur_idx++;
            }

            edgelist_graph.add_edge(edgelist_graph.get_zero_indexed_from_file(node1),
                                    edgelist_graph.get_zero_indexed_from_file(node2),
                                    0);
        }

        file.close();
        return edgelist_graph;
    }

    graph<int, int> read_csv_imdb(const std::string &filename, const std::string &dataset) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("In read_csv: could not open file " + filename);
        }

        graph<int, int> csv_graph;
        csv_graph.set_dataset(dataset);
        int cur_idx = 0;

        std::string line;

        while (std::getline(file, line)) {
            std::vector<std::string> row;
            std::stringstream ss(line);
            std::string value;

            while (std::getline(ss, value, ',')) {
                row.push_back(value);
            }

            if (!csv_graph.find_file_label(row[0])) {
                csv_graph.add_node(cur_idx, 0);
                csv_graph.set_file_to_zero_indexed(row[0], cur_idx);
                cur_idx += 1;
            }
            if (!csv_graph.find_file_label(row[1])) {
                csv_graph.add_node(cur_idx, 0);
                csv_graph.set_file_to_zero_indexed(row[1], cur_idx);
                cur_idx += 1;
            }

            csv_graph.add_edge(csv_graph.get_zero_indexed_from_file(row[0]), csv_graph.get_zero_indexed_from_file(row[1]), 0);
        }

        file.close();
        return csv_graph;
    }


} // namespace IO