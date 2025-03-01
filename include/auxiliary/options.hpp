#ifndef GEDC_MY_OPTIONS_HPP
#define GEDC_MY_OPTIONS_HPP

#include <string>
#include <utility>

class options {

public:
        std::string solutionFilePath_;
        std::string dataset_name_;
        std::string formulation_name_;
        std::string output_fname_;
        std::string output_lpfile_;
        std::string log_fname_;
        std::string output_;
        /// for PROTEIN the filename differs from graph id
        std::string G_filename_;
        std::string G_id_;
        std::string G_num_nodes_;
        std::string G_num_edges_;
        /// for PROTEIN the filename differs from graph id
        std::string H_filename_;
        std::string H_id_;
        std::string H_num_nodes_;
        std::string H_num_edges_;
        bool relaxed_ = false;
        int threads_ = 1;
        /// setting for different branching priorities
        int setting_ = 0;               // soon deprecated branch priority setting only for pop formulations,
        int branchPrio_ = 0;            // branch priority setting
        int branchDirection_ = 0;       // which direction to branch (to 0 or to 1, look at the executables/setBranchPriority() for which direction is which)
        bool includeObjCoefIntoFracDegreeBranch = false;  // should the obj function coefficient be included as tiebreak when two variables have same fractional degree
        int seed_ = 0;
        int n_vars_full_ = 0;
        int n_cons_full_ = 0;
        int numNZ_ = 0;
        int status_ = 0;
        double timelimit_ = 0;
        double lpbound_ = 0.0;
        double root_dual_ = -1.0;
        double root_primal_ = -1.0;
        double root_gap_ = -1.0;
        double root_time_ = -1.0;
        double objval_ = 0.0;
        double bbnodecount_ = 0.0;
        double iterCount_ = 0.0;
        double constant_ = 0.0;
        double mipgap_ = 0.0;
        double time_ = 0.0;
        double final_dualbound_ = 0.0;
        /// if class 2 (F1) inequalities should be added
        bool addC2F1_ = false;
        /// if class 2 (F2) inequalities should be added
        bool addC2F2_ = false;
        bool addC3_ = false;
        bool addC3onlySub_ = false;
        bool disablePresolving = false;
        bool disableSymmetry = false;
        bool Reoptimize_ = false;
        double lprelval = 0.0;
        double lpreltime = 0.0;
        double ILP_lprelval = -1.0;
        double ILP_lpreltime = -1.0;


        options() = default;

        /// set all options that are available before gurobi is built
        options(std::string formulation_name, std::string output_fname, std::string G_filename, std::string G_id,
                std::string G_num_nodes, std::string G_num_edges, std::string H_filename, std::string H_id,
                std::string H_num_nodes, std::string H_num_edges, int setting, int seed)
                : formulation_name_(std::move(formulation_name)), output_fname_(std::move(output_fname)),
                  G_filename_(std::move(G_filename)), G_id_(std::move(G_id)), G_num_nodes_(std::move(G_num_nodes)), G_num_edges_(std::move(G_num_edges)),
                  H_filename_(std::move(H_filename)), H_id_(std::move(H_id)), H_num_nodes_(std::move(H_num_nodes)), H_num_edges_(std::move(H_num_edges)),
                  setting_(setting), seed_(seed) {}


};

#endif //GEDC_MY_OPTIONS_HPP
