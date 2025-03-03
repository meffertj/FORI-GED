#ifndef GEDC_F2_HPP
#define GEDC_F2_HPP

#include "gurobi_c++.h"
#include "auxiliary/graph.hpp"
#include "utils.hpp"

#include "auxiliary/io.hpp"
#include "root_relaxation_cb.hpp"
#include "no_branching_cb.hpp"
#include "auxiliary/options.hpp"


template<typename T, typename U>
class F2 {
private:

        bool testing_ = false;
        double constant_ = 0.0;
        bool relax_ = false;
        char VAR_TYPE_ = GRB_BINARY;
        int threads_ = 1;
        int seed_ = 1;

public:

        F2() = default;

        explicit F2(options &opt) : opt_(opt) {}


        inline void relax() {
                relax_ = true;
                opt_.relaxed_ = true;
        }

        inline void set_testing() { testing_ = true; }

        inline void set_threads(int threads) { threads_ = threads; }

        inline void set_seed(int seed) { seed_ = seed; }

        std::string log_name_;
        std::string output_fname_;
        double timelimit_ = 0;
        options &opt_;

        // ------------------- template function implementations -------------------------------------------------------


        std::string ged(graph<T, U> &G, graph<T, U> &H) {

                auto n_g = G.number_of_nodes(); /** Knotenmenge von G */
                auto n_h = H.number_of_nodes(); /** Knotenmenge von H */
                auto m_g = G.number_of_edges(); /** Kantenmenge von G */
                auto m_h = H.number_of_edges(); /** Kantenmenge von H */

                try {
                        GRBEnv *env;
                        GRBVar **node_sub = new GRBVar*[n_g];
                        for (int i = 0; i < n_g; i++) {
                                node_sub[i] = new GRBVar[n_h];
                        }
                        GRBVar **edge_sub = new GRBVar *[m_g];
                        for(int i = 0; i < m_g; i++) {
                                edge_sub[i] = new GRBVar[m_h];
                        }

                        env = new GRBEnv();
                        GRBModel *model = new GRBModel(*env);
                        model->set(GRB_StringAttr_ModelName, "F2_"+opt_.dataset_name_ + "_"+opt_.G_id_ + "_" +opt_.H_id_);
                        model->set(GRB_StringParam_LogFile, opt_.log_fname_ + ".log");
                        model->set(GRB_IntParam_Threads, opt_.threads_);

                        model->set(GRB_DoubleParam_TimeLimit, opt_.timelimit_);


                        model->set(GRB_IntParam_Seed, opt_.seed_);

                        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);


                        int kl = 0;

                        std::vector<std::vector<double>> c_ik(n_g, std::vector<double>(n_h, 0));
                        std::vector<double> c_ie(n_g, 0);
                        std::vector<double> c_ek(n_h, 0);
                        std::vector<std::vector<double>> c_ijkl(m_g, std::vector<double>(m_h, 0));
                        std::vector<double> c_ije(m_g, 0);
                        std::vector<double> c_ekl(m_h, 0);


                        // cost-function
                        G.cost_function(G, H, c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);

                        // handle weird infinity cost
                        if (G.get_dataset() == "CMU") {
                                for (int i = 0; i < n_g; i++) {
                                        c_ie[i] = 0;
                                }
                                for (int k = 0; k < n_h; k++) {
                                        c_ek[k] = 0;
                                }
                        }

                        // für jeden knoten in n_g gibt es n_h viele mögliche mappings
                        for (int i = 0; i < n_g; i++) {
                                for (int k = 0; k < n_h; k++)
                                        /** coefficient of substitution variable combines cost for deletion and insertion to implicitly model deletion and insertion */
                                        node_sub[i][k] = model->addVar(0, 1, (c_ik[i][k] - c_ie[i] - c_ek[k]), GRB_BINARY,
                                                                       "map_node_" + std::to_string(i) + "_to_" + std::to_string(k));
                        }


                        // für jede kante in m_g gibt es m_h viele mögliche mappings
                        for (int ij = 0; ij < m_g; ij++) {
                                for (kl = 0; kl < m_h; kl++)
                                        /** coefficient of substitution variable combines cost for deletion and insertion to implicitly model deletion and insertion */
                                        edge_sub[ij][kl] = model->addVar(0, 1, (c_ijkl[ij][kl] - c_ije[ij] - c_ekl[kl]),
                                                                         GRB_BINARY, "map_edge_" + std::to_string(ij) + "_(" + std::to_string(G.get_edge(ij).first) + "," +
                                                                                         std::to_string(G.get_edge(ij).second) + ")_to_" + std::to_string(kl) + "_(" + std::to_string(H.get_edge(kl).first) + "," +
                                                                                         std::to_string(H.get_edge(kl).second) + ")");
                        }


                        // Constraint node_sub i to k
                        GRBLinExpr le;
                        for (int i = 0; i < n_g; i++) {
                                le = 0;
                                for (int k = 0; k < n_h; k++)
                                        le += node_sub[i][k];
                                model->addConstr(le, GRB_LESS_EQUAL, 1, "Node_Sub");
                        }

                        // Constraint node_sub k to i
                        GRBLinExpr le1;
                        for (int k = 0; k < n_h; k++) {
                                le1 = 0;
                                for (int i = 0; i < n_g; i++)
                                        le1 += node_sub[i][k];
                                model->addConstr(le1, GRB_LESS_EQUAL, 1, "Node_Sub");
                        }

                        // Constraint edge_sub ij to kl
                        // kante ij kann nur auf kl gemapped werden wenn i auf k oder l, oder j auf k oder l gemapped ist
                        GRBLinExpr le2;
                        GRBLinExpr le3;
                        for (int ij = 0; ij < m_g; ij++) {
                                for (int k = 0; k < n_h; k++) {
                                        le2 = 0;
                                        le3 = node_sub[G.get_edge(ij).first][k] + node_sub[G.get_edge(ij).second][k];
                                        for (kl = 0; kl < m_h; kl++) {
                                                if (k == H.get_edge(kl).first || k == H.get_edge(kl).second) {
                                                        le2 += edge_sub[ij][kl];
                                                }
                                        }
                                        model->addConstr(le2, GRB_LESS_EQUAL, le3, "Topological_(" + std::to_string(G.get_edge(ij).first) + "," +
                                                                                   std::to_string(G.get_edge(ij).second) + ")" + "_" + std::to_string(k));
                                }
                        }

                        model->update();

                        for (int i = 0; i < n_g; i++)
                                constant_ += c_ie[i];

                        for (int k = 0; k < n_h; k++)
                                constant_ += c_ek[k];

                        for (int ij = 0; ij < m_g; ij++)
                                constant_ += c_ije[ij];

                        for (kl = 0; kl < m_h; kl++)
                                constant_ += c_ekl[kl];
                        model->addVar(1.0,1.0,constant_,GRB_BINARY,"constant");
                        model->write(opt_.output_lpfile_);


                        model->optimize();

                        opt_.constant_ = constant_;
                        opt_.bbnodecount_ = model->get(GRB_DoubleAttr_NodeCount);
                        opt_.mipgap_ = model->get(GRB_DoubleAttr_MIPGap);
                        opt_.objval_ = model->get(GRB_DoubleAttr_ObjVal);
                        opt_.final_dualbound_ = model->get(GRB_DoubleAttr_ObjBound);
                        opt_.time_ = model->get(GRB_DoubleAttr_Runtime);
                        opt_.status_ = model->get(GRB_IntAttr_Status);
                        opt_.n_vars_full_ = model->get(GRB_IntAttr_NumVars);
                        opt_.n_cons_full_ = model->get(GRB_IntAttr_NumConstrs);
                        opt_.numNZ_ = model->get(GRB_IntAttr_NumNZs);
                        opt_.iterCount_ = model->get(GRB_DoubleAttr_IterCount);


                        // std::cout << "F2 optimization time: " << opt_.time_  << " | MIPGap " << opt_.mipgap_
                        // << " | ObjVal " << opt_.objval_ << " | BBN " << opt_.bbnodecount_ << std::endl;

                        if(opt_.Reoptimize_) {
                                model->reset();
                                model->relax();

                                model->optimize();
                                opt_.lpreltime = model->get(GRB_DoubleAttr_Runtime);
                                opt_.lprelval = model->get(GRB_DoubleAttr_ObjVal);
                        }

                        delete[] node_sub;
                        delete[] edge_sub;
                        delete model;
                        delete env;

                        opt_.output_ = IO::create_output(opt_);

                        return "OK";
                }
                catch (GRBException &e) {
                        std::cout << "Error code = " << e.getErrorCode() << std::endl;
                        std::cout << e.getMessage() << std::endl;
                }
                exit(1);
        }

};


#endif //GEDC_F2_HPP