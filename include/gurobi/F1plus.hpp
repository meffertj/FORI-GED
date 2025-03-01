#ifndef GEDC_F1plus_HPP
#define GEDC_F1plus_HPP

#include <random>

#include "gurobi_c++.h"

#include "auxiliary/graph.hpp"
#include "root_relaxation_cb.hpp"
#include "no_branching_cb.hpp"
#include "auxiliary/io.hpp"
#include "auxiliary/options.hpp"


template<typename T, typename U>
class F1plus {
private:
        bool testing_ = false;
        bool relax_ = false;
        double constant_ = 0.0;
        char VAR_TYPE_ = GRB_BINARY;
        int threads_ = 1;
        int seed_ = 1;

public:
        F1plus() = default;

        explicit F1plus(options &opt) : opt_(opt) {}


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
                        GRBVar **node_sub = new GRBVar *[n_g];
                        for(int i = 0; i < n_g; i++) {
                          node_sub[i] = new GRBVar[n_h];
                        }
                        GRBVar *node_del = new GRBVar[n_g];
                        GRBVar *node_ins = new GRBVar [n_h];
                        GRBVar **edge_sub = new GRBVar *[m_g];
                        for(int i = 0; i < m_g; i++) {
                          edge_sub[i] = new GRBVar[m_h];
                        }
                        GRBVar *edge_del = new GRBVar[m_g];
                        GRBVar *edge_ins = new GRBVar[m_h];

                        env = new GRBEnv();
                        GRBModel *model = new GRBModel(*env);
                        model->set(GRB_StringAttr_ModelName, "F1plus_"+opt_.dataset_name_ + "_"+opt_.G_id_ + "_" +opt_.H_id_);
                        model->set(GRB_StringParam_LogFile, opt_.log_fname_ + ".log");
                        model->set(GRB_IntParam_Threads, opt_.threads_);
                        if (timelimit_ != 0) {
                                model->set(GRB_DoubleParam_TimeLimit, opt_.timelimit_);
                        }

                        model->set(GRB_IntParam_Seed, opt_.seed_);

                        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

                        int kl = 0;


                        std::vector<std::vector<double>> c_ik(n_g, std::vector<double>(n_h, 0));
                        std::vector<double> c_ie(n_g, 0);
                        std::vector<double> c_ek(n_h, 0);
                        std::vector<std::vector<double>> c_ijkl(m_g, std::vector<double>(m_h, 0));
                        std::vector<double> c_ije(m_g, 0);
                        std::vector<double> c_ekl(m_h, 0);


                        // depending on the graph type this chooses the corresponding cost function
                        G.cost_function(G, H, c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);


                        // für jeden knoten in n_g gibt es n_h viele mögliche mappings
                        for (int i = 0; i < n_g; i++) {
                                for (int k = 0; k < n_h; k++) {
                                        node_sub[i][k] = model->addVar(0, 1, c_ik[i][k], GRB_BINARY, "map_node_" + std::to_string(i) + "_to_" + std::to_string(k));
                                }
                        }

                        for (int i = 0; i < n_g; i++) {
                                node_del[i] = model->addVar(0, 1, c_ie[i], GRB_BINARY, "delete_node_" + std::to_string(i));
                                if (c_ie[i] == std::numeric_limits<float>::infinity()) {
                                        node_del[i].set(GRB_DoubleAttr_UB, 0);
                                }
                        }

                        for (int k = 0; k < n_h; k++) {
                                node_ins[k] = model->addVar(0, 1, c_ek[k], GRB_BINARY, "insert_node_" + std::to_string(k));
                                if (c_ek[k] == std::numeric_limits<float>::infinity()) {
                                        node_ins[k].set(GRB_DoubleAttr_UB, 0);
                                }
                        }


                        // für jede kante in m_g gibt es m_h viele mögliche mappings
                        for (int ij = 0; ij < m_g; ij++) {
                                for (kl = 0; kl < m_h; kl++)
                                        edge_sub[ij][kl] = model->addVar(0, 1, c_ijkl[ij][kl],
                                                                         GRB_BINARY, "map_edge_" + std::to_string(ij) + "_" + std::to_string(G.get_edge(ij).first) + "-" +
                                                                                         std::to_string(G.get_edge(ij).second) + "_to_" + std::to_string(kl) + "_" + std::to_string(H.get_edge(kl).first) + "-" +
                                                                                         std::to_string(H.get_edge(kl).second));
                        }

                        for (int ij = 0; ij < m_g; ij++)
                                edge_del[ij] = model->addVar(0, 1, c_ije[ij], GRB_BINARY, "delete_edge_" + std::to_string(ij));

                        for (kl = 0; kl < m_h; kl++)
                                edge_ins[kl] = model->addVar(0, 1, c_ekl[kl], GRB_BINARY, "insert_edge_" + std::to_string(kl));


                        // Constraint node i gets mapped to k, or gets deleted
                        GRBLinExpr le;
                        for (int i = 0; i < n_g; i++) {
                                le = 0;
                                le += node_del[i];
                                for (int k = 0; k < n_h; k++)
                                        le += node_sub[i][k];
                                model->addConstr(le, GRB_EQUAL, 1, "node_" + std::to_string(i) + "_mapped_or_deleted");
                        }

                        // Constraint node k gets mapped to i, or gets inserted
                        GRBLinExpr le1;
                        for (int k = 0; k < n_h; k++) {
                                le1 = 0;
                                le1 += node_ins[k];
                                for (int i = 0; i < n_g; i++)
                                        le1 += node_sub[i][k];
                                model->addConstr(le1, GRB_EQUAL, 1, "node_" + std::to_string(k) + "_mapped_or_inserted");
                        }


                        // Constraint edge gets mapped or deleted
                        GRBLinExpr le2;
                        for (int ij = 0; ij < m_g; ij++) {
                                le2 = 0;
                                le2 += edge_del[ij];
                                for (kl = 0; kl < m_h; kl++)
                                        le2 += edge_sub[ij][kl];
                                model->addConstr(le2, GRB_EQUAL, 1, "edge_id" + std::to_string(ij) + "_" +
                                                                    std::to_string(G.get_edge(ij).first) + "_to_" + std::to_string(G.get_edge(ij).second) + "_mapped_or_deleted");
                        }

                        // Constraint edge gets mapped or inserted
                        GRBLinExpr le3;
                        for (kl = 0; kl < m_h; kl++) {
                                le3 = 0;
                                le3 += edge_ins[kl];
                                for (int ij = 0; ij < m_g; ij++)
                                        le3 += edge_sub[ij][kl];
                                model->addConstr(le3, GRB_EQUAL, 1, "edge_id" + std::to_string(kl) + "_" +
                                                                    std::to_string(H.get_edge(kl).first) + "_to_" + std::to_string(H.get_edge(kl).second) + "_mapped_or_inserted");
                        }


                        GRBLinExpr le8;
                        GRBLinExpr le9;
                        for (int ij = 0; ij < m_g; ij++) {
                                for (int k = 0; k < n_h; k++) {
                                        le8 = 0;
                                        le9 = node_sub[G.get_edge(ij).first][k] + node_sub[G.get_edge(ij).second][k];
                                        for (kl = 0; kl < m_h; kl++) {
                                                if (k == H.get_edge(kl).first || k == H.get_edge(kl).second) {
                                                        le8 += edge_sub[ij][kl];
                                                }
                                        }
                                        model->addConstr(le8, GRB_LESS_EQUAL, le9, "Topological_(" + std::to_string(G.get_edge(ij).first) + "," +
                                                                                   std::to_string(G.get_edge(ij).second) + ")" + "_" + std::to_string(k));
                                }
                        }
                        GRBLinExpr le10;
                        GRBLinExpr le11;
                        for (int ij = 0; ij < m_g; ij++) {
                                for (int k = 0; k < n_h; k++) {
                                        le10 = 0;
                                        le11 = node_sub[G.get_edge(ij).first][k] + node_sub[G.get_edge(ij).second][k];
                                        for (kl = 0; kl < m_h; kl++) {
                                                if (k == H.get_edge(kl).first || k == H.get_edge(kl).second) {
                                                        le10 += edge_sub[ij][kl];
                                                }
                                        }
                                        model->addConstr(le10, GRB_LESS_EQUAL, le11, "Topological_(" + std::to_string(G.get_edge(ij).first) + "," +
                                                                                   std::to_string(G.get_edge(ij).second) + ")" + "_" + std::to_string(k));
                                }
                        }

                        // Constraint edge_sub kl to ij
                        GRBLinExpr le12;
                        GRBLinExpr le13;
                        for (int kl = 0; kl < m_h; kl++) {
                              for (int i = 0; i < n_g; i++){
                                le12 = 0;
                                le13 = node_sub[i][H.get_edge(kl).first] + node_sub[i][H.get_edge(kl).second];
                                for (int ij = 0; ij < m_g; ij++) {
                                  if(i == G.get_edge(ij).first || i == G.get_edge(ij).second){
                                    le12 += edge_sub[ij][kl];
                                  }
                                }
                                      model->addConstr(le12, GRB_LESS_EQUAL, le13, "Topological_" + std::to_string(i) + "_(" +
                                                                                         std::to_string(H.get_edge(kl).first)
                                                                               + "," + std::to_string(H.get_edge(kl).second) + ")");
                              }

                        }
                        model->write(opt_.output_lpfile_.c_str());

                        model->optimize();

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

//                        std::cout << "F1+ optimization time: " << opt_.time_  << " | MIPGap " << opt_.mipgap_
//                       << " | ObjVal " << opt_.objval_ << " | BBN " << opt_.bbnodecount_ << std::endl;
                        if(opt_.Reoptimize_) {
                                model->reset();
                                model->relax();

                                model->optimize();
                                opt_.lpreltime = model->get(GRB_DoubleAttr_Runtime);
                                opt_.lprelval = model->get(GRB_DoubleAttr_ObjVal);
                        }

                        delete[] node_sub;
                        delete[] node_del;
                        delete[] node_ins;
                        delete[] edge_sub;
                        delete[] edge_del;
                        delete[] edge_ins;
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


#endif //GEDC_F1plus_HPP