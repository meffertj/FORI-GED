//
// Created by julian on 25.11.24.
//

#ifndef GEDC_GEDLIB_COSTS_HPP
#define GEDC_GEDLIB_COSTS_HPP

#define GXL_GEDLIB_SHARED
#include "src/env/ged_env.hpp"

#include "auxiliary/graph.hpp"

template<typename T, typename U>
class getGEDLIBcosts {

private:
    graph<T, U> *graph1_;
    graph<T, U> *graph2_;

public:
    // getGEDLIBcosts(graph<std::string,int> &G, graph<std::string ,int> *graph1)
    getGEDLIBcosts(graph<T, U> *graph1,
                   graph<T, U> *graph2)
            : graph1_(graph1), graph2_(graph2) {}

    /// @brief Sets edit costs directly in the graph class, depending on the Dataset
    void getEditCosts() {
            if (graph1_->get_dataset() != graph2_->get_dataset()) {
                    throw std::runtime_error("Graphs are not from the same dataset!");
            }
            if (graph1_->get_dataset() == "mutagenicity" || graph1_->get_dataset() == "aids") {
                    getChem2EditCosts(graph1_, graph2_);
            } else if (graph1_->get_dataset() == "protein") {
                    getProteinsEditCosts(graph1_, graph2_);
            } else {
                    throw std::runtime_error("Cost function not implemented for this dataset!");
            }
    }

private:
    void getChem2EditCosts(graph<std::pair<int, std::string>, std::tuple<int, int, int>> *graph1, graph<std::pair<int, std::string>, std::tuple<int, int, int>> *graph2) {
            throw std::runtime_error("This is a Protein graph not a Mutagenicity or AIDS graph! -> Called wrong getGEDLIBeditcost() function");
    }

    void getChem2EditCosts(graph<std::string, int> *graph1, graph<std::string, int> *graph2) {
            if (graph1->get_dataset() != "mutagenicity" && graph1->get_dataset() != "aids") {
                    throw std::runtime_error("Graph is not a Mutagenicity neither an AIDS graph!");
            }

            ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
            env.set_edit_costs(ged::Options::EditCosts::CHEM_2);

            std::vector<std::vector<double>> c_ik(graph1->number_of_nodes(), std::vector<double>(graph2->number_of_nodes(), 0));
            std::vector<double> c_ie(graph1->number_of_nodes(), 0);
            std::vector<double> c_ek(graph2->number_of_nodes(), 0);
            std::vector<std::vector<double>> c_ijkl(graph1->number_of_edges(), std::vector<double>(graph2->number_of_edges(), 0));
            std::vector<double> c_ije(graph1->number_of_edges(), 0);
            std::vector<double> c_ekl(graph2->number_of_edges(), 0);


            for (int i = 0; i < graph1->number_of_nodes(); i++) {
                    std::map<std::string, std::string> gxlLabelNode1;
                    const auto labelNode1 = graph1->get_node_label(i);
                    gxlLabelNode1.insert({"chem", labelNode1});
                    c_ie[i] = env.node_del_cost(gxlLabelNode1);

                    for (int k = 0; k < graph2->number_of_nodes(); k++) {
                            std::map<std::string, std::string> gxlLabelNode2;
                            const auto labelNode2 = graph2->get_node_label(k);
                            gxlLabelNode2.insert({"chem", labelNode2});

                            c_ik[i][k] = env.node_rel_cost(gxlLabelNode1, gxlLabelNode2);
                    }
            }
            for (int k = 0; k < graph2->number_of_nodes(); k++) {
                    std::map<std::string, std::string> gxlLabelNode2;
                    const auto labelNode2 = graph2->get_node_label(k);
                    gxlLabelNode2.insert({"chem", labelNode2});

                    c_ek[k] = env.node_ins_cost(gxlLabelNode2);

            }

            // für edges brauchen wir frequency, type0 und type1 (abhängig von frequency)
            for (int ij = 0; ij < graph1->number_of_edges(); ij++) {
                    std::map<std::string, std::string> gxlLabelEdge1;
                    const auto labelEdge1 = graph1->get_edge_label(ij);
                    gxlLabelEdge1.insert({"valence", to_string(labelEdge1)});

                    c_ije[ij] = env.edge_del_cost(gxlLabelEdge1);

                    for (int kl = 0; kl < graph2->number_of_edges(); kl++) {
                            std::map<std::string, std::string> gxlLabelEdge2;
                            const auto labelEdge2 = graph2->get_edge_label(kl);
                            gxlLabelEdge2.insert({"valence", to_string(labelEdge2)});

                            c_ijkl[ij][kl] = env.edge_rel_cost(gxlLabelEdge1, gxlLabelEdge2);
                    }
            }
            for (int kl = 0; kl < graph2->number_of_edges(); kl++) {
                    std::map<std::string, std::string> gxlLabelEdge2;
                    const auto labelEdge2 = graph2->get_edge_label(kl);
                    gxlLabelEdge2.insert({"valence", to_string(labelEdge2)});

                    c_ekl[kl] = env.edge_ins_cost(gxlLabelEdge2);
            }

            graph1->setGEDLIBeditcosts(c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);
            graph2->setGEDLIBeditcosts(c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);
    }


    void getProteinsEditCosts(graph<std::string, int> *graph1, graph<std::string, int> *graph2) {
            throw std::runtime_error("This is a Mutagenicity graph not a Proteins graph! -> Called wrong getGEDLIBeditcost() function");
    }

    void getProteinsEditCosts(graph<std::pair<int, std::string>, std::tuple<int, int, int>> *graph1, graph<std::pair<int, std::string>, std::tuple<int, int, int>> *graph2) {
            if (graph1->get_dataset() != "protein") {
                    throw std::runtime_error("Graph is not a Proteins graph!");
            }

            ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
            env.set_edit_costs(ged::Options::EditCosts::PROTEIN);

            std::vector<std::vector<double>> c_ik(graph1->number_of_nodes(), std::vector<double>(graph2->number_of_nodes(), 0));
            std::vector<double> c_ie(graph1->number_of_nodes(), 0);
            std::vector<double> c_ek(graph2->number_of_nodes(), 0);
            std::vector<std::vector<double>> c_ijkl(graph1->number_of_edges(), std::vector<double>(graph2->number_of_edges(), 0));
            std::vector<double> c_ije(graph1->number_of_edges(), 0);
            std::vector<double> c_ekl(graph2->number_of_edges(), 0);


            for (int i = 0; i < graph1->number_of_nodes(); i++) {
                    std::map<std::string, std::string> gxlLabelNode1;
                    const auto labelNode1 = graph1->get_node_label(i);
                    gxlLabelNode1.insert({"type", to_string(labelNode1.first)});
                    gxlLabelNode1.insert({"sequence", labelNode1.second});
                    c_ie[i] = env.node_del_cost(gxlLabelNode1);
                    for (int k = 0; k < graph2->number_of_nodes(); k++) {
                            std::map<std::string, std::string> gxlLabelNode2;
                            const auto labelNode2 = graph2->get_node_label(k);
                            gxlLabelNode2.insert({"type", to_string(labelNode2.first)});
                            gxlLabelNode2.insert({"sequence", labelNode2.second});

                            c_ik[i][k] = env.node_rel_cost(gxlLabelNode1, gxlLabelNode2);
                    }
            }
            for (int k = 0; k < graph2->number_of_nodes(); k++) {
                    std::map<std::string, std::string> gxlLabelNode2;
                    const auto labelNode2 = graph2->get_node_label(k);
                    gxlLabelNode2.insert({"type", to_string(labelNode2.first)});
                    gxlLabelNode2.insert({"sequence", labelNode2.second});

                    c_ek[k] = env.node_ins_cost(gxlLabelNode2);
            }

            // für edges brauchen wir frequency, type0 und type1 (abhängig von frequency)
            for (int ij = 0; ij < graph1->number_of_edges(); ij++) {
                    std::map<std::string, std::string> gxlLabelEdge1;
                    const auto labelEdge1 = graph1->get_edge_label(ij);
                    gxlLabelEdge1.insert({"frequency", to_string(std::get<0>(labelEdge1))});
                    gxlLabelEdge1.insert({"type0", to_string(std::get<1>(labelEdge1))});
                    gxlLabelEdge1.insert({"type1", to_string(std::get<2>(labelEdge1))});

                    c_ije[ij] = env.edge_del_cost(gxlLabelEdge1);
                    for (int kl = 0; kl < graph2->number_of_edges(); kl++) {
                            std::map<std::string, std::string> gxlLabelEdge2;
                            const auto labelEdge2 = graph2->get_edge_label(kl);
                            gxlLabelEdge2.insert({"frequency", to_string(std::get<0>(labelEdge2))});
                            gxlLabelEdge2.insert({"type0", to_string(std::get<1>(labelEdge2))});
                            gxlLabelEdge2.insert({"type1", to_string(std::get<2>(labelEdge2))});

                            c_ijkl[ij][kl] = env.edge_rel_cost(gxlLabelEdge1, gxlLabelEdge2);
                    }
            }
            for (int kl = 0; kl < graph2->number_of_edges(); kl++) {
                    std::map<std::string, std::string> gxlLabelEdge2;
                    const auto labelEdge2 = graph2->get_edge_label(kl);
                    gxlLabelEdge2.insert({"frequency", to_string(std::get<0>(labelEdge2))});
                    gxlLabelEdge2.insert({"type0", to_string(std::get<1>(labelEdge2))});
                    gxlLabelEdge2.insert({"type1", to_string(std::get<2>(labelEdge2))});

                    c_ekl[kl] = env.edge_ins_cost(gxlLabelEdge2);
            }

            graph1->setGEDLIBeditcosts(c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);
            graph2->setGEDLIBeditcosts(c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);
    }
};


#endif //GEDC_GEDLIB_COSTS_HPP
