#ifndef GEDC_IO_HPP
#define GEDC_IO_HPP

#include <iostream>
#include <vector>
#include <sstream>

#include "json.hpp"
//#include "gurobi_c++.h"

#include "base.hpp"
#include "graph.hpp"
#include "options.hpp"

namespace IO {

//        std::string create_output(GRBModel *model, options &opt);

        std::string create_output(const std::string &log_name, const std::string &formulation_name, const std::string &g1_id,
                                  node G1_no_of_nodes, node G1_no_of_edges, const std::string &g2_id, node G2_no_of_nodes,
                                  node G2_no_of_edges, double n_variables, double n_constraints, double bb_node_count,
                                  double ged_val, double time, int NumNZs, double IterCount, double first_root,
                                  double second_root, double constant, int setting, double mipgap);

        std::string create_output(options &opt);

        void setFilename(const std::string &datasetName, const std::string &formulation, options &opt);

        std::vector<std::string> split_string(std::string const &s);

        std::string get_filename(const std::string &path);

        // std::vector<int_label_edge> read_mutag_A(const std::string &mutag_a, const std::string &edge_labels, const std::vector<node> &indicator_vector);

        // std::vector<node> read_indicator(std::string const &indicator_file);

        // void int_label_edgelist_to_files(const std::string &name, const std::vector<int_label_edge> &edges, const std::vector<node> &indicator);

        // std::vector<Edge> read_edgelist(std::string const &filename);

        std::vector<node> get_labels(std::string const &filename);

        std::string indicator_to_string(std::string const &filename);

        graph<std::string, int> read_imdb(int graph_id);

        graph<std::string, int> read_imdb_bonna(int graph_id);

        graph<int, int> read_csv(const std::string &filename, const std::string &dataset);

        graph<int, int> read_edgelist(const std::string &filename, const std::string &dataset);

        graph<int, int> read_edgelist(const std::string &filename, const std::string &dataset);

        graph<int, int> read_csv_imdb(const std::string &filename, const std::string &dataset);




}; // namespace IO

#endif //GEDC_IO_HPP
