#ifndef GEDC_GRAPH_HPP
#define GEDC_GRAPH_HPP

#include <unordered_map>
#include <limits>
#include <cmath>
#include <numeric>
#include <cassert>

#include "base.hpp"

/**
 *  Node IDs are assumed to be continuous
 *
 * @tparam T node attribute type
 * @tparam U edge attribute type
 */
template<typename T, typename U>
class graph {
        using node = u_int32_t;
        using edge = std::pair<node, node>;

        std::string graph_id_;
        std::string dataset_;
        std::vector<node> node_list_; /** list of nodes */
        std::unordered_map<node, T> node_labels_; /** maps node to label */
        std::unordered_map<node, node> node_map_; /** maps nodeID in this graph to original nodeID in Mutagenicity file */
        std::unordered_map<u_int64_t, U> edge_labels_; /** maps edge to label */
        std::unordered_map<std::string, node> file_label_to_zero_based_;
        std::unordered_map<node, std::string> zero_based_to_file_label_;

        std::vector<std::vector<double>> GEDLIB_node_sub_cost_;


private:
        std::vector<double> GEDLIB_node_del_cost_;
        std::vector<double> GEDLIB_node_ins_cost_;
        std::vector<std::vector<double>> GEDLIB_edge_sub_cost_;
        std::vector<double> GEDLIB_edge_del_cost_;
        std::vector<double> GEDLIB_edge_ins_cost_;


        bool GEDLIB_costs_are_set_ = false;

        /// @brief takes edge as pair of two nodes, returns key for edge_labels_ map (is unaffected by template type)
        inline u_int64_t key(edge edge) {
                return (u_int64_t) edge.first << 32 | (u_int64_t) edge.second;
        }

protected:
        std::vector<std::vector<node>> adjacencylist_;
        std::vector<std::pair<node, node>> edge_list_;
        node n_ = 0; /** number of nodes in graph */
        node m_ = 0; /** number of edges in graph */

public:
        graph<T, U>() = default;


        /**
         * creates graph with node_count nodes, assigning string label "0" to each node
         * @param node_count
         */
        explicit graph<T, U>(node node_count) {
                adjacencylist_.resize(node_count);
                for (node i = 0; i < node_count; i++) {
                        add_node(i, "0");
                }
        }

        [[nodiscard]] node number_of_nodes() const { return n_; }

        [[nodiscard]] node number_of_edges() const { return m_; }


        void set_dataset(std::string name) { dataset_ = std::move(name); }

        [[nodiscard]] std::string get_dataset() const { return dataset_; }


        [[nodiscard]] T get_node_label(node node) const {
                return node_labels_.at(node);
        }

        /// @brief input is key(edge)
        [[nodiscard]] U get_edge_label(edge edge1) const {
                const auto edge_key = key(edge1);
                return edge_labels_.at(edge_key);
        }

        [[nodiscard]] U get_edge_label(idx edgeidx) {
                auto edge_key = key(get_edge(edgeidx));
                return edge_labels_.at(edge_key);
        }


        void set_graph_id(const std::string &graph_id) { graph_id_ = graph_id; }

        [[nodiscard]] std::string get_graph_id() const { return graph_id_; }


        /// we use GEDLIB to calculate edit costs in scip_proteins.cpp, which we simply copy here in order to leave the graph class mostly untouched
        void setGEDLIBeditcosts(std::vector<std::vector<double>> &GEDLIB_node_sub_cost,
                                std::vector<double> &GEDLIB_node_del_cost,
                                std::vector<double> &GEDLIB_node_ins_cost,
                                std::vector<std::vector<double>> &GEDLIB_edge_sub_cost,
                                std::vector<double> &GEDLIB_edge_del_cost,
                                std::vector<double> &GEDLIB_edge_ins_cost) {
                GEDLIB_node_sub_cost_ = GEDLIB_node_sub_cost;
                GEDLIB_node_del_cost_ = GEDLIB_node_del_cost;
                GEDLIB_node_ins_cost_ = GEDLIB_node_ins_cost;
                GEDLIB_edge_sub_cost_ = GEDLIB_edge_sub_cost;
                GEDLIB_edge_del_cost_ = GEDLIB_edge_del_cost;
                GEDLIB_edge_ins_cost_ = GEDLIB_edge_ins_cost;
                GEDLIB_costs_are_set_ = true;
        }

        const std::vector<std::vector<double>> &getGedlibNodeSubCost() const {
                return GEDLIB_node_sub_cost_;
        }

        const std::vector<double> &getGedlibNodeInsCost() const {
                return GEDLIB_node_ins_cost_;
        }

        const std::vector<double> &getGedlibNodeDelCost() const {
                return GEDLIB_node_del_cost_;
        }

        const std::vector<std::vector<double>> &getGedlibEdgeSubCost() const {
                return GEDLIB_edge_sub_cost_;
        }

        const std::vector<double> &getGedlibEdgeDelCost() const {
                return GEDLIB_edge_del_cost_;
        }

        const std::vector<double> &getGedlibEdgeInsCost() const {
                return GEDLIB_edge_ins_cost_;
        }





        // ------------------- template function implementations -------------------------------------------------------


        bool has_node(node node) const {
                if (node >= n_) {
                        return false;
                }
                if (std::find(node_list_.begin(), node_list_.end(), node) != node_list_.end()) {
                        return true;
                }
                return false;
        }

        void add_node(node node, T label) {
                n_ = std::max<u_int32_t>(n_, node + 1);
                adjacencylist_.resize(n_);

                node_list_.push_back(node);
                node_labels_.insert({node, label});
        }


        ///@brief nodeIDs in Mutagenicity file start at 1, so we need to map them to 0-based nodeIDs
        void set_original_nodeID(node new_nodeID, node original_nodeID) {
                if (new_nodeID >= n_)
                        throw std::runtime_error("node index out of bounds");
                node_map_.insert({new_nodeID, original_nodeID});
        }


        node get_original_nodeID(node nodeID) {
                if (nodeID >= n_)
                        throw std::runtime_error("node index out of bounds");
                return node_map_.at(nodeID);
        }


        void set_file_to_zero_indexed(const std::string& file_nodeID, node zero_indexed_node_id) {
                if (zero_indexed_node_id >= n_)
                        throw std::runtime_error("set_file_to_zero_indexed: called before node was added.");
                file_label_to_zero_based_.insert({file_nodeID, zero_indexed_node_id});
                zero_based_to_file_label_.insert({zero_indexed_node_id, file_nodeID});
        }


        node get_zero_indexed_from_file(const std::string& file_ID) {
                if (file_label_to_zero_based_.find(file_ID) == file_label_to_zero_based_.end())
                        throw std::runtime_error("get_zero_indexed_from_file: original file ID wasn't added to map");
                return file_label_to_zero_based_.at(file_ID);
        }

         std::string get_fileID_from_zero_indexed(node nodeID) {
                if (nodeID >= n_ or zero_based_to_file_label_.find(nodeID) == zero_based_to_file_label_.end())
                        throw std::runtime_error("get_fileID_from_zero_indexed: Node ID not found in map");
                return zero_based_to_file_label_.at(nodeID);
        }

        /// @brief returns false if find returns end()
        bool find_file_label(const std::string& file_ID) {
                return !(file_label_to_zero_based_.find(file_ID) == file_label_to_zero_based_.end());
        }


        bool has_edge(node node1, node node2) const {
                if (node1 >= n_ or node2 >= n_) {
                        return false;
                }
                if (adjacencylist_[node1].size() > adjacencylist_[node2].size()) {
                        std::swap(node1, node2);
                }
                if (std::any_of(adjacencylist_[node1].begin(), adjacencylist_[node1].end(), [node2](node n) { return n == node2; })) {
                        return true;
                }


                return false;
        }


        /// @brief does not add parallel-edges
        void add_edge(node node1, node node2, U label) {
                if (node1 > n_ or node2 > n_) {
                        throw std::runtime_error("node doesnt exist, file should contain every node if it is part of an edge");
                }
                if (has_edge(node1, node2)) {
                        return;
                }
                if (node1 > node2) {
                        std::swap(node1, node2);
                }

                adjacencylist_[node1].push_back(node2);
                adjacencylist_[node2].push_back(node1);

                auto this_edge = std::pair<node, node>{node1, node2};
                edge_list_.push_back(this_edge);
                edge_labels_.insert({key(this_edge), label});
                m_++;
        }


        node get_node(idx index) {
                if (index >= node_list_.size())
                        throw std::runtime_error("node index out of bounds");
                return node_list_[index];
        }


        std::pair<node, node> get_edge(idx index) {
                if (index >= edge_list_.size())
                        throw std::runtime_error("edge index out of bounds");
                return edge_list_[index];
        }

        std::vector<std::pair<node, node>> get_edgelist() {
                return this->edge_list_;
        }


        std::vector<node> get_neighbors(node node) {
                if (node >= adjacencylist_.size())
                        throw std::runtime_error("node index out of bounds");
                return adjacencylist_[node];
        }


        node get_degree(node node) {
                if (node >= adjacencylist_.size())
                        throw std::runtime_error("node index out of bounds");
                return adjacencylist_[node].size();
        }

        /// @brief returns V_G \ S
        std::vector<node> nodeset_minus(const std::vector<node> &S) const {
                std::vector<node> ret_vector(n_);
                std::iota(ret_vector.begin(), ret_vector.end(), 0);
                ret_vector.erase(std::remove_if(ret_vector.begin(), ret_vector.end(),
                                                [&S](node entry) { return std::find(S.begin(), S.end(), entry) != S.end(); }));
                return ret_vector;
        }

        /// @brief assumes that edge exists, so only call on "get_neighbors" entries, or after checking has_edge()
        node get_edge_id(node node1, node node2) const {
                if (node1 > node2) {
                        std::swap(node1, node2);
                }
                if (!has_edge(node1, node2)) {
                        throw std::runtime_error("edge doesnt exist");
                }
                for (int ij = 0; ij < edge_list_.size(); ij++) { // i believe edges get added with first node being the smaller id
                        if (edge_list_[ij].first == node1 and edge_list_[ij].second == node2) {
                                return ij;
                        }
                }
                exit(1); // should never happen
        }

        std::vector<node> get_non_neighbors(node node1) {
                if (node1 >= adjacencylist_.size())
                        throw std::runtime_error("node index out of bounds");
                std::vector<node> non_neighbors;
                auto neighbors = adjacencylist_[node1];

                for (int i = 0; i < n_; i++) {
                        if (i == node1)
                                continue;
                        auto it = std::find(neighbors.begin(), neighbors.end(), i);
                        if (it == neighbors.end()) {
                                non_neighbors.push_back(i);
                        }
                }
                return non_neighbors;
        }


        /// @brief mutagenicity and imdb_multi graphs are both graph<string,int> but need different cost_functions
        /// mutagenicity: graph edit distance contest setting 2 for categorical attributes
        /// imdb_multi: graph edit distance contest setting substitution cost: 0, insertion and deletion cost 1
        void cost_function(graph<std::string, int> &G, graph<std::string, int> &H,
                           std::vector<std::vector<double>> &c_ik, std::vector<double> &c_ie, std::vector<double> &c_ek,
                           std::vector<std::vector<double>> &c_ijkl, std::vector<double> &c_ije, std::vector<double> &c_ekl) {
                if (G.get_dataset() == "imdb_multi") {
                        for (int i = 0; i < G.number_of_nodes(); i++) {
                                for (int k = 0; k < H.number_of_nodes(); k++) {
                                        if (G.get_node_label(i) != H.get_node_label(k)) {
                                                c_ik[i][k] = 1;
                                        }
                                }
                        }
                        for (auto &entry: c_ie)
                                entry = 1;
                        for (auto &entry: c_ek)
                                entry = 1;

                        for (int ij = 0; ij < G.number_of_edges(); ij++) {
                                for (int kl = 0; kl < H.number_of_edges(); kl++) {
                                        if (G.get_edge_label(ij) != H.get_edge_label(kl)) {
                                                c_ijkl[ij][kl] = 1;
                                        }
                                }
                        }
                        for (auto &entry: c_ije)
                                entry = 1;
                        for (auto &entry: c_ekl)
                                entry = 1;
                } else if (G.get_dataset() == "mutagenicity" || G.get_dataset() == "aids") {
                        if (not GEDLIB_costs_are_set_) {
                                throw std::runtime_error("Cost function needs to be set manually using \"setGEDLIBeditcosts()\" for the Mutagenicity dataset!");
                        }

                        c_ik = GEDLIB_node_sub_cost_;
                        c_ie = GEDLIB_node_del_cost_;
                        c_ek = GEDLIB_node_ins_cost_;
                        c_ijkl = GEDLIB_edge_sub_cost_;
                        c_ije = GEDLIB_edge_del_cost_;
                        c_ekl = GEDLIB_edge_ins_cost_;
                }
        }

        void cost_function(graph<int, int> &G, graph<int, int> &H,
                           std::vector<std::vector<double>> &c_ik, std::vector<double> &c_ie, std::vector<double> &c_ek,
                           std::vector<std::vector<double>> &c_ijkl, std::vector<double> &c_ije, std::vector<double> &c_ekl) {
                for (int i = 0; i < G.number_of_nodes(); i++) {
                        for (int k = 0; k < H.number_of_nodes(); k++) {
                                if (G.get_node_label(i) != H.get_node_label(k)) {
                                        c_ik[i][k] = 1;
                                }
                        }
                }
                for (auto &entry: c_ie)
                        entry = 1;
                for (auto &entry: c_ek)
                        entry = 1;

                for (int ij = 0; ij < G.number_of_edges(); ij++) {
                        for (int kl = 0; kl < H.number_of_edges(); kl++) {
                                if (G.get_edge_label(ij) != H.get_edge_label(kl)) {
                                        c_ijkl[ij][kl] = 1;
                                }
                        }
                }
                for (auto &entry: c_ije)
                        entry = 1;
                for (auto &entry: c_ekl)
                        entry = 1;

        }


        /// @brief this is cost_function for CMU-HOUSE-A
        /// nodes are labeled with x,y coordinates. node insertion and deletion cost is infinity. edge substitution cost is 0 deletion and insertion cost is 0.5
        void cost_function(graph<std::pair<double, double>, float> &G, graph<std::pair<double, double>, float> &H,
                           std::vector<std::vector<double>> &c_ik, std::vector<double> &c_ie, std::vector<double> &c_ek,
                           std::vector<std::vector<double>> &c_ijkl, std::vector<double> &c_ije,
                           std::vector<double> &c_ekl) {

                for (int i = 0; i < G.number_of_nodes(); i++) {
                        for (int k = 0; k < H.number_of_nodes(); k++) {
                                auto node1 = G.get_node_label(i);
                                auto node2 = H.get_node_label(k);
                                c_ik[i][k] = std::sqrt(std::pow((node1.first - node2.first), 2.0) + std::pow((node1.second - node2.second), 2.0)); // L2 norm
                        }
                }
                for (auto &entry: c_ie)
                        entry = std::numeric_limits<double>::infinity();
                for (auto &entry: c_ek)
                        entry = std::numeric_limits<double>::infinity();

                for (int ij = 0; ij < G.number_of_edges(); ij++) {
                        for (int kl = 0; kl < H.number_of_edges(); kl++) {
                                c_ijkl[ij][kl] = 0;
                        }
                }
                for (auto &entry: c_ije)
                        entry = 0.5;
                for (auto &entry: c_ekl)
                        entry = 0.5;
        }

        /// @brief cost function for protein graphs
        void cost_function(graph<std::pair<int, std::string>, std::tuple<int, int, int>> &G, graph<std::pair<int, std::string>, std::tuple<int, int, int>> &H,
                           std::vector<std::vector<double>> &c_ik, std::vector<double> &c_ie, std::vector<double> &c_ek,
                           std::vector<std::vector<double>> &c_ijkl, std::vector<double> &c_ije,
                           std::vector<double> &c_ekl) {

                if (not GEDLIB_costs_are_set_) {
                        throw std::runtime_error("Cost function needs to be set manually using \"setGEDLIBeditcosts()\" for the Proteins dataset!");
                }

                c_ik = GEDLIB_node_sub_cost_;
                c_ie = GEDLIB_node_del_cost_;
                c_ek = GEDLIB_node_ins_cost_;
                c_ijkl = GEDLIB_edge_sub_cost_;
                c_ije = GEDLIB_edge_del_cost_;
                c_ekl = GEDLIB_edge_ins_cost_;
        }
};

#endif //GEDC_GRAPH_HPP
