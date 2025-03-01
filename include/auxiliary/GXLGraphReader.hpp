#ifndef NETWORKIT_IO_GXL_GRAPH_READER_HPP_
#define NETWORKIT_IO_GXL_GRAPH_READER_HPP_

#include "base.hpp"
#include "graph.hpp"

class GXLGraphReader {
public:
        GXLGraphReader() = default;

        /// @brief mutagenicity graphs have categorical integer labels on edges and nodes
        static graph<std::string, int> read_mutagenicity(const std::string &path);


        /// @brief nodes are labeled with x,y coordinates, and delaunay triangulation as edges
        static graph<std::pair<double, double>, float> read_CMU(const std::string &path);

        /**
         * excerpt from "IAM Graph Database Repository for Graph Based Pattern Recognition and Machine Learning"
         * The protein data set consists of graphs representing proteins originally used
            in [3]. The graphs are constructed from the Protein Data Bank [28] and labeled with their corresponding enzyme class labels from the BRENDA enzyme
            database [29]. The proteins database consists of six classes (EC 1, EC 2, EC 3,
            EC 4, EC 5, EC 6), which represent proteins out of the six enzyme commission
            top level hierarchy (EC classes). The proteins are converted into graphs by representing the secondary structure elements of a protein with nodes and edges of
            an attributed graph. Nodes are labeled with their type (helix, sheet, or loop) and
            their amino acid sequence (e.g. TFKEVVRLT). Every node is connected with
            an edge to its three nearest neighbors in space. Edges are labeled with their type
            and the distance they represent in angstroms. In Fig. 7 six images of proteins of
            all six classes are given.
         * @param path
         * @return
         */
        /// @brief nodes are labeled with aminoacid sequences, multiple edge labels: frequency:{1,2}, type0:{1,2,3,4}, distance0:{double}, type1:{4,5}, distance1:{int}
        static graph<std::pair<int, std::string>, std::tuple<int, int, int>> read_Proteins(const std::string &path);

        static graph<std::string, int> read_AIDS(const std::string &path);
};

#endif // NETWORKIT_IO_GXL_GRAPH_READER_HPP_