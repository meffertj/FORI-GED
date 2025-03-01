#include <cstdint>
#include <iostream>
#include <fstream>
#include <stack>


#include "auxiliary/graph.hpp"

#include "auxiliary/GXLGraphReader.hpp"

/**
 * in mutagenicity graphs, nodeIDs start from one
 * nodes are labeled with a string, edges are labeled with an int
 * @param path
 * @return
 */
graph<std::string, int> GXLGraphReader::read_mutagenicity(const std::string &path) {

        std::ifstream graphFile(path);
        if (!graphFile.is_open()) {
                throw std::runtime_error("could not open file " + path);
        }
        std::string line;

        auto syntaxCheck = [](const std::vector<std::string> &tokens) {
                std::stack<std::string> openTags;
                std::string cur;
                for (auto token: tokens) {
                        if (token[0] != '/') {
                                auto end = token.find(' ');
                                cur = token.substr(0, end);
                                if (cur == "node" || cur == "edge" || cur == "graph" || cur == "attr" ||
                                    cur == "string" || cur == "int" ||
                                    cur == "Float" || cur == "Double") {
                                        openTags.push(cur);
                                }
                        } else {
                                if (openTags.empty())
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered opening tag without closing tag" + token);
                                if (openTags.top() != token.substr(1))
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered closing tag for " + token.substr(1) +
                                                " but expected closing tag for " +
                                                openTags.top());
                                openTags.pop();
                        }
                }
        };

        /**
         * returns vector of tokens (contents between < and >)
         */
        auto tokenize = [](std::string &line, idx i) {
                std::vector<std::string> tokens;
                idx start;
                idx end = i;
                while (end < line.size()) {
                        if (line[end] == '<') {
                                ++end;
                                start = end;
                                while (line[end] != '>') {
                                        if (line[end] == '<') {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered another < before closing > " + line);
                                        } else if (end == line.size() - 1) {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before closing > " + line);
                                        }
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                                ++end;
                        } else if (line[end] != '<' && end < line.size()) {
                                start = end;
                                while (line[end] != '<') {
                                        if (line[end] == '>')
                                                throw std::runtime_error("In tokenize: encountered > before opening < " + line);
                                        else if (end == line.size() - 1)
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before opening < " + line);
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                        }
                }
                return tokens;
        };

        /**
         * assumes attributed nodes
         * supports lines of the form: <node id="XYZ"><attr name="XYZ"><type>XYZ</type></attr></node>
         * matching the format found in e.g. Mutagenicity dataset
         * (multiple attributes per node as in CMU, or Proteins is not supported)
         * double as attribute get casted to float
         *
         */
        auto parseNode = [&](graph<std::string, int> &G, std::vector<std::string> tokens) {
                std::string original_nodeID;
                std::string attrName;
                std::string attrType = tokens[2];
                std::string nodeAttr = tokens[3];

                if (tokens[1].find("attr") == std::string::npos)
                        throw std::runtime_error("In parseNode: line containing node element does not have attr element");

                auto startID = tokens[0].find('"', 0);
                auto endID = tokens[0].find('"', startID + 1);
                original_nodeID = tokens[0].substr(startID + 1, endID - startID - 1);

                auto startName = tokens[1].find('"', 0);
                auto endName = tokens[1].find('"', startName + 1);
                attrName = tokens[1].substr(startName + 1, endName - startName - 1);

                if (attrName != "chem")
                        throw std::runtime_error("In parseNode: attribute name is not chem");

                node new_nodeID = std::stoul(original_nodeID) - 1;
                G.add_node(new_nodeID, nodeAttr);
                G.set_original_nodeID(new_nodeID, std::stoul(original_nodeID));
        };

        /**
         * wenn man das mapping durch die GED rekonstruieren will muss man die nodeIDs aus dem file aus den hier gesetzten rekonstruieren
         */
        auto parseEdge = [&](graph<std::string, int> &G, std::vector<std::string> tokens) {
                node original_uID;
                node original_vID;
                std::string attrName;
                int attrValue = std::stoi(tokens[3]);
                std::string attrType = tokens[2];
                if (attrType != "int")
                        throw std::runtime_error("In parseEdge: Edge attribute type should be int in Mutagenicity dataset");

                auto start_attr_name = tokens[1].find('"', 0);
                auto end_attr_name = tokens[1].find('"', start_attr_name + 1);
                attrName = tokens[1].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                if (attrName != "valence") {
                        throw std::runtime_error("In parseEdge: Edge attribute name should be valence in Mutagenicity dataset");
                }

                auto startNode1 = tokens[0].find('"', 0);
                auto endNode1 = tokens[0].find('"', startNode1 + 1);
                auto startNode2 = tokens[0].find('"', endNode1 + 1);
                auto endNode2 = tokens[0].find('"', startNode2 + 1);
                original_uID = std::stoul(tokens[0].substr(startNode1 + 1, endNode1 - startNode1 - 1));
                original_vID = std::stoul(tokens[0].substr(startNode2 + 1, endNode2 - startNode2 - 1));

                auto node1 = original_uID - 1;
                auto node2 = original_vID - 1;
                G.add_edge(node1, node2, attrValue);
        };

        auto parseHead = [](std::ifstream &graphFile, std::string &line, graph<std::string, int> &G) {
                std::string graphID;

                std::getline(graphFile, line);
                if (line.find("xml") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected xml tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("gxl") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected gxl tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("graph") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected graph tag, file format not .gxl");
                } else if (line.find("id=") != std::string::npos) {
                        auto start = line.find('"', 0);
                        auto end = line.find('"', start + 1);
                        graphID = line.substr(start + 1, end - start - 1);
                        G.set_graph_id(graphID);
                } else {
                        throw std::runtime_error("In parseHead: expected graph id, file format not .gxl");
                }
        };

        /**
         * maybe implement converter that fits each gxl element into one line (e.g. CMU housing data set splits these into multiple lines)
         */
        auto parseFile = [&]() {
                graph<std::string, int> G;
                parseHead(graphFile, line, G);
                for (std::string tmpline; getline(graphFile, tmpline);) {
                        if (tmpline.find("node") != std::string::npos) {
                                auto tokens = tokenize(tmpline, 0);
                                syntaxCheck(tokens);
                                parseNode(G, tokens);
                        } else if (tmpline.find("edge") != std::string::npos) {
                                auto tokens = tokenize(tmpline, 0);
                                syntaxCheck(tokens);
                                parseEdge(G, tokens);
                        }
                }
                return G;

        };
        graph<std::string, int> G = parseFile();
        G.set_dataset("mutagenicity");
        graphFile.close();

        return G;
}

/**
 * nodes are unlabeled, edge attributes are euclidian distance between nodes
 * reads the CMU dataset from the gxl file
 * @param path
 * @return
 */
graph<std::pair<double, double>, float> GXLGraphReader::read_CMU(const std::string &path) {

        std::ifstream graphFile(path);
        if (!graphFile.is_open()) {
                throw std::runtime_error("In read_CMU: could not open file " + path);
        }
        std::string line;

        auto syntaxCheck = [](const std::vector<std::string> &tokens) {
                std::stack<std::string> openTags;
                std::string cur;
                for (auto token: tokens) {
                        if (token[0] != '/') {
                                auto end = token.find(' ');
                                cur = token.substr(0, end);
                                if (cur == "node" || cur == "edge" || cur == "graph" || cur == "attr" ||
                                    cur == "string" || cur == "int" ||
                                    cur == "Float" || cur == "Double") {
                                        openTags.push(cur);
                                }
                        } else {
                                if (openTags.empty())
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered opening tag without closing tag" + token);
                                if (openTags.top() != token.substr(1))
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered closing tag for " + token.substr(1) +
                                                " but expected closing tag for " +
                                                openTags.top());
                                openTags.pop();
                        }
                }
        };

        /**
         * returns vector of tokens (contents between < and >)
         */
        auto tokenize = [](std::string &line, idx i) {
                std::vector<std::string> tokens;
                idx start;
                idx end = i;
                while (end < line.size()) {
                        if (line[end] == '<') {
                                ++end;
                                start = end;
                                while (line[end] != '>') {
                                        if (line[end] == '<') {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered another < before closing > " + line);
                                        } else if (end == line.size() - 1) {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before closing > " + line);
                                        }
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                                ++end;
                        } else if (line[end] != '<' && end < line.size()) {
                                start = end;
                                while (line[end] != '<') {
                                        if (line[end] == '>')
                                                throw std::runtime_error("In tokenize: encountered > before opening < " + line);
                                        else if (end == line.size() - 1)
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before opening < " + line);
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                        }
                }
                return tokens;
        };

        /**
         * nodes are labeled with x,y coordinates
         * supports lines of the form:
         * <node id="1">
         * <attr name="x"><Double>197.048387</Double></attr><attr name="y"><Double>342.887097</Double></attr></node>
         * matching the format found in CMU dataset
         *
         */
        auto parseNode = [&](graph<std::pair<double, double>, float> &G, std::vector<std::string> tokens) {
                std::string original_nodeID;
                std::string attrName1;
                std::string attrName2;
                std::string attrType1 = tokens[2];
                std::string attrType2 = tokens[7];
                double nodeAttr1 = std::stod(tokens[3]);
                double nodeAttr2 = std::stod(tokens[8]);

                if (tokens[1].find("attr") == std::string::npos) {
                        throw std::runtime_error("In parseNode: line containing node element does not have attr element for first attribute");
                }
                if (tokens[6].find("attr") == std::string::npos) {
                        throw std::runtime_error("In parseNode: line containing node element does not have attr element for second attribute");
                }
                if (attrType1 != "Double")
                        throw std::runtime_error("In parseNode: Node attribute type should be Double in CMU dataset (first node attribute)");
                if (attrType2 != "Double")
                        throw std::runtime_error("In parseNode: Node attribute type should be Double in CMU dataset (second node attribute)");

                auto startID = tokens[0].find('"', 0);
                auto endID = tokens[0].find('"', startID + 1);
                original_nodeID = tokens[0].substr(startID + 1, endID - startID - 1);

                auto startName = tokens[1].find('"', 0);
                auto endName = tokens[1].find('"', startName + 1);
                attrName1 = tokens[1].substr(startName + 1, endName - startName - 1);

                if (attrName1 != "x")
                        throw std::runtime_error("In parseNode: attribute name is not x");

                startName = tokens[6].find('"', 0);
                endName = tokens[6].find('"', startName + 1);
                attrName2 = tokens[6].substr(startName + 1, endName - startName - 1);

                if (attrName2 != "y")
                        throw std::runtime_error("In parseNode: attribute name is not y");

                node new_nodeID = std::stoul(original_nodeID) - 1;
                G.add_node(new_nodeID, std::pair<double, double>{nodeAttr1, nodeAttr2});
                G.set_original_nodeID(new_nodeID, std::stoul(original_nodeID));
        };

        /**
         * wenn man das mapping durch die GED rekonstruieren will muss man die nodeIDs aus dem file aus den hier gesetzten rekonstruieren
         */
        auto parseEdge = [&](graph<std::pair<double, double>, float> &G, std::vector<std::string> tokens) {
                node original_uID;
                node original_vID;
                std::string attrName;
                float attrValue = std::stof(tokens[3]);
                std::string attrType = tokens[2];
                if (attrType != "Float")
                        throw std::runtime_error("In parseEdge: Edge attribute type should be Float in CMU dataset");

                auto start_attr_name = tokens[1].find('"', 0);
                auto end_attr_name = tokens[1].find('"', start_attr_name + 1);
                attrName = tokens[1].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                if (attrName != "dist") {
                        throw std::runtime_error("In parseEdge: Edge attribute name should be dist in CMU dataset");
                }

                auto startNode1 = tokens[0].find('"', 0);
                auto endNode1 = tokens[0].find('"', startNode1 + 1);
                auto startNode2 = tokens[0].find('"', endNode1 + 1);
                auto endNode2 = tokens[0].find('"', startNode2 + 1);
                original_uID = std::stoul(tokens[0].substr(startNode1 + 1, endNode1 - startNode1 - 1));
                original_vID = std::stoul(tokens[0].substr(startNode2 + 1, endNode2 - startNode2 - 1));

                auto node1 = original_uID - 1;
                auto node2 = original_vID - 1;
                G.add_edge(node1, node2, attrValue);
        };

        auto parseHead = [](graph<std::pair<double, double>, float> &G, std::ifstream &graphFile, std::string &line) {
                std::string graphID;

                std::getline(graphFile, line);
                if (line.find("xml") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected xml tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("gxl") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected gxl tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("graph") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected graph tag, file format not .gxl");
                } else if (line.find("id=") != std::string::npos) {
                        auto start = line.find('"', 0);
                        auto end = line.find('"', start + 1);
                        graphID = line.substr(start + 1, end - start - 1);
                        G.set_graph_id(graphID);
                } else {
                        throw std::runtime_error("In parseHead: expected graph id, file format not .gxl");
                }
        };

        /**
         * in CMU node id and edge id are in seperate lines than attributes, thats why we erase the \r
         */
        auto parseFile = [&]() {
                graph<std::pair<double, double>, float> G;
                parseHead(G, graphFile, line);
                for (std::string tmpline; getline(graphFile, tmpline);) {
                        if (tmpline.find("node") != std::string::npos) {
                                std::string concatenate = tmpline;
                                getline(graphFile, tmpline);
                                concatenate += tmpline;
                                concatenate.erase(std::remove(concatenate.begin(), concatenate.end(), '\r'), concatenate.end());

                                auto tokens = tokenize(concatenate, 0);
                                syntaxCheck(tokens);
                                parseNode(G, tokens);
                        } else if (tmpline.find("edge") != std::string::npos) {
                                std::string concatenate = tmpline;
                                getline(graphFile, tmpline);
                                concatenate += tmpline;
                                concatenate.erase(std::remove(concatenate.begin(), concatenate.end(), '\r'), concatenate.end());

                                auto tokens = tokenize(concatenate, 0);
                                syntaxCheck(tokens);
                                parseEdge(G, tokens);
                        }
                }
                return G;

        };

        graph<std::pair<double, double>, float> G = parseFile();
        G.set_dataset("CMU");
        graphFile.close();

        return G;
}

graph<std::pair<int, std::string>, std::tuple<int, int, int>> GXLGraphReader::read_Proteins(const std::string &path) {

        std::ifstream graphFile(path);
        if (!graphFile.is_open()) {
                throw std::runtime_error("In read_Proteins: could not open file " + path);
        }
        std::string line;

        auto syntaxCheck = [](const std::vector<std::string> &tokens) {
                std::stack<std::string> openTags;
                std::string cur;
                for (auto token: tokens) {
                        if (token[0] != '/') {
                                auto end = token.find(' ');
                                cur = token.substr(0, end);
                                if (cur == "node" || cur == "edge" || cur == "graph" || cur == "attr" ||
                                    cur == "string" || cur == "int" ||
                                    cur == "Float" || cur == "Double" || cur == "double") {
                                        openTags.push(cur);
                                }
                        } else {
                                if (openTags.empty())
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered opening tag without closing tag" + token);
                                if (openTags.top() != token.substr(1))
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered closing tag for " + token.substr(1) +
                                                " but expected closing tag for " +
                                                openTags.top());
                                openTags.pop();
                        }
                }
        };

        /**
         * returns vector of tokens (contents between < and >)
         */
        auto tokenize = [](std::string &line, idx i) {
                std::vector<std::string> tokens;
                idx start;
                idx end = i;
                while (end < line.size()) {
                        if (line[end] == '<') {
                                ++end;
                                start = end;
                                while (line[end] != '>') {
                                        if (line[end] == '<') {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered another < before closing > " + line);
                                        } else if (end == line.size() - 1) {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before closing > " + line);
                                        }
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                                ++end;
                        } else if (line[end] != '<' && end < line.size()) {
                                start = end;
                                while (line[end] != '<') {
                                        if (line[end] == '>')
                                                throw std::runtime_error("In tokenize: encountered > before opening < " + line);
                                        else if (end == line.size() - 1)
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before opening < " + line);
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                        }
                }
                return tokens;
        };

        /**
         * nodes are labeled with amino acid sequences
         * supports lines of the form:
         * <node id="1"><attr name="type"><int>0</int></attr><attr name="aaLength"><int>11</int></attr><attr name="sequence"><int>NVLIEDLKWRG</int></attr></node>
         * <node id="2"><attr name="type"><int>0</int></attr><attr name="aaLength"><int>12</int></attr><attr name="sequence"><int>DEQGIEDLLNKE</int></attr></node>
         * matching the format found in the "Proteins" dataset
         */
        auto parseNode = [&](graph<std::pair<int, std::string>, std::tuple<int, int, int>> &G, std::vector<std::string> tokens) {
                std::string original_nodeID;
                std::string attrName;

                std::string attrType = tokens[2];
                if (attrType != "int")
                        throw std::runtime_error("In parseNode: node attr in protein should be int");

                if (tokens[1].find("attr") == std::string::npos)
                        throw std::runtime_error("In parseNode: line containing node element does not have attr element");


                std::pair<int, std::string> nodeAttr = {std::stoi(tokens[3]), tokens[13]};


                auto startID = tokens[0].find('"', 0);
                auto endID = tokens[0].find('"', startID + 1);
                original_nodeID = tokens[0].substr(startID + 1, endID - startID - 1);

                auto startName = tokens[1].find('"', 0);
                auto endName = tokens[1].find('"', startName + 1);
                attrName = tokens[1].substr(startName + 1, endName - startName - 1);

                if (attrName != "type")
                        throw std::runtime_error("In parseNode: attribute name is not type");

                node new_nodeID = std::stoul(original_nodeID) - 1;
                G.add_node(new_nodeID, nodeAttr);
                G.set_original_nodeID(new_nodeID, std::stoul(original_nodeID));
        };

        /**
         * wenn man das mapping durch die GED rekonstruieren will muss man die nodeIDs aus dem file aus den hier gesetzten rekonstruieren
         */
        auto parseEdge = [&](graph<std::pair<int, std::string>, std::tuple<int, int, int>> &G, std::vector<std::string> tokens) {
                node original_uID;
                node original_vID;
                std::string attrName;

                std::tuple<int, int, int> attrValue;

                std::string attrType = tokens[2];
                if (attrType != "int")
                        throw std::runtime_error("In parseEdge: Edge attribute type should be int in Protein dataset");


                auto start_attr_name = tokens[1].find('"', 0);
                auto end_attr_name = tokens[1].find('"', start_attr_name + 1);
                attrName = tokens[1].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                if (attrName != "frequency") {
                        throw std::runtime_error("In parseEdge: First edge attributes name should be frequency in Protein dataset");
                }


                start_attr_name = tokens[6].find('"', 0);
                end_attr_name = tokens[6].find('"', start_attr_name + 1);
                attrName = tokens[6].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                if (attrName != "type0") {
                        throw std::runtime_error("In parseEdge: Second edge attributes name should be type0 in Protein dataset");
                }
                if (tokens[7] != "double")
                        throw std::runtime_error("In parseEdge: Edge attribute type of type0 should be double in Protein dataset");


                start_attr_name = tokens[11].find('"', 0);
                end_attr_name = tokens[11].find('"', start_attr_name + 1);
                attrName = tokens[11].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                if (attrName != "distance0") {
                        throw std::runtime_error("In parseEdge: Third edge attributes name should be distance0 in Protein dataset");
                }
                if (tokens[12] != "double")
                        throw std::runtime_error("In parseEdge: Edge attribute type of distance0 should be double in Protein dataset");


                // if frequency is 2 then there are two sets of edge attributes
                if (std::stoi(tokens[3]) == 2) {
                        start_attr_name = tokens[16].find('"', 0);
                        end_attr_name = tokens[16].find('"', start_attr_name + 1);
                        attrName = tokens[16].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                        if (attrName != "type1") {
                                throw std::runtime_error("In parseEdge: Fourth edge attributes name should be type1 in Protein dataset");
                        }
                        if (tokens[17] != "double")
                                throw std::runtime_error("In parseEdge: Edge attribute type of type1 should be int in Protein dataset");


                        start_attr_name = tokens[21].find('"', 0);
                        end_attr_name = tokens[21].find('"', start_attr_name + 1);
                        attrName = tokens[21].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                        if (attrName != "distance1") {
                                throw std::runtime_error("In parseEdge: Fifth edge attributes name should be distance1 in Protein dataset");
                        }
                        if (tokens[22] != "double")
                                throw std::runtime_error("In parseEdge: Edge attribute type should of distance1 be double in Protein dataset");

                }

                auto startNode1 = tokens[0].find('"', 0);
                auto endNode1 = tokens[0].find('"', startNode1 + 1);
                auto startNode2 = tokens[0].find('"', endNode1 + 1);
                auto endNode2 = tokens[0].find('"', startNode2 + 1);
                original_uID = std::stoul(tokens[0].substr(startNode1 + 1, endNode1 - startNode1 - 1));
                original_vID = std::stoul(tokens[0].substr(startNode2 + 1, endNode2 - startNode2 - 1));

                auto node1 = original_uID - 1;
                auto node2 = original_vID - 1;

                if (std::stoi(tokens[3]) == 1) {
                        attrValue = {std::stoi(tokens[3]), std::stoi(tokens[8]), -1};
                } else if (std::stoi(tokens[3]) == 2) {
                        attrValue = {std::stoi(tokens[3]), std::stoi(tokens[8]), std::stoi(tokens[18])};
                }

                G.add_edge(node1, node2, attrValue);
        };

        auto parseHead = [](std::ifstream &graphFile, std::string &line, graph<std::pair<int, std::string>, std::tuple<int, int, int>> &G) {
                std::string graphID;

                std::getline(graphFile, line);
                if (line.find("xml") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected xml tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("gxl") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected gxl tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("graph") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected graph tag, file format not .gxl");
                } else if (line.find("id=") != std::string::npos) {
                        auto start = line.find('"', 0);
                        auto end = line.find('"', start + 1);
                        graphID = line.substr(start + 1, end - start - 1);
                        G.set_graph_id(graphID);
                } else {
                        throw std::runtime_error("In parseHead: expected graph id, file format not .gxl");
                }
        };

        /**
         *
         */
        auto parseFile = [&]() {
                graph<std::pair<int, std::string>, std::tuple<int, int, int>> G;
                parseHead(graphFile, line, G);
                for (std::string tmpline; getline(graphFile, tmpline);) {
                        if (tmpline.find("node") != std::string::npos) {
                                auto tokens = tokenize(tmpline, 0);
                                syntaxCheck(tokens);
                                parseNode(G, tokens);
                        } else if (tmpline.find("edge") != std::string::npos) {
                                auto tokens = tokenize(tmpline, 0);
                                syntaxCheck(tokens);
                                parseEdge(G, tokens);
                        }
                }
                return G;

        };

        graph<std::pair<int, std::string>, std::tuple<int, int, int>> G = parseFile();
        G.set_dataset("protein");
        graphFile.close();

        return G;
}



graph<std::string, int> GXLGraphReader::read_AIDS(const std::string &path) {

       std::ifstream graphFile(path);
        if (!graphFile.is_open()) {
                throw std::runtime_error("could not open file " + path);
        }
        std::string line;

        auto syntaxCheck = [](const std::vector<std::string> &tokens) {
                std::stack<std::string> openTags;
                std::string cur;
                for (auto token: tokens) {
                        if (token[0] != '/') {
                                auto end = token.find(' ');
                                cur = token.substr(0, end);
                                if (cur == "node" || cur == "edge" || cur == "graph" || cur == "attr" ||
                                    cur == "string" || cur == "int" ||
                                    cur == "float" || cur == "Double") {
                                        openTags.push(cur);
                                }
                        } else {
                                if (openTags.empty())
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered opening tag without closing tag" + token);
                                if (openTags.top() != token.substr(1))
                                        throw std::runtime_error(
                                                "In syntaxCheck: encountered closing tag for " + token.substr(1) +
                                                " but expected closing tag for " +
                                                openTags.top());
                                openTags.pop();
                        }
                }
        };

        /**
         * returns vector of tokens (contents between < and >)
         */
        auto tokenize = [](std::string &line, idx i) {
                std::vector<std::string> tokens;
                idx start;
                idx end = i;
                while (end < line.size()) {
                        if (line[end] == '<') {
                                ++end;
                                start = end;
                                while (line[end] != '>') {
                                        if (line[end] == '<') {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered another < before closing > " + line);
                                        } else if (end == line.size() - 1) {
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before closing > " + line);
                                        }
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                                ++end;
                        } else if (line[end] != '<' && end < line.size()) {
                                start = end;
                                while (line[end] != '<') {
                                        if (line[end] == '>')
                                                throw std::runtime_error("In tokenize: encountered > before opening < " + line);
                                        else if (end == line.size() - 1)
                                                throw std::runtime_error(
                                                        "In tokenize: encountered end of line before opening < " + line);
                                        ++end;
                                }
                                tokens.push_back(line.substr(start, end - start));
                        }
                }
                return tokens;
        };

        /**
         * assumes attributed nodes
         * supports lines of the form: <node id="XYZ"><attr name="XYZ"><type>XYZ</type></attr></node>
         * matching the format found in e.g. Mutagenicity dataset
         * (multiple attributes per node as in CMU, or Proteins is not supported)
         * double as attribute get casted to float
         *
         */
        auto parseNode = [&](graph<std::string, int> &G, std::vector<std::string> tokens) {
                std::string original_nodeID;
                std::string attrName;
                std::string attrType = tokens[7];
                std::string nodeAttr = tokens[8];

                if (tokens[6].find("attr") == std::string::npos)
                        throw std::runtime_error("In parseNode: line containing node element does not have attr element");

                auto startID = tokens[0].find('"', 0);
                auto endID = tokens[0].find('"', startID + 1);
                original_nodeID = tokens[0].substr(startID + 1, endID - startID - 1);
                original_nodeID.erase(std::remove(original_nodeID.begin(), original_nodeID.end(), '_'), original_nodeID.end());
                auto startName = tokens[6].find('"', 0);
                auto endName = tokens[6].find('"', startName + 1);
                attrName = tokens[6].substr(startName + 1, endName - startName - 1);

                if (attrName != "chem")
                        throw std::runtime_error("In parseNode: attribute name is not chem");

                node new_nodeID = std::stoul(original_nodeID) - 1;
                G.add_node(new_nodeID, nodeAttr);
                G.set_original_nodeID(new_nodeID, std::stoul(original_nodeID));
        };

        /**
         * wenn man das mapping durch die GED rekonstruieren will muss man die nodeIDs aus dem file aus den hier gesetzten rekonstruieren
         */
        auto parseEdge = [&](graph<std::string, int> &G, std::vector<std::string> tokens) {
                node original_uID;
                node original_vID;
                std::string toclean_original_uID;
                std::string toclean_original_vID;

                std::string attrName;
                int attrValue = std::stoi(tokens[3]);
                std::string attrType = tokens[2];
                if (attrType != "int")
                        throw std::runtime_error("In parseEdge: Edge attribute type should be int in Mutagenicity dataset");

                auto start_attr_name = tokens[1].find('"', 0);
                auto end_attr_name = tokens[1].find('"', start_attr_name + 1);
                attrName = tokens[1].substr(start_attr_name + 1, end_attr_name - start_attr_name - 1);
                if (attrName != "valence") {
                        throw std::runtime_error("In parseEdge: Edge attribute name should be valence in Mutagenicity dataset");
                }

                auto startNode1 = tokens[0].find('"', 0);
                auto endNode1 = tokens[0].find('"', startNode1 + 1);
                auto startNode2 = tokens[0].find('"', endNode1 + 1);
                auto endNode2 = tokens[0].find('"', startNode2 + 1);
                toclean_original_uID = tokens[0].substr(startNode1 + 1, endNode1 - startNode1 - 1);
                toclean_original_uID.erase(std::remove(toclean_original_uID.begin(), toclean_original_uID.end(), '_'),
                        toclean_original_uID.end());

                original_uID = std::stoul(toclean_original_uID);

                toclean_original_vID = tokens[0].substr(startNode2 + 1, endNode2 - startNode2 - 1);
                toclean_original_vID.erase(std::remove(toclean_original_vID.begin(), toclean_original_vID.end(), '_'),
                        toclean_original_vID.end());

                original_vID = std::stoul(toclean_original_vID);

                auto node1 = original_uID - 1;
                auto node2 = original_vID - 1;
                G.add_edge(node1, node2, attrValue);
        };

        auto parseHead = [](std::ifstream &graphFile, std::string &line, graph<std::string, int> &G) {
                std::string graphID;

                std::getline(graphFile, line);
                if (line.find("xml") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected xml tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("gxl") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected gxl tag, file format not .gxl");
                }
                std::getline(graphFile, line);
                if (line.find("graph") == std::string::npos) {
                        throw std::runtime_error("In parseHead: expected graph tag, file format not .gxl");
                } else if (line.find("id=") != std::string::npos) {
                        auto start = line.find('"', 0);
                        auto end = line.find('"', start + 1);
                        graphID = line.substr(start + 1, end - start - 1);
                        G.set_graph_id(graphID);
                } else {
                        throw std::runtime_error("In parseHead: expected graph id, file format not .gxl");
                }
        };

        /**
         * maybe implement converter that fits each gxl element into one line (e.g. CMU housing data set splits these into multiple lines)
         */
        auto parseFile = [&]() {
                graph<std::string, int> G;
                parseHead(graphFile, line, G);
                for (std::string tmpline; getline(graphFile, tmpline);) {
                        if (tmpline.find("node") != std::string::npos) {
                                auto tokens = tokenize(tmpline, 0);
                                syntaxCheck(tokens);
                                parseNode(G, tokens);
                        } else if (tmpline.find("edge") != std::string::npos) {
                                auto tokens = tokenize(tmpline, 0);
                                syntaxCheck(tokens);
                                parseEdge(G, tokens);
                        }
                }
                return G;

        };
        graph<std::string, int> G = parseFile();
        G.set_dataset("aids");
        graphFile.close();

        return G;
}



