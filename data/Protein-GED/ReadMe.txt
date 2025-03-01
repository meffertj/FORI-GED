CONTENTS OF THIS FILE
---------------------
   
 * Reference
 * Authors
 * Database Content

Reference
---------------------
A Graph Database Repository and Performance Evaluation Metrics for Graph Edit Distance, submitted to GBR2015.


Authors
---------------------
Zeina Abu-Aisheh, 3rd year PhD student at François Rabelais University, Tours, France.
Romain Raveaux, Assistant Professor at François Rabelais University, Tours, France.
Jean-Yves Ramel, Professor at François Rabelais University, Tours, France.  


Database Content
---------------------
 - Protein-cost-function

   Two meta parameters are included tau_{vertex} and tau_{edge} where tau_{vertex} denotes a vertex deletion or insertion whereas tau_{edge}
   denotes an edge deletion or insertion. Both tau_{vertex} and tau_{edge} are non-negative parameters. A third meta parameter alpha is integrated to control whether
   the edit operation cost on the vertices or on the edges is more important. For all the subsets of Protein,  tau_{vertex}, tau_{edge} and alpha are set to
   11, 1 and 0.75 respectively. The cost function is written in Java and is found in the Protein-GED folder.
 
 - Protein-low-level-info
   
  The pairwise comparisons of each of the aforementioned subsets is performed. For each comparison, the following information is added:
  1- The name of each pair of graphs d(g_1,g_2).
  2- The number of vertices of each pair of graphs.
  3- The number of edges of each pair of graphs.
  4- The name of a graph edit distance (GED) method P that succeeds at finding the best distance for d(g_1,g_2).
  5- A parameter used for the method BeamSearch GED which represents the size of the stack (1, 10 or 100 in the experiments).
  6- The best distance that was found by method P.
  7- A boolean value that tells whether the found solution is optimal or not.
  8- The classes of graphs g_1 and g_2.
  9- The matching sequence found by method P.
  10- For instance Node:3->4=0.0 means that substituting vertex 3 on graph g_1 with vertex 4 on graph g_2 costs 0.0. Moreover, Edge:1_<>2->1_<>5=4.0 means that substituting
      1_<>2 on graph g_1 with edge 1_<>5 on graph g_2 costs 4.0. On the other hand, Node:3->eps=4.0 means that deleting vertex 3 costs 4.0 while vertex:eps->4=6.0 means that
      inserting vertex 4 costs 6.0.	 
 
 - Protein-subsets
   Contains 4 subsets (Protein-20, Protein-30, Protein-40 and Protein-MIX) from the Protein database, each of which has 10 graphs. 100 pairwise matchings were conducted per subset.
   Protein-MIX represents 10 graphs that were picked up from Protein-20, Protein-30 and Protein-40 and were put in this mixed database.

 - Information about Protein

   In Protein, graphs are constructed from the Protein Data Bank and labeled with their corresponding enzyme class labels from the BRENDA enzyme database. The
   proteins database consists of six classes (EC 1, EC 2, EC 3, EC 4, EC 5, EC 6), which represent proteins out of the six enzyme commission top level hierarchy 
  (EC classes). The proteins are converted into graphs by representing the secondary structure elements of a protein with vertexs and edges of an attributed graph. vertexs are labeled with their 
   type (helix, sheet, or loop) and their amino acid sequence (e.g. TFKEVVRLT). Every vertex is connected with an edge to its three nearest neighbors in space. Edges are labeled with their
   type and the distance they rep- resent in angstroms. 600 proteins are uniformly distributed over 100 classes. The resulting data set is split into a training, 
   a validation and a test set of size 200 each. The classication rate achieved on this data set is 65.5%.

 