import os
import networkit as nk

def readbins(basefolder:str):
    z = open(basefolder+"/imdb_graphs_statistics.csv", 'w')
    z.write("Graph,Folder,V,E\n")
    for type in ["testcsv", "traincsv"]:
        for file in os.listdir(basefolder+"/"+type):
            graphname = file[:file.find(".csv")]
            elr = nk.graphio.EdgeListReader(',', 0)
            graph = elr.read(basefolder+"/"+type+"/"+file)
            z.write(f"{graphname},{type},{graph.numberOfNodes()},{graph.numberOfEdges()}\n")
    z.close()



basefolder = "."
readbins(basefolder)