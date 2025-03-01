import os

path = "./"
dirs = os.listdir(path)
for file in dirs:
    if file.endswith(".gxl"):
        f = open(file, 'r')

        gxl = ""

        gxl += f.readline()
        templine = f.readline()
        f.close()
        gxl += templine[templine.find("<gxl>"):]
        gxl = gxl[:gxl.find("<gxl>")+5] + "\n" + gxl[gxl.find("<gxl>")+5:]
        gxl = gxl[:gxl.find("</graph>")] + "\n" + gxl[gxl.find("</graph>"):]
        gxl = gxl[:gxl.find("</gxl>")] + "\n" + gxl[gxl.find("</gxl>"):]
        idx = gxl.find("><node")
        while idx != -1:
            gxl = gxl[:idx+1] + "\n" + gxl[idx+1:]
            idx = gxl.find("><node")
        idx = gxl.find("><edge")
        while idx != -1:
            gxl = gxl[:idx+1] + "\n" + gxl[idx+1:]
            idx = gxl.find("><edge")
        f = open(file, 'w')
        f.write(gxl)
        f.close()

