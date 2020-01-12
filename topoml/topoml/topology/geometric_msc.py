# Standard library imports

# Third party imports

# Local application imports
import os

print("%%%%%%%%%%%%%%%%%%%%%% ", os.getcwd())

class MSCNode:
    def __init__(self):
        self.arcs = []

    def read_from_line(self, line):
        tmplist = line.split(",")
        self.cellid = int(tmplist[0])
        self.index = None #int(tmplist[1])
        #self.value = float(tmplist[2])
        #self.boundary = int(tmplist[3])
        self.xy = (float(tmplist[1]), float(tmplist[2]))

    def add_arc(self, arc):
        self.arcs.append(arc)


class MSCArc:
    def __init(self):
        self.nodes = []

    def __group_xy(self, lst):
        for i in range(0, len(lst), 2):
            yield tuple(lst[i : i + 2])

    def read_from_line(self, line):
        tmplist = line.split(",")
        self.index = int(tmplist[0])
        self.node_ids = [int(tmplist[1]), int(tmplist[2])]
        self.line = [
            i for i in self.__group_xy([float(i) for i in tmplist[3:]])
        ] #read the rest of the the points in the arc as xy tuples

class GeoMSC:
    def __init__(self):
        self.nodes = {}
        self.arcs = []

    def read_from_file(self, fname_base):
        nodesname = fname_base + ".nodes.txt"
        arcsname = fname_base + ".arcs.txt"
        node_file = open(nodesname, "r")
        nodes_lines = node_file.readlines()
        node_file.close()
        for l in nodes_lines:
            n = MSCNode()
            n.read_from_line(l)
            self.nodes[n.cellid] = n
        arcs_file = open(arcsname, "r")
        arcs_lines = arcs_file.readlines()
        arcs_file.close()
        for l in arcs_lines:
            a = MSCArc()
            a.read_from_line(l)
            n1 = self.nodes[a.node_ids[0]]
            n2 = self.nodes[a.node_ids[1]] 
            n1.index = a.index
            n2.index = a.index
            n1.add_arc(a)
            n2.add_arc(a)
            a.nodes = [n1, n2]
            self.arcs.append(a)
