
class MSCNode :

    def __init__(self):
        self.arcs = []

    def read_from_line(self, line) :
        tmplist = line.split()
        self.cellid = int(tmplist[0])
        self.index = int(tmplist[1])
        self.value = float(tmplist[2])
        self.boundary = int(tmplist[3])
        self.xy = (float(tmplist[4]), float(tmplist[5]))

    def add_arc(self, arc) :
        self.arcs.append(arc)

class MSCArc :

    def __init(self):
        self.nodes = []

    def __group_xy(self, lst):
        for i in range(0, len(lst), 3):
            yield tuple(lst[i:i+2])

    def read_from_line(self, line) :
        tmplist = line.split()
        self.index = int(tmplist[0])
        self.node_ids = [int(tmplist[1]), int(tmplist[2])]
        self.line = [i for i in self.__group_xy([float(i) for i in tmplist[6:]])]

class MSC :

    def __init__(self):
        self.nodes = {}
        self.arcs = []

    def read_from_file(self, fname_base) :
        nodesname = fname_base + ".nodes.txt"
        arcsname = fname_base + ".arcs.txt"
        node_file = open(nodesname, "r")
        nodes_lines = node_file.readlines()
        node_file.close()
        for l in nodes_lines :
            n = MSCNode()
            n.read_from_line(l)
            self.nodes[n.cellid] = n
        arcs_file = open(arcsname, "r")
        arcs_lines = arcs_file.readlines()
        arcs_file.close()
        for l in arcs_lines :
            a = MSCArc()
            a.read_from_line(l)
            n1 = self.nodes[a.node_ids[0]]
            n2 = self.nodes[a.node_ids[1]]
            n1.add_arc(a)
            n2.add_arc(a)
            a.nodes = [n1, n2]
            self.arcs.append(a)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:
        fname = "Screen711x706.raw"

    test_msc = MSC()
    test_msc.read_from_file(fname)
    for cellid, n in test_msc.nodes.items() :
        print(str(n.cellid) + " == " + str(cellid))
    for a in test_msc.arcs :
        print( str(a.index) + ": " + str(a.line[0:4] )+ "...")
