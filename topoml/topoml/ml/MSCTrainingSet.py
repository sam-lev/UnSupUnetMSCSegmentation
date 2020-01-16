class SupervisedMSC:
    def __init__(self, msc=None, geomsc=None):
        self.msc = msc
        self.geomsc = geomsc
        self.arcs = None
        self.arc_map = None

    def map_labeling(self, msc=None, geomsc=None, labeled_data=None):
        if not msc:
            self.msc = msc
        if not geomsc:
            self.geomsc = geomsc
            self.msc = geomsc
        # This ensures smaller arcs take precedence in the event of
        # pixel overlap
        sorted_arcs = sorted(self.msc.arcs, key=lambda arc: len(arc.line))
        arc_points = []
        self.arcs = []
        for arc in sorted_arcs:
            if (not self.use_ridge_arcs and is_ridge_arc(arc, self.msc)) or (
                    not self.use_valley_arcs and is_valley_arc(arc, self.msc)
            ):
                continue
            index = tuple(arc.node_ids) + (len(arc.line),)#make_arc_id(arc)
            arc_points.extend(arc.line)
            self.arc_map.extend([index] * len(arc.line))
            self.arcs.append(arc)

