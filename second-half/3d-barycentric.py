import numpy as np
import math

class Simplex:
    """
    A simple simplex class to hold data and stuff. It's not super
    efficient or elegant, basically this just exists so that I don't
    have to use foreach and stuff in TikZ.
    """
    def __init__(self, verts=None, r=1, m=2, color=None):
        """
        verts is a list of numpy arrays representing positions of
        vertices

        n is the number of vertices

        r is the distance from each vertex to the barycenter
        """

        self.verts = verts
        self.simplices = None
        self.m = m
        self.color=color

        if self.verts is None:
            # Initialize vertices. Multiply by r to scale.
            verts = r * np.array([np.array([0,1,0]),
                                  np.array([-1*math.sqrt(3)/2, -.5,0]),
                                  np.array([math.sqrt(3)/2, -.5,0]),
                                  np.array([0,0,1.5])])

            self.verts = verts

            # Define 0 coordinate to be the barycenter
            self.barycenter = np.mean(verts, axis=0)
        else:
            self.barycenter = np.mean(self.verts,axis=0)
        if m:
            self.subdivide(m)

    def subdivide(self, m):
        """
        m is the depth to subdivide to
        """
        if m == 0:
            return
        else:
            color_dict = {
                0 : "red",
                1 : "green",
                2 : "blue",
                3 : "orange"
            }

            b = self.barycenter

            # Technically, I could automate the generation of all of the
            # following by using itertools and stuff --- although this would
            # bring me great satisfaction, I seem to be unable to visualize
            # things in 4D space, so sadly, it's hard to justify. Maybe I'll
            # come back to it?

            v0, v1, v2, v3 = self.verts

            # This might be more verbose than it's worth, but I want to be able
            # to select the proper maximal simplices without doing too much
            # hard-coding of it. b0_dict associates to each i = 0,...,3 the
            # i-th 0-simplex barycenter (the i-th vertex)
            b0_dict = {
                "0" : tuple(v0),
                "1" : tuple(v1),
                "2" : tuple(v2),
                "3" : tuple(v3)
            }

            v01, v02, v12 = .5*(v0 + v1), .5*(v0 + v2), .5*(v1 + v2)
            v03, v13, v23 = .5*(v0 + v3), .5*(v1 + v3), .5*(v2 + v3)

            # Get all 1-simplex barycenters
            b1_dict = {
                "01" : tuple(v01),
                "02" : tuple(v02),
                "03" : tuple(v03),
                "12" : tuple(v12),
                "13" : tuple(v13),
                "23" : tuple(v23)
            }

            v012 = .25 * (v0 + v1 + v2)
            v013 = .25 * (v0 + v1 + v3)
            v023 = .25 * (v0 + v2 + v3)
            v123 = .25 * (v1 + v2 + v3)

            # Get all 2-simplex barycenters
            b2_dict = {
                "012" : tuple(v012),
                "013" : tuple(v013),
                "023" : tuple(v023),
                "123" : tuple(v123)
            }

            simplex_list = []

            # Just get all the indices as strings

            for i in range(4):
                bi = b0_dict[str(i)]

                for j in range(4):
                    if i == j:
                        continue

                    # Eek don't judge me
                    ip, jp = sorted([i,j])
                    bj = b1_dict[str(ip) + str(jp)]

                    for k in range(4):
                        if i == k or j == k:
                            continue
                        # aaaaaaaaaa
                        ip, jp, kp = sorted([i,j,k])

                        bk = b2_dict[str(ip) + str(jp) + str(kp)]

                        newverts = np.array([bi, bj, bk, b])
                        if self.color is not None:
                            newcolor = self.color
                        else:
                            newcolor = color_dict[i]
                        simplex_list += [Simplex(
                            verts=newverts,
                            m=m-1,
                            color=newcolor
                        )]
            # print(len(simplex_list))
            self.simplices = simplex_list

            for simplex in self.simplices:
                simplex.subdivide(m-1)

    def get_TikZ_sets(self, nodes, s_edges, d_edges):
        """
        We want to store nodes and edges in dictionaries so that we can more
        easily eliminate redundancy during the drawing process and stuff.

        nodes is a set of node tuples

        s_edges is the set of edges to draw with solid lines
        d_edges is the set of edges to draw with dotted lines
        """
        # Something something hashable types
        v0, v1, v2, v3 = [tuple(vert) for vert in self.verts]

        # Add these nodes to the global set
        nodes |= {v0, v1, v2, v3}

        if self.simplices:
            # Add the edges we need to draw later. A single tuple represents
            # one edge to draw.

            # These edges are solid because they were present in the previous
            # subdivision
            # c = self.color
            # if self.color:
            #     s_edges |= {(v0, v1, c), (v0, v2, c), (v0, v3, c), (v1, v2, c),
            #                 (v1, v3, c), (v2, v3, c)}
            s_edges |= {(v0,v1), (v0,v2), (v0,v3), (v1,v2), (v1,v3), (v2,v3)}


            for simplex in self.simplices:
                # Recurse
                nodes, s_edges, d_edges = simplex.get_TikZ_sets(
                    nodes,
                    s_edges,
                    d_edges
                )
        else:
            # Only draw lines that weren't there already (everything but v0 --
            # v1). v2 is the barycenter of the triangular face, v3 that of the
            # tetrahedron.
            d_edges |= {(v0, v2), (v0, v3), (v1, v2), (v1, v3), (v2, v3)}

        return nodes, s_edges, d_edges



    def get_TikZ(self):
        out_str = ""
        nodes, s_edges, d_edges = self.get_TikZ_sets(set(), set(), set())

        print(s_edges)

        s_str = "% draw the solid edges\n\\draw[very thin]"
        for va, vb in s_edges:
            s_str += f" {va} -- {vb}"
        s_str += ";\n"

        d_str = "% draw the dotted edges\n\\draw[densely dotted]"
        for va, vb in d_edges:
            d_str += f" {va} -- {vb}"
        d_str += ";\n"

        n_str = "% draw the nodes in\n"
        for i, node in enumerate(nodes):
            n_str += f"\\node ({i}) at {node} " + "{};\n"

        out_str += s_str + d_str + n_str

        return out_str

def draw(simplex):
    m = simplex.m
    filename = f"./figures/3d-barycentric-{m}.tex"
    preamble = "\\documentclass{standalone}\n\\usepackage{tikz}\n\\usepackage{tikz-3dplot}\n\\begin{document}\n\\tdplotsetmaincoords{70}{80}\n\\begin{tikzpicture}[tdplot_main_coords,scale=" + str(m+1) + ",every node/.style={circle, draw=black, fill=white, inner sep=0pt, minimum size=3pt}]\n"
    draw_code = simplex.get_TikZ()

    v0, v1, v2, v3 = simplex.verts
    v01, v02, v12 = .5*(v0 + v1), .5*(v0 + v2), .5*(v1 + v2)
    v03, v13, v23 = .5*(v0 + v3), .5*(v1 + v3), .5*(v2 + v3)

    # I know how awful this is. I really just want to get color drawing working
    # and be done with it though.
    v0, v1, v2, v3, v01, v02, v12, v03, v13, v23 = [tuple(v) for v in [v0, v1,
                                                                       v2, v3,
                                                                       v01,
                                                                       v02,
                                                                       v12,
                                                                       v03,
                                                                       v13, v23]]

    draw_code += f"\\fill[red!30!white, opacity=.5] {v0} -- {v1} -- {v2} -- {v0} -- cycle;\n"
    draw_code += f"\\fill[blue!30!white, opacity=.5] {v0} -- {v1} -- {v3} -- {v0} -- cycle;\n"
    draw_code += f"\\fill[green!30!white, opacity=.5] {v0} -- {v2} -- {v3} -- {v0} -- cycle;\n"
    draw_code += f"\\fill[orange!30!white, opacity=.5] {v1} -- {v2} -- {v3} -- {v1} -- cycle;\n"

    b = simplex.barycenter


    postamble = "\\end{tikzpicture}\n\\end{document}"
    out_str = preamble + draw_code + postamble
    print(out_str)

    with open(filename, "w") as my_file:
        my_file.write(out_str)

    return

tetrahedron = Simplex(r=1.5,m=1)
draw(tetrahedron)
