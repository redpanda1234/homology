import numpy as np
import math

class Simplex:
    """
    A simple simplex class to hold data and stuff. It's not super
    efficient or elegant, basically this just exists so that I don't
    have to use foreach and stuff in TikZ.
    """
    def __init__(self, verts=None, r=1, m=2):
        """
        verts is a list of numpy arrays representing positions of
        vertices

        n is the number of vertices

        r is the distance from each vertex to the barycenter
        """
        self.verts = verts
        self.simplices = None
        self.m = m

        if self.verts is None:
            # Initialize vertices. Multiply by r to scale.
            verts = r * np.array([np.array([0,1]),
                                  np.array([-1*math.sqrt(3)/2, -.5]),
                                  np.array([math.sqrt(3)/2, -.5])])

            self.verts = verts

            # Define 0 coordinate to be the barycenter
            self.barycenter = np.array([0.0,0.0])
        else:
            # print("\n\n")
            self.barycenter = np.mean(self.verts,axis=0)
            # print("barycenter is ", self.barycenter, "for verts", self.verts)
        if m:
            self.subdivide(m)

    def subdivide(self, m):
        """
        m is the depth to subdivide to
        """
        if m == 0:
            return
        else:
            b = self.barycenter
            # print(b)
            v0, v1, v2 = self.verts[0], self.verts[1], self.verts[2]
            v01, v02, v12 = .5*(v0 + v1), .5*(v0 + v2), .5*(v1 + v2)

            s0_01 = Simplex(verts=np.array([v0, v01, b]),m=m-1)
            s0_02 = Simplex(verts=np.array([v0, v02, b]),m=m-1)

            s1_01 = Simplex(verts=np.array([v1, v01, b]),m=m-1)
            s1_12 = Simplex(verts=np.array([v1, v12, b]),m=m-1)

            s2_02 = Simplex(verts=np.array([v2, v02, b]),m=m-1)
            s2_12 = Simplex(verts=np.array([v2, v12, b]),m=m-1)

            self.simplices = [s0_01, s0_02, s1_01, s1_12, s2_02, s2_12]
            for simplex in self.simplices:
                simplex.subdivide(m-1)

    def get_TikZ(self):
        out_str = ""
        if self.simplices:
            for i, vert in enumerate(self.verts):
                out_str += f"\\node ({i}) at " + str(tuple(vert)) + "{};\n"

            out_str += "\\draw[very thin] (0) -- (1) -- (2) -- (0) -- cycle;\n"
            # if self.simplices[0].simplices is None:
            #     for i, vert in enumerate(self.verts):
            #         out_str += f"\\node ({i}) at " + str(tuple(vert)) + "{};\n"

            #     out_str += "\\draw[very thin] (0) -- (1) -- (2) -- cycle;\n"
            for simplex in self.simplices:
                out_str += simplex.get_TikZ()

        else:
            # print(self.verts)
            for i, vert in enumerate(self.verts):
                out_str += f"\\node ({i}) at " + str(tuple(vert)) + "{};\n"
            out_str += f"\\draw[densely dotted] (1) -- (2) (0) -- (2);\n"

        return out_str

def draw(simplex):
    m = simplex.m
    filename = f"./figures/barycentric-{m}.tex"
    preamble = "\\documentclass{standalone}\n\\usepackage{tikz}\n\\begin{document}\\begin{tikzpicture}[scale=" + str(m+1) + ",every node/.style={circle, draw=black, fill=white, inner sep=0pt, minimum size=3pt}]\n"
    draw_code = simplex.get_TikZ()
    postamble = "\\end{tikzpicture}\n\\end{document}"
    out_str = preamble + draw_code + postamble
    print(out_str)

    with open(filename, "w") as my_file:
        my_file.write(out_str)

    return

for m in range(1,6):
    triangle = Simplex(r=1.5, m=m)
    draw(triangle)
