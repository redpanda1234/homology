import numpy as np
import math


class Simplex:
    """
    A simple simplex class to hold data and stuff. It's not super
    efficient or elegant, basically this just exists so that I don't
    have to use foreach and stuff in TikZ.
    """

    def __init__(self, verts=None, r=1, m=2, offset_percent=None):
        """
        verts is a list of numpy arrays representing positions of
        vertices

        n is the number of vertices

        r is the distance from each vertex to the barycenter

        offset_percent helps fragment the barycentric subdivision
        """
        self.verts = verts
        self.simplices = None
        self.m = m
        self.offset_percent = offset_percent

        if self.verts is None:
            # Initialize vertices. Multiply by r to scale.
            verts = r * np.array(
                [
                    np.array([0, 1]),
                    np.array([-1 * math.sqrt(3) / 2, -0.5]),
                    np.array([math.sqrt(3) / 2, -0.5]),
                ]
            )

            self.verts = verts

            # Define 0 coordinate to be the barycenter
            self.barycenter = np.mean(verts, axis=0)
        else:
            self.barycenter = np.mean(self.verts, axis=0)
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

            v0, v1, v2 = self.verts

            b0_dict = {"0": tuple(v0), "1": tuple(v1), "2": tuple(v2)}

            v01, v02, v12 = 0.5 * (v0 + v1), 0.5 * (v0 + v2), 0.5 * (v1 + v2)

            # Get all 1-simplex barycenters
            b1_dict = {
                "01": tuple(v01),
                "02": tuple(v02),
                "12": tuple(v12),
            }

            simplex_list = []

            # Just get all the indices as strings

            for i in range(3):
                bi = b0_dict[str(i)]

                for j in range(3):
                    if i == j:
                        continue

                    # Eek don't judge me
                    ip, jp = sorted([i, j])
                    bj = b1_dict[str(ip) + str(jp)]

                    newverts = np.array([bi, bj, b])

                    # new barycenter
                    newb = np.mean(newverts, axis=0)

                    if self.offset_percent is not None:
                        offset = self.offset_percent * (newb - b)
                        for vert in newverts:
                            vert += offset

                    simplex_list += [
                        Simplex(
                            verts=newverts, m=m - 1, offset_percent=self.offset_percent
                        )
                    ]
            # print(len(simplex_list))
            self.simplices = simplex_list

            for simplex in self.simplices:
                simplex.subdivide(m - 1)

    def get_TikZ(self):
        out_str = ""
        if self.offset_percent:
            if self.simplices:
                for simplex in self.simplices:
                    out_str += simplex.get_TikZ()
            else:
                for i, vert in enumerate(self.verts):
                    out_str += (
                        f"\\node ({i}) at "
                        + str(tuple([foo.item() for foo in vert]))
                        + "{};\n"
                    )
                out_str += f"\\draw[very thin] (0) -- (1);\n"
                out_str += f"\\draw[densely dotted] (1) -- (2) (0) -- (2);\n"

        elif self.simplices:
            for i, vert in enumerate(self.verts):
                out_str += (
                    f"\\node ({i}) at "
                    + str(tuple([foo.item() for foo in vert]))
                    + "{};\n"
                )

            out_str += "\\draw[very thin] (0) -- (1) -- (2) -- (0) -- cycle;\n"
            for simplex in self.simplices:
                out_str += simplex.get_TikZ()

        else:
            for i, vert in enumerate(self.verts):
                out_str += (
                    f"\\node ({i}) at "
                    + str(tuple([foo.item() for foo in vert]))
                    + "{};\n"
                )
            out_str += f"\\draw[densely dotted] (1) -- (2) (0) -- (2);\n"

        return out_str


def draw(simplex):
    m = simplex.m
    if simplex.offset_percent is not None:
        filename = f"../figures/offset-barycentric-{m}.tex"
    else:
        filename = f"../figures/barycentric-{m}.tex"
    preamble = (
        "\\documentclass{standalone}\n\\usepackage{tikz}\n\\begin{document}\\begin{tikzpicture}[scale="
        + str(m + 1)
        + ",every node/.style={circle, draw=black, fill=white, inner sep=0pt, minimum size=3pt}]\n"
    )
    draw_code = simplex.get_TikZ()
    postamble = "\\end{tikzpicture}\n\\end{document}"
    out_str = preamble + draw_code + postamble
    print(out_str)

    with open(filename, "w") as my_file:
        my_file.write(out_str)

    return


for m in range(1, 5):
    triangle = Simplex(r=1.5, m=m)
    draw(triangle)
