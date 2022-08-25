from ..lattice import Lattice, np


class TriangularLattice(Lattice):
    a = 30
    b = 26

    def __init__(self):
        a, b = self.a, self.b
        Lattice.__init__(self, np.array([[a, 0],
                                         [-a, 0],
                                         [-a/2, b],
                                         [a/2, -b],
                                         [a/2, b],
                                         [-a/2, -b]
                                         ],
                                        dtype=int),
                         [0, np.pi/3, -np.pi/3, 2*np.pi/3, -2*np.pi/3, np.pi],
                         "triangular")

    @classmethod
    def construct_max_dense_core(cls, n_nodes):
        d = 9 + 8 * (n_nodes + 1)
        p = int(np.ceil((-3 + np.sqrt(d)) / 2))
        nodes = []
        nodes += [[cls.a * i, 0] for i in range(p - 1)]
        nodes += [[-cls.a//2 + cls.a * i, -cls.b] for i in range(p)]
        nodes_left = n_nodes - 2 * p + 1
        x0 = 0
        y0 = -cls.b * 2
        for k in range(p - 1):
            if nodes_left == 0:
                break
            n_nodes_in_layer = (p - 1) - k
            n_nodes_in_layer = min(n_nodes_in_layer, nodes_left)
            nodes += [[x0 + cls.a * i, y0 - k * cls.b] for i in range(n_nodes_in_layer)]
            nodes_left -= n_nodes_in_layer
            x0 += cls.a // 2
        return np.array(nodes)
