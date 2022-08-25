from ..lattice import Lattice
from ...math import EquationSolver
import numpy as np


class FCC3DLattice(Lattice):
    def __init__(self):
        Lattice.__init__(self, np.array([[1,1,0],
                                         [1,-1,0],
                                         [-1,1,0],
                                         [-1,-1,0],
                                         [0,1,1],
                                         [0,-1,1],
                                         [1,0,1],
                                         [-1,0,1],
                                         [0, 1, -1],
                                         [0, -1, -1],
                                         [1, 0, -1],
                                         [-1, 0, -1]
                                         ]),
                         [0, np.pi / 2, -np.pi / 2, np.pi],
                         "3D_FCC")

    @classmethod
    def construct_max_dense_core(cls, n_nodes):
        def build_layer(n):
            l_nodes = []
            for k in range(n):
                l_n = k + 1
                x_range = np.arange(0, 2 * l_n, 2) - k
                l_nodes += np.vstack(([x_range],
                                      np.full((1, l_n), -k))).T.tolist()

            if n > 1:
                app_part = np.array(l_nodes[:-n])
                app_part[:, 1] -= n
                l_nodes += app_part.tolist()
            l_nodes = np.hstack((l_nodes,
                                 np.zeros((n ** 2, 1))))
            return l_nodes

        ps = EquationSolver.solve_equation([2, 15, 13, 6 * (1 - n_nodes)])
        if len(ps) > 1:
            raise NotImplementedError(f"Revise the case for N={n_nodes}: ps = {ps}")
        p = int(np.ceil(ps[0]))
        nodes = []
        # nodes_left = n_nodes
        layer_nodes = build_layer(p)
        layer_nodes[:, :2] -= np.full(layer_nodes[:, :2].shape, layer_nodes[:, :2].mean(axis=0))

        nodes += layer_nodes.tolist()
        nodes_left = n_nodes - len(nodes)

        for i in range(p + 1):
            if nodes_left == 0:
                break
            s = p + 1 - i
            layer_nodes = build_layer(s)
            layer_nodes[:, 2] = -1 - i
            layer_nodes[:, :2] -= np.full(layer_nodes[:, :2].shape, layer_nodes[:, :2].mean(axis=0))
            layer_nodes = layer_nodes[:nodes_left]
            nodes += layer_nodes.tolist()
            nodes_left -= layer_nodes.shape[0]

        return np.array(nodes)
