from ..lattice import Lattice, np


class SquareLattice(Lattice):
    def __init__(self):
        Lattice.__init__(self, np.array([[0, 1],
                                         [0, -1],
                                         [1, 0],
                                         [-1, 0]]), [0, np.pi/2, -np.pi/2, np.pi],
                         "square")
