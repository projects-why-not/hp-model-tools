import numpy as np


class Lattice:
    def __init__(self, vectors, rotation_angles, kind):
        self.vectors = vectors.tolist()
        self.rotation_angles = rotation_angles
        self.kind = kind

    @classmethod
    def construct_max_dense_core(cls, n_nodes):
        raise Exception("Not implemented!")
