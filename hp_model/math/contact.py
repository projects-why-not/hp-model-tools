import numpy as np


class ContactMath:
    def __init__(self):
        pass

    @classmethod
    def compute_contact_num(cls, vertices, basis):
        n = 0
        for i in range(len(vertices) - 1):
            for j in range(i + 1, len(vertices)):
                if np.round(vertices[i] - vertices[j]).tolist() in basis:
                    n += 1
        return n
