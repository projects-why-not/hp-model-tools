from .....lattice import PseudoTriangularLattice, FCC3DLattice
from .....hcore import HCore
import numpy as np


def _generate_2DFCC(n_row_pts, offsets):
    lattice = PseudoTriangularLattice()
    points = np.empty((0,2))
    i = 0
    while n_row_pts[i] == 0:
        i += 1
    layer = np.vstack((np.arange(n_row_pts[i]), np.full(n_row_pts[i], len(n_row_pts) - i))).T
    points = np.vstack((points, layer))
    i += 1
    while i < len(n_row_pts):
        if n_row_pts[i] == 0:
            i += 1
            continue
        layer = np.vstack((np.arange(n_row_pts[i]),
                           np.full(n_row_pts[i],
                                   len(n_row_pts) - i))).T
        layer[:,0] += np.sum(offsets[:i])
        points = np.vstack((points, layer))
        i += 1
    return HCore(points, lattice)


def _generate_3DFCC(n_row_pts, offsets):
    # print("3D FCC generation started")
    # print("n_row_pts:\n", n_row_pts)
    # print("\noffsets:\n", offsets)

    lattice = FCC3DLattice()
    points = np.empty((0, 3))
    n_layers, n_rows_in_layer = [np.sqrt(len(n_row_pts)).astype(int)] * 2
    matr_pts_in_row = np.array([n_row_pts[(i, j)]
                                for i in range(n_layers)
                                for j in range(n_rows_in_layer)]).reshape((n_layers, n_rows_in_layer))

    # print("M points in row:\n", matr_pts_in_row)

    i = 0
    while matr_pts_in_row[i].sum() == 0:
        i += 1

    layer = np.vstack([np.vstack((i % 2 + np.arange(0, n_row_pts * 2, 2),
                                  np.full(n_row_pts, row_j),
                                  np.full(n_row_pts, i))).T
                       for row_j, n_row_pts in enumerate(matr_pts_in_row[i])])
    points = np.vstack((points, layer))
    i += 1

    while i < len(matr_pts_in_row):
        if matr_pts_in_row[i].sum() == 0:
            i += 1
            continue
        layer_off = offsets[i - 1]
        dy = sum([off[0] for off in offsets[:i]])
        print(f"Creating layer. off: {layer_off}; pts: {matr_pts_in_row[i]}")
        layer = np.vstack([np.vstack((i % 2 + np.arange(0, n_row_pts * 2, 2) - 2 * sum(layer_off[1][:-(i % 2) + row_j // 2 + 1]),
                                      np.full(n_row_pts, row_j) + dy, # layer_off[0],
                                      np.full(n_row_pts, i))).T
                           for row_j, n_row_pts in enumerate(matr_pts_in_row[i])
                           #                            if n_row_pts > 0
                           ])
        points = np.vstack((points, layer))
        i += 1

    # TODO: make rotation by -pi/4
    phi = -np.pi / 4
    points[:, :2] -= points[:, :2].mean(axis=0)
    M = np.array([[np.cos(phi), -np.sin(phi)],
                  [np.sin(phi), np.cos(phi)]])
    points[:, :2] = (M @ points[:, :2].T).T
    points[:, :2] *= np.sqrt(2) / 2

    # points = np.round(points)

    return HCore(points, lattice)


def generate_core(n_row_pts, offsets, is_2d=True):
    if is_2d:
        return _generate_2DFCC(n_row_pts, offsets.reshape(list(offsets.shape)[:-1]))
    else:
        return _generate_3DFCC(n_row_pts, offsets)
