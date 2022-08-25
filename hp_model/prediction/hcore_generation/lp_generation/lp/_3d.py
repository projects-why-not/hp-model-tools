from pulp import *
import numpy as np


def make_3D_FCC_core_problem(N, Kmax=None):
    M = 1e3
    #     P = N
    P = np.ceil(N ** (1 / 3)).astype(int) + 1  # TODO: maybe increase?

    problem = LpProblem("basic_3D_construction",
                        LpMaximize)
    #     n_l = N # np.ceil(N / 2).astype(int)

    # n points in layer (i), row (j)
    xs = LpVariable.dicts("x",
                          ((i, j) for i in range(P) for j in range(P)),
                          lowBound=0,
                          cat="Integer")
    # n contacts in layer (i) between rows (j), (j - 1)
    ws = LpVariable.dicts("w",
                          ((i, j, k) for i in range(P) for j in range(P) for k in range(2)),
                          lowBound=0,
                          cat="Integer")
    # z_{i,j} = 1 if row (j) of layer (i) is non-empty, otherwise 0
    zs = LpVariable.dicts("z",
                          ((i, j) for i in range(P) for j in range(P)),
                          cat="Binary")
    # v_{i,j} = min(x_{i,j}, x_{i,j-1})
    vs = LpVariable.dicts("v",
                          ((i, j, k) for i in range(P) for j in range(P) for k in range(2)),
                          lowBound=0,
                          cat="Integer")
    # p_{i,j} = 1 if x_{i,j} == x_{i, j - 1}
    ps = LpVariable.dicts("p",
                          ((i, j, k) for i in range(P) for j in range(P) for k in range(2)),
                          # lowBound=0,
                          cat="Binary")
    # t_{i,j} = 1 if x_{i,j} == x_{i, j - 1} == 0 (??? TODO: CHECK!)
    ts = LpVariable.dicts("t",
                          ((i, j, k) for i in range(P) for j in range(P) for k in range(2)),
                          # lowBound=0,
                          cat="Binary")
    qs = LpVariable.dicts("q",
                          ((i, j) for i in range(P) for j in range(2, P)),
                          # lowBound=0,
                          cat="Integer")

    inter_layer = LpVariable("bl_n",
                             lowBound=0,
                             cat="Integer")
    intra_layer_intra_row = LpVariable("ilir_n",
                                       lowBound=0,
                                       cat="Integer")
    intra_layer_inter_row = LpVariable("ilbr_n",
                                       lowBound=0,
                                       cat="Integer")
    K = LpVariable("K",
                   lowBound=0,
                   cat="Integer")

    # consistency
    problem += lpSum([xs[(i, j)] for i in range(P) for j in range(P)]) - N == 0

    # inter-layer contacts
    for i in range(1, P):
        for j in range(i % 2, P):
            for direct in [-1, 1]:
                k = max(0, direct)
                if j + direct < 0 or j + direct >= P:
                    #                     problem += ws[(i,j,k)] == 0
                    continue
                # min(x_i, x_{i - 1}) <-> v_i
                problem += vs[(i, j, k)] - xs[(i, j)] <= 0
                problem += vs[(i, j, k)] - xs[(i - 1, j + direct)] <= 0

                # {1 if x_i == x_{i - 1}}
                problem += xs[(i, j)] + xs[(i - 1, j + direct)] - 2 * vs[(i, j, k)] - ps[(i, j, k)] >= 0

                # if x_i == x_{i - 1} == 0, prev. fails. so:
                problem += 2 * (xs[(i, j)] + xs[(i - 1, j + direct)]) - 2 * vs[(i, j, k)] - M * ts[(i, j, k)] <= 0

                # w(x_i, x_{i - 1})
                problem += ws[(i, j, k)] == 2 * vs[(i, j, k)] - 1 + ps[(i, j, k)] + 1 - ts[(i, j, k)]
    problem += inter_layer == lpSum([ws[(i, j, k)]
                                     for i in range(1, P)
                                     for j in range(i % 2, P)
                                     for k in range(2)
                                     if not (j + (-1 + 2 * k) < 0 or j + (-1 + 2 * k) >= P)])

    # TARGET: intra-layer, inter-row
    for i in range(P):
        for j in range(2, P):
            problem += qs[(i, j)] <= xs[(i, j)]
            problem += qs[(i, j)] <= xs[(i, j - 2)]
    problem += intra_layer_inter_row == lpSum([qs[(i, j)] for i in range(P) for j in range(2, P)])
    #     problem += intra_layer_inter_row == lpSum([ws[(i,j)] for i in range(P) for j in range(1, P)])

    # TARGET: intra-layer, intra-row
    # z_i = {1 if x_i >= 1; 0 if x_i == 0}
    for i in range(P):
        for j in range(P):
            problem += M * zs[(i, j)] >= xs[(i, j)]  # zs[i] == (n_l - i) * xs[i]
    problem += intra_layer_intra_row == N - lpSum([zs[(i, j)] for i in range(P) for j in range(P)])

    # TARGET
    problem += K == intra_layer_inter_row + intra_layer_intra_row + inter_layer

    if Kmax is not None:
        problem += K <= Kmax

    problem += K

    return problem
