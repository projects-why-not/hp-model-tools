from pulp import *


def make_2D_FCC_core_problem(N, Kmax=None):
    M = 1e3

    problem = LpProblem("basic_2D_construction",
                        LpMaximize)
    n_l = N # np.ceil(N / 2).astype(int)

    xs = LpVariable.dicts("x",
                          (i for i in range(n_l)),
                          lowBound=0,
                          cat="Integer")
    ws = LpVariable.dicts("w",
                          (i for i in range(1, n_l)),
                          lowBound=0,
                          cat="Integer")
    zs = LpVariable.dicts("z",
                          (i for i in range(n_l)),
                          cat="Binary")
    vs = LpVariable.dicts("v",
                          (i for i in range(1, n_l)),
                          lowBound=0,
                          cat="Integer")
    ps = LpVariable.dicts("p",
                          (i for i in range(1, n_l)),
                          # lowBound=0,
                          cat="Binary")
    ts = LpVariable.dicts("t",
                          (i for i in range(1, n_l)),
                          # lowBound=0,
                          cat="Binary")

    inter = LpVariable("inter_n",
                       lowBound=0,
                       cat="Integer")
    intra = LpVariable("intra_n",
                       lowBound=0,
                       cat="Integer")
    K = LpVariable("K",
                   lowBound=0,
                   cat="Integer")


    # consistency
    problem += lpSum([xs[i] for i in range(n_l)]) - N == 0

    for i in range(1, n_l):
        # min(x_i, x_{i - 1}) <-> v_i
        problem += vs[i] - xs[i] <= 0
        problem += vs[i] - xs[i - 1] <= 0

        # {1 if x_i == x_{i - 1}}
        problem += xs[i] + xs[i - 1] - 2 * vs[i] - ps[i] >= 0

        # if x_i == x_{i - 1} == 0, prev. fails. so:
        problem += 2 * (xs[i] + xs[i - 1]) - 2 * vs[i] - M * ts[i] <= 0

        # w(x_i, x_{i - 1})
        problem += ws[i] == 2 * vs[i] - 1 + ps[i] + 1 - ts[i]

    # inter-contacts
    problem += inter == lpSum([ws[i] for i in range(1, n_l)])

    # z_i = {1 if x_i >= 1; 0 if x_i == 0}

    for i in range(n_l):
        problem += M * zs[i] >= xs[i] # zs[i] == (n_l - i) * xs[i]

    # intra-contacts
    problem += intra == N - lpSum([zs[i] for i in range(n_l)])

    # target
    problem += K == inter + intra

    if Kmax is not None:
        problem += K <= Kmax

    problem += K

    return problem
