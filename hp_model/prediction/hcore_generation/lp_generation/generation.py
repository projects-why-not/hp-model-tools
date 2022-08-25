from .lp import make_3D_FCC_core_problem, make_2D_FCC_core_problem
from .lp_to_hcore import lp_to_cores


def make_cores(hp_sequence, K_max, fit_3d):
    N = hp_sequence.count("H")
    if fit_3d:
        prob = make_3D_FCC_core_problem(N, K_max)
    else:
        prob = make_2D_FCC_core_problem(N, K_max)
    prob.solve()
    print("OBJECTIVE:", prob.objective.value())
    cores = lp_to_cores(prob.variablesDict())

    return cores
