from .csp import HCoreCSP


def fit_protein(hp_sequence, hcore):
    prob = HCoreCSP(hp_sequence, hcore, hcore.lattice)
    return prob.solve()
