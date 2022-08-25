from constraint import *
from ....math import np


class ConstraintBuilder:
    def __init__(self, hp_sequence, hcore, lattice):
        self._hcore_nodes = hcore.vertices.astype(int)
        self._hcore_node_depths = hcore.compute_buriedness()
        self._hp_sequence = hp_sequence
        self._basis = lattice.vectors

    def get_constraints(self):
        seq_analysis = self.check_feasibility(verbose=True)
        if not seq_analysis["verdict"]:
            raise Exception("H-Core and HP-sequence do not fit!")

    def check_feasibility(self, verbose=False):
        # TODO: distance to prev also!
        hp_dists = {"H-P": [],
                    "P-H": []}

        for i in range(len(self._hp_sequence)):
            if self._hp_sequence[i] == "P":
                tgt_monomer = "H"
                key = "P-H"
            else:
                tgt_monomer = "P"
                key = "H-P"

            d_to_prev = self._hp_sequence[:i + 1][::-1].find(tgt_monomer)
            d_to_next = self._hp_sequence[i:].find(tgt_monomer)
            d_to_prev = d_to_prev if d_to_prev > 0 else len(self._hp_sequence)
            d_to_next = d_to_next if d_to_next > 0 else len(self._hp_sequence)
            d_to_nearest = min(d_to_next,
                               d_to_prev)
            hp_dists[key] = hp_dists[key] + [d_to_nearest]

        if max(self._hcore_node_depths) > max(hp_dists["H-P"]):
            if verbose:
                return {"verdict": False,
                        "max_core_depth": max(self._hcore_node_depths),
                        "max_H-P_dist": max(hp_dists["H-P"]),
                        "max_P-H_dist": max(hp_dists["P-H"])}
            return False

        # TODO: other filtration criteria?

        if verbose:
            return {"verdict": True,
                    "max_core_depth": max(self._hcore_node_depths),
                    "max_H-P_dist": max(hp_dists["H-P"]),
                    "max_P-H_dist": max(hp_dists["P-H"])}
        return True


class HCoreLinearCSP:
    def __init__(self, sequence, hcore, lattice):
        self._sequence = sequence
        self._core_nodes = hcore.vertices.tolist()
        self._core_node_depths = hcore.compute_buriedness()
        self._core_nodes_str_set = set([f"{node[0]}_{node[1]}" for node in self._core_nodes])
        if len(self._core_nodes) != len([m for m in self._sequence if m == "H"]):
            raise ValueError("H-Core Power and H number in sequence do not match!")
        self._lat_basis = lattice.vectors
        self._seq_analysis = ConstraintBuilder(sequence, hcore, lattice).check_feasibility(True)
        self._problem = self._make_problem()

    def _make_problem(self):
        problem = Problem()

        # all node enumeration
        all_nodes = list(self._core_nodes)
        p_node_dists = []
        for i, node in enumerate(self._core_nodes):
            if self._core_node_depths[i] > 1:
                continue
            for v in self._lat_basis:
                for dist in range(1, self._seq_analysis["max_P-H_dist"] + 1):
                    new_node = (np.array(node) + dist * v).tolist()
                    if new_node not in all_nodes:
                        all_nodes += [new_node]
                        p_node_dists += [dist]

        # variable creation
        # TODO: avoid infeasible combinations
        for i, monomer in enumerate(self._sequence):
            problem.addVariable(f"x_{i}", [-len(self._sequence), len(self._sequence)])
            problem.addVariable(f"y_{i}", [-len(self._sequence), len(self._sequence)])
            for j, node_coord in enumerate(all_nodes):
                problem.addVariable(f"s_{i}_{j}",
                                    [0, 1])

                if (monomer == "H" and j >= len(self._core_nodes)) or (monomer == "P" and j < len(self._core_nodes)):
                    problem.addConstraint(lambda a: a == 0,
                                          (f"s_{i}_{j}"))
            problem.addConstraint(lambda x, *a: x - sum([all_nodes[j][0] * a[j]]) == 0,
                                  tuple([f"x_{i}"] + [f"s_{i}_{j}" for j in range(len(all_nodes))]))
            problem.addConstraint(lambda y, *a: y - sum([all_nodes[j][0] * a[j]]) == 0,
                                  tuple([f"y_{i}"] + [f"s_{i}_{j}" for j in range(len(all_nodes))]))

        for i in range(len(self._sequence)):
            problem.addConstraint(lambda *a: sum(a) == 1,
                                  (f"s_{i}_{j}" for j in range(len(all_nodes))))

        for j in range(len(all_nodes)):
            problem.addConstraint(lambda *a: sum(a) <= 1,
                                  (f"s_{i}_{j}" for i in range(len(self._sequence))))

        for i in range(1, len(self._sequence)):
            problem.addConstraint(lambda xi, yi, xim1, yim1: -1 <= (xi - xim1) + (yi - xim1) <= 1,
                                  (f"x_{i}", f"y_{i}", f"x_{i - 1}", f"y_{i - 1}"))
            problem.addConstraint(lambda xi, yi, xim1, yim1: -1 <= (xim1 - xi) + (yim1 - xi) <= 1,
                                  (f"x_{i}", f"y_{i}", f"x_{i - 1}", f"y_{i - 1}"))

        return problem

    def solve(self, ax=None, return_all=False):
        # print(f"entered solve. problem = {self._problem}")
        if return_all:
            return self._problem.getSolutions()

        solution = self._problem.getSolution()
        if solution is None:
            return {}
        if ax is not None:
            points = np.array([solution[f"L_{i}"] for i in range(len(self._sequence))])
            ax.scatter(*points[np.array(list(self._sequence)) != "H"].T,
                       # c=["red" if self._sequence[i] == "H" else "green" for i in range(len(self._sequence))],
                       facecolors="white", edgecolors="black",
                       s=100,
                       zorder=1)
            ax.scatter(*points[np.array(list(self._sequence)) == "H"].T,
                       c="black", # ["red" if self._sequence[i] == "H" else "green" for i in range(len(self._sequence))],
                       # facecolors="none", edgecolors="black",
                       s=100,
                       zorder=1)
            ax.plot(*points.T, c="black", linewidth=1, zorder=0)

        return solution
