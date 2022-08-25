from constraint import *
import numpy as np


class ConstraintBuilder:
    def __init__(self, hp_sequence, hcore, lattice):
        self._hcore_nodes = hcore.vertices  # .astype(int)
        # self._hcore_nodes -= self._hcore_nodes.min(axis=0)
        self._hcore_node_depths = hcore.compute_buriedness()
        self._hp_sequence = hp_sequence
        self._basis = lattice.vectors

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

    def get_constraints(self):
        seq_analysis = self.check_feasibility(verbose=True)
        if not seq_analysis["verdict"]:
            # return None
            raise Exception("H-Core and HP-sequence do not fit!")

        def __constraint_no_overlay(p, q):
            return p != q

        constraints = []
        var_ranges = {}

        h_nodes = {d: [] for d in np.unique(self._hcore_node_depths)}
        for i in range(len(self._hcore_nodes)):
            h_nodes[self._hcore_node_depths[i]] = h_nodes[self._hcore_node_depths[i]] + [self._hcore_nodes[i].tolist()]
        p_nodes = {d: [] for d in np.unique(self._hcore_node_depths)}
        for offset in range(1, seq_analysis["max_P-H_dist"] + 1):
            off_nodes = [(node + offset * np.array(v)).astype(int).tolist() for node in self._hcore_nodes for v in self._basis]
            off_nodes = [node for node in off_nodes if node not in (self._hcore_nodes.tolist() if offset == 1 else p_nodes[offset - 1])]
            off_nodes = [list(node) for node in list(set([tuple(node) for node in off_nodes]))]
            p_nodes[offset] = off_nodes

        for i in range(len(self._hp_sequence)):
            var_name = f"L_{i}"

            if i > 0:
                constraints += [(lambda p, q: [np.round(p[l] - q[l]) for l in range(len(p))] in self._basis,
                                 (f"L_{i}", f"L_{i - 1}"))]

            for j in range(0, i - 1):
                constraints += [(__constraint_no_overlay,  # lambda x_i, y_i, x_j, y_j: (x_i != x_j) and (y_i != y_j),
                                (f"L_{i}", f"L_{j}"))]

            if self._hp_sequence[i] == "H":
                d_to_prev_p = self._hp_sequence[:i + 1][::-1].find("P")
                d_to_next_p = self._hp_sequence[i:].find("P")
                d_to_prev_p = d_to_prev_p if d_to_prev_p > 0 else len(self._hp_sequence)
                d_to_next_p = d_to_next_p if d_to_next_p > 0 else len(self._hp_sequence)
                d_to_nearest_p = min(d_to_next_p,
                                     d_to_prev_p)
                var_range = []
                for k in range(1, d_to_nearest_p + 1):
                    if k not in h_nodes.keys():
                        break
                    var_range += h_nodes[k]
                var_ranges[var_name] = var_range

            else:
                d_to_prev_h = self._hp_sequence[:i + 1][::-1].find("H")
                d_to_next_h = self._hp_sequence[i:].find("H")
                d_to_prev_h = d_to_prev_h if d_to_prev_h > 0 else len(self._hp_sequence)
                d_to_next_h = d_to_next_h if d_to_next_h > 0 else len(self._hp_sequence)
                d_to_nearest_h = min(d_to_next_h, d_to_prev_h)
                var_range = []
                for k in range(1, d_to_nearest_h + 1):
                    var_range += p_nodes[k]
                var_ranges[var_name] = var_range

        return {"ranges": var_ranges,
                "constraints": constraints}


class HCoreCSP:
    def __init__(self, sequence, hcore, lattice):
        self._sequence = sequence
        self._core_nodes = hcore.vertices.tolist()
        self._core_node_depths = hcore.compute_buriedness()
        self._core_nodes_str_set = set([f"{node[0]}_{node[1]}" for node in self._core_nodes])

        if len(self._core_nodes) != len([m for m in self._sequence if m == "H"]):
            raise ValueError("H-Core Power and H number in sequence do not match!")
        self._lat_basis = lattice.vectors # .tolist()
        # self._init_configuration()

        self._problem = Problem()
        config = ConstraintBuilder(sequence, hcore, lattice).get_constraints()
        for k,v in config["ranges"].items():
            self._problem.addVariable(k, v)
        for constr in config["constraints"]:
            self._problem.addConstraint(constr[0],
                                        constr[1])

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

    def generate_H_constraints(self):
        node_groups_by_depth = {}
        for i in range(len(self._core_node_depths)):
            if self._core_node_depths[i] not in node_groups_by_depth:
                node_groups_by_depth[self._core_node_depths[i]] = []
            node_groups_by_depth[self._core_node_depths[i]] = node_groups_by_depth[self._core_node_depths[i]] + [self._core_nodes[i]]
