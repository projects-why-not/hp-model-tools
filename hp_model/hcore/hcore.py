# coding=utf-8
from ..math import RotationMath, ContactMath
import numpy as np
# from .symmetry import SymmetryInspection


class HCore:
    def __init__(self, points, lattice):
        self.vertices = points  # .astype(int)
        self.lattice = lattice
        self.basis = lattice.vectors
        self.n_contacts = ContactMath.compute_contact_num(self.vertices,
                                                          self.basis)

    def compute_vertex_powers(self):
        neigh_nums = []
        for node in self.vertices:
            v_n = 0
            for neigh_dir in self.basis:
                # if (self.vertices == node + neigh_dir).any(axis=0):
                if (node + neigh_dir).tolist() in self.vertices.tolist():
                    v_n += 1
            neigh_nums += [v_n]
        return np.array(neigh_nums)

    def compute_buriedness(self):
        def update_dists(M, node_neigh_nums, cur_p_ind):
            next_p0_inds = []
            for v in self.basis:
                next_node = self.vertices[cur_p_ind] + v
                next_node_ind = np.where(np.abs(self.vertices - np.full(self.vertices.shape, [next_node])).sum(axis=1) < 1e-10)[0]

                # if next_node.tolist() not in self.vertices.tolist():
                if next_node_ind.shape[0] == 0:
                    continue
                next_node_ind = next_node_ind[0]  # self.vertices.tolist().index(next_node.tolist())
                if node_neigh_nums[next_node_ind] < len(self.basis):
                    if M[next_node_ind] > 1:
                        M[next_node_ind] = 1
                        next_p0_inds += [next_node_ind]
                    # else - we have already fully studied this node
                else:
                    if M[next_node_ind] > M[cur_p_ind] + 1:
                        M[next_node_ind] = M[cur_p_ind] + 1
                        next_p0_inds += [next_node_ind]
            for ind in next_p0_inds:
                update_dists(M, node_neigh_nums, ind)

        node_n_neighs = self.compute_vertex_powers()
        p0 = np.argmin(node_n_neighs)
        depths = np.full(len(self.vertices), np.inf)
        depths[p0] = 1
        update_dists(depths, node_n_neighs, p0)

        return depths

    def check_connectivity(self):
        def wfs(cur_node_ind, visited_nodes):
            if cur_node_ind in visited_nodes:
                return
            next_nodes = []
            for v in self.basis:
                new_node = self.vertices[cur_node_ind] + v
                if new_node.tolist() in self.vertices.tolist():
                    new_node_ind = self.vertices.tolist().index(new_node.tolist())
                    if new_node_ind in visited_nodes:
                        continue
                    next_nodes += [new_node_ind]
            visited_nodes += [cur_node_ind]
            for ind in next_nodes:
                wfs(ind, visited_nodes)

        wfs_vertices = []
        wfs(0, wfs_vertices)
        # print(wfs_vertices)
        return len(wfs_vertices) == len(self.vertices)

    def dec_n_contacts(self):
        tgt_n_contacts = self.n_contacts - 1
        if tgt_n_contacts < len(self.vertices) - 1:
            raise Exception("Minimal number of contacts reached!")
        new_core = None
        for i, node in enumerate(self.vertices):
            new_nodes = np.delete(np.array(self.vertices), i, axis=0)
            for v in self.basis:
                new_core = HCore(np.vstack((new_nodes, [node + v])), self.lattice)
                if new_core.n_contacts == tgt_n_contacts and new_core.check_connectivity():
                    return new_core
        raise Exception("Failed to decrease n_contacts by 1!")

    def dec_n_contacts_no_ends(self, select_max_degree=True):
        tgt_n_contacts = self.n_contacts - 1
        if tgt_n_contacts < len(self.vertices) - 1:
            raise Exception("Minimal number of contacts reached!")
        new_core = None

        vert_powers = self.compute_vertex_powers()
        vert_depths = self.compute_buriedness()
        vert_rm_candidates = []
        for i, node in enumerate(self.vertices):
            if vert_depths[i] > 1:
                continue
            if vert_powers[i] < 3:
                continue
            is_candidate = True
            for j, neigh_node in enumerate(self.vertices):
                if ((node - neigh_node).tolist() in self.basis) and (vert_powers[j] < 3):
                    is_candidate = False
                    break
            if is_candidate:
                vert_rm_candidates += [i]

        if len(vert_rm_candidates) == 0:
            raise Exception("No way to create next structure without end vertices!")

        vert_rm_candidates.sort(key=lambda x: (-1 if select_max_degree else 1) * vert_powers[x])
        for i in vert_rm_candidates:
            # i = vert_rm_candidates[0]
            new_nodes = np.delete(np.array(self.vertices), i, axis=0)
            for k, node in enumerate(new_nodes):
                for v in self.basis:
                    new_core = HCore(np.vstack((new_nodes, [node + v])), self.lattice)
                    if new_core.n_contacts == tgt_n_contacts and new_core.check_connectivity():
                        return new_core

        raise Exception("Failed to decrease n_contacts by 1!")

    @classmethod
    def from_step_sequence(cls, sequence, lattice):
        raise Exception("TODO: make step sequence converter!")
