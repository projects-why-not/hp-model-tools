from .hcore import HCore
from ..math import np
from .utils.symmetry import SymmetryInspection


class HCoreGenerator:
    def __init__(self, lattice, N, K, N_in=-1):
        self.__lattice = lattice
        self.__N = N
        self.__K = K
        self.__N_in = N_in
        self._cores = []
        self._next_return_core_ind = -1
        self._next_core_generate_ind = -1

    def __iter__(self):
        core0 = HCore(self.__lattice.construct_max_dense_core(self.__N), self.__lattice)
        while core0.n_contacts > self.__K:
            core0 = core0.dec_n_contacts_no_ends()
        self._cores = [core0]
        self._next_return_core_ind = 0
        self._next_core_generate_ind = 0
        return self

    def __next__(self):
        def select_exclude_candidates(core_nodes, prev_nodes, length):
            if len(prev_nodes) == length:
                return prev_nodes.reshape((-1,2))
            neighs = [[prev_nodes[-1] + bas_v] for bas_v in self.__lattice.vectors
                      if np.where((core_nodes == prev_nodes[-1] + bas_v).all(axis=1))[0].shape[0] > 0
                      and np.where((prev_nodes == prev_nodes[-1] + bas_v).all(axis=1))[0].shape[0] == 0]
            candidates = [select_exclude_candidates(core_nodes,
                                                    np.vstack((prev_nodes, neigh)),
                                                    length).reshape((-1, length, 2))
                          for neigh in neighs]
            if len(candidates) == 0:
                return np.array([])
            elif len(candidates) == 1:
                return candidates[0]
            candidates = np.vstack(candidates)
            return candidates

        def get_analogues(core, cur_analogues):     # , ind_tree, cur_ind):
            analogues = []
            # n_added = 0
            vert_set = set(["_".join(v.astype(str).tolist()) for v in core.vertices])
            for node in core.vertices:
                candidates = select_exclude_candidates(core.vertices,
                                                       np.array([node.tolist()]),
                                                       1)
                # print(candidates)
                for node_grp in candidates:
                    node_grp = np.array(node_grp).reshape((-1, 2))
                    node_grp_vert_set = set(["_".join(v.astype(str).tolist()) for v in node_grp])
                    verts_no_node_grp = [v for v in vert_set if v not in node_grp_vert_set]
                    verts_no_node_grp = np.array([list(map(float, v.split("_"))) for v in verts_no_node_grp])
                    grp_init = np.array(node_grp)
                    grp_init[:, 0] -= np.min(grp_init[:, 0])
                    grp_init[:, 1] -= np.min(grp_init[:, 1])

                    for vert in verts_no_node_grp:
                        for vect in self.__lattice.vectors:
                            new_nodes = grp_init + vect + vert
                            did_break = False
                            for new_node in new_nodes:
                                if new_node.tolist() in verts_no_node_grp.tolist():
                                    did_break = True
                                    break
                            if did_break:
                                continue

                            # print("\tREACHED NEW CONF")
                            new_conf = HCore(np.vstack((verts_no_node_grp, new_nodes)), self.__lattice)
                            if not new_conf.check_connectivity():
                                continue
                            # print("\t\tCHECKED CONNECTIVITY")
                            new_conf_powers = new_conf.classify_nodes().tolist()
                            if 1 in new_conf_powers:
                                continue
                            if max(new_conf_powers) != max(self._cores[0].classify_nodes()) or \
                                    new_conf_powers.count(max(new_conf_powers)) != self._cores[0].classify_nodes().tolist().count(max(new_conf_powers)):
                                continue

                            new_conf_cont_num = new_conf.n_contacts
                            if new_conf_cont_num == self.__K:
                                did_break = False
                                for analogue in cur_analogues + analogues:
                                    if SymmetryInspection.are_equal(analogue,
                                                                    new_conf,
                                                                    self.__lattice.rotation_angles):
                                        did_break = True
                                        break
                                if did_break:
                                    continue
                                # analogues += [new_conf]
                                analogues.append(new_conf)
                                # if cur_ind not in ind_tree.keys():
                                #     ind_tree[cur_ind] = []
                                # ind_tree[cur_ind] = ind_tree[cur_ind] + [len(cur_analogues) - 1]
                                # n_added += 1
                                # print("ADDED")

                                # get_analogues(new_conf.vertices,
                                #                                cur_analogues, #  + analogues,
                                #                                ind_tree,
                                #                                len(cur_analogues) - 1)
                                # analogues += next_analogues
            return analogues

        if self._next_core_generate_ind >= len(self._cores):
            raise StopIteration

        if self._next_return_core_ind < len(self._cores):
            self._next_return_core_ind += 1
            return self._cores[self._next_return_core_ind - 1]

        self._cores += get_analogues(self._cores[self._next_core_generate_ind],
                                     self._cores)
        self._next_core_generate_ind += 1
        if self._next_return_core_ind < len(self._cores):
            self._next_return_core_ind += 1
            return self._cores[self._next_return_core_ind - 1]
        else:
            raise StopIteration

    # def enumerate_cores(self, N, K, N_in=-1):

    # def generate_max_dense(self, n):
    #     if self.__lattice.kind == "square":
    #         a = np.ceil(np.sqrt(n)).astype(int)
    #         b = np.floor(n / a).astype(int)
    #         core_points = np.zeros((num_h_monomers, 2))
    #         for i in range(a):
    #             core_points[b * i:b * (i + 1), 1] = i
    #             core_points[b * i:b * (i + 1), 0] = np.arange(b)[:min(b, num_h_monomers - b * i)]
    #         return HCore(core_points, self.__lattice)
    #     elif self.__lattice.kind == "triangular":
    #         pass
    #
    #     return None
    #
    # def generate_by_seq_parts(self, hp_sequence, max_p_split_num):
    #     pass

