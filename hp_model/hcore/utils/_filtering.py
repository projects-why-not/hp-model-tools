from math import np


class CoreFilterer:
    def __init__(self, lattice):
        self.lattice = lattice

    def check_compatibility(self, hcore, hp_sequence):
        # result = {"verdict": False}
        def buried_to_surface_dfs(cur_path, v_neigh_nums, max_len):
            if v_neigh_nums[cur_path[-1]] < len(self.lattice.vectors):
                return len(cur_path)
            if len(cur_path) == max_len:
                return None
            min_l = max_len
            for neigh_dir in self.lattice.vectors:
                next_node_ind = hcore.vertices.tolist().index((hcore.vertices[cur_path[-1]] + neigh_dir).tolist())
                if next_node_ind in cur_path:
                    continue
                path = buried_to_surface_dfs(cur_path + [next_node_ind], v_neigh_nums, max_len)
                if path is not None and path < min_l:
                    min_l = path
            return min_l

        # MARK: studying buried points
        v_n_neighs = hcore.classify_nodes()
        constraints = []
        # print(v_n_neighs, np.where(v_n_neighs == len(self.lattice.vectors)))
        for p_ind in np.where(v_n_neighs == len(self.lattice.vectors))[0]:
            path_lens = []
            for neigh_dir in self.lattice.vectors:
                path_lens += [buried_to_surface_dfs([p_ind,
                                                     hcore.vertices.tolist().index((hcore.vertices[p_ind] + neigh_dir).tolist())],
                                                    v_n_neighs,
                                                    len(hp_sequence))]
            path_lens = sorted(path_lens)
            if ("H" * (path_lens[0] + path_lens[1] - 1) not in hp_sequence) \
                    and (hp_sequence.find("H" * path_lens[0]) != 0) \
                    and (hp_sequence[::-1].find("H" * path_lens[0]) != 0):
                return {"verdict": False}

        return {"verdict": True,
                "constraints": constraints}
