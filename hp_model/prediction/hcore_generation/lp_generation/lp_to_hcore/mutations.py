from .contacts import get_n_contacts_bw_layers
import numpy as np
from itertools import product


def get_offsets(k, m, w, lattice_kind="triangular"):
    if k == 0 or m == 0:
        # may be any, but let us consider 0 only
        return [0]
        # return np.arange(-max(k,m), max(k,m) + 1)
    ds = []
    for d in range(-m + 1, k + 1):
        if get_n_contacts_bw_layers(k, m, d, lattice_kind) == w:
            ds += [d]
    return ds


def get_core_mutations(n_row_pts, k_bw_rows):
    def get_possible_offsets(constraints):
        offsets = set(get_offsets(*constraints[0]))
        for constr in constraints[1:]:
            offsets = offsets & set(get_offsets(*constr))
        return list(offsets)

    if type(list(k_bw_rows.keys())[0][0]) == int:
        # MARK: if 2D lattice
        offsets = [get_possible_offsets([[n_row_pts[r]
                                          for r in pair][::-1] + [K_pair]])
                   for pair, K_pair in k_bw_rows.items()]
        offsets = [row_off if n_row_pts[i + 1] > 0 and n_row_pts[i] > 0 else [0] for i, row_off in enumerate(offsets)]
        mutations = np.array(np.meshgrid(*offsets)).T.reshape(-1, len(n_row_pts) - 1, 1)

    else:
        # MARK: if 3D lattice
        # MARK: offset of layer (i) from (i-1) by Oy. Offset of rows by Ox
        n_layers, n_rows_in_layer = [np.sqrt(len(n_row_pts)).astype(int)] * 2
        offsets = []
        matr_pts_in_row = np.array([n_row_pts[(i, j)]
                                    for i in range(n_layers)
                                    for j in range(n_rows_in_layer)]).reshape((n_layers, n_rows_in_layer))
        print("M:\n", matr_pts_in_row)

        # first_layer_passed = False
        empty_layers_passed = False
        for i in range(n_layers):
            nonempty_rows = np.where(matr_pts_in_row[i] > 0)[0]
            if len(nonempty_rows) == 0:
                if empty_layers_passed:
                    offsets += [[[0,
                                  np.full((len(list(range(i % 2, n_rows_in_layer, 2))),
                                           1),
                                          0).tolist()]]]
                continue
            empty_layers_passed = True

            # MARK: layer kind. May be 0 or 1, defines odevity of layer row indices
            layer_kind = nonempty_rows[0] % 2
            # MARK: flag whether the current layer is the first
            is_first_layer = len(offsets) == 0

            layer_offsets = []

            # MARK: if first layer, only dy = 0 is considered. Otherwise - [-n_rows, r_rows]
            if is_first_layer:
                dy_range = [0]
                # first_layer_passed = True
            else:
                dy_range = range(-n_rows_in_layer * 2,
                                 n_rows_in_layer * 2 + 1,
                                 2)

            # TODO: find all combinations of horizontal offsets first!
            # layer_horiz_offsets = []
            # for j in range(layer_kind, n_rows_in_layer, 2):
            #     # MARK: if first layer, only offset 0 is allowed for first row
            #     if j == layer_kind and is_first_layer:
            #         layer_horiz_offsets += [[0]]
            #         continue
            #
            #     if ((i, j), (i, j - 2)) in k_bw_rows.keys():
            #         constr = [n_row_pts[(i,j)],
            #                   n_row_pts[(i, j - 2)],
            #                   k_bw_rows[((i, j), (i, j - 2))],
            #                   "quadratic"]
            #         row_horiz_offsets = get_possible_offsets([constr])
            #         layer_horiz_offsets += [row_horiz_offsets]
            #     else:
            #         layer_horiz_offsets += [np.arange(-n_row_pts[(i,j)],
            #                                           n_row_pts[(i,j)] + 1)]
            #
            # print(matr_pts_in_row[i])
            # print(layer_horiz_offsets)
            #
            # if is_first_layer:
            #     offsets += [[0, layer_horiz_offsets]]
            #     continue
            #
            # # n_layer_inter_contacts = sum([sum([v for k,v in k_bw_rows.items()
            # #                                    if (i,j) in k and ((k[0][0] == i and k[1][0] == i - 1) or (k[1][0] == i and k[0][0] == i - 1))])
            # #                               for j in range(layer_kind, n_rows_in_layer, 2)])
            #
            # n_layers_inter_contacts = sum([v for k, v in k_bw_rows.items()
            #                                if ((k[0][0] == i and k[1][0] == i - 1) or (k[1][0] == i and k[0][0] == i - 1))])
            #
            # print(n_layers_inter_contacts)
            #
            # for dy in dy_range:
            #     n_dy_inter_contacts = 0
            #



                # constr = []
                # for pair, K_pair in k_bw_rows.items():
                #     if (i, j) not in pair:
                #         continue
                #     constr = []
                #     other_row_pair_ind = 1 - list(pair).index((i, j))
                #     other_row = pair[other_row_pair_ind]
                #     if other_row[0] != i:
                #         # MARK: if other row lies in another layer
                #         continue
                #
                #     if other_row[1] > j:
                #         # MARK: if other row lies behind the current
                #         continue
                #
                #     lat_type = "quadratic"
                #     constr += [n_row_pts[r] for r in pair]
                #     constr += [K_pair, lat_type]
                #     break




            # TODO: then, for each combination of offsets (dy, horizontal), check if it is feasible by prev. layer!


            #
            # for dy in dy_range:
            #     row_offsets = []
            #     dy_invalid = False
            #
            #     for j in range(layer_kind, n_rows_in_layer, 2):
            #         # MARK: if first layer, only offset 0 is allowed for first row
            #         if j == layer_kind and is_first_layer:
            #             row_offsets += [[0]]
            #             continue
            #
            #         constraints = []
            #         for pair, K_pair in k_bw_rows.items():
            #             if (i, j) not in pair:
            #                 continue
            #             constr = []
            #             other_row_pair_ind = 1 - list(pair).index((i, j))
            #             other_row = pair[other_row_pair_ind]
            #             if other_row[0] > i:
            #                 # MARK: if other row lies above the current
            #                 continue
            #
            #             if other_row[0] == i:
            #                 # MARK: if other row is in current layer
            #                 if other_row[1] > j:
            #                     # MARK: if other row lies behind the current
            #                     continue
            #
            #                 lat_type = "quadratic"
            #                 constr += [n_row_pts[r] for r in pair]
            #             else:
            #                 # MARK: if other row is in previous layer
            #                 lat_type = "triangular"
            #                 other_row = (other_row[0], other_row[1] + dy)
            #                 if other_row[1] < 0 or other_row[1] >= n_rows_in_layer:
            #                     # continue
            #                     dy_invalid = True
            #                     break
            #                 constr = list(pair)
            #                 # constr += [(i, j)] * 2
            #                 # constr[other_row_pair_ind] = other_row
            #                 constr = [n_row_pts[r] for r in constr]
            #             constr += [K_pair, lat_type]
            #             constraints += [constr]
            #         if len(constraints) == 0 or dy_invalid:
            #             dy_invalid = True
            #             break
            #         off_row = get_possible_offsets(constraints)
            #         if len(off_row) > 0:
            #             row_offsets += [off_row]
            #         else:
            #             dy_invalid = True
            #             break
            #     if dy_invalid:
            #         continue
            #     layer_offsets += [[dy, row_offsets]]
            #
            # if len(layer_offsets) > 0:
            #     offsets += [layer_offsets]
            # else:
            #     offsets += [[[0,
            #                   np.full((len(list(range(i % 2, n_rows_in_layer, 2))),
            #                            1),
            #                           0).tolist()]]]



        # raise NotImplementedError("Review mutation generation!")

        for i in range(1, n_layers):
            if np.sum(matr_pts_in_row[i - 1]) == 0 or np.sum(matr_pts_in_row[i]) == 0:
                offsets += [[[0,
                              np.full((len(list(range(i % 2, n_rows_in_layer, 2))),
                                       1),
                                      0).tolist()]]]
                continue

            layer_offsets = []
            for dy in range(-n_rows_in_layer * 2,
                            n_rows_in_layer * 2 + 1,
                            2):
                row_offsets = []
                dy_invalid = False
                for j in range(i % 2, n_rows_in_layer, 2):
                    constraints = []
                    for pair, K_pair in k_bw_rows.items():
                        constr = []
                        if (i, j) not in pair:
                            continue
                        other_row_pair_ind = 1 - list(pair).index((i, j))
                        other_row = pair[other_row_pair_ind]
                        if other_row[0] > i:
                            continue
                        lat_type = "triangular"
                        if other_row[0] == i:
                            # same layer
                            constr += [n_row_pts[r] for r in pair]
                            lat_type = "quadratic"
                        else:
                            other_row = (other_row[0], other_row[1] + dy)
                            if other_row[1] < 0 or other_row[1] >= n_rows_in_layer:
                                # continue
                                dy_invalid = True
                                break
                            constr += [(i, j)] * 2
                            constr[other_row_pair_ind] = other_row
                            constr = [n_row_pts[r] for r in constr]
                        constr += [K_pair, lat_type]
                        constraints += [constr]
                    if len(constraints) == 0 or dy_invalid:
                        dy_invalid = True
                        break

                    off_row = get_possible_offsets(constraints)
                    if len(off_row) > 0:
                        row_offsets += [off_row]
                    else:
                        dy_invalid = True
                        break
                if dy_invalid:
                    continue

                layer_offsets += [[dy, row_offsets]]
            if len(layer_offsets) > 0:
                offsets += [layer_offsets]
            else:
                offsets += [[[0,
                              np.full((len(list(range(i % 2, n_rows_in_layer, 2))),
                                       1),
                                      0).tolist()]]]

        print("\n\n\n****** offsets to mutations:")
        print("offsets:\n", offsets)


        mutations = [[[off_comb[0], np.array(np.meshgrid(*off_comb[1])).T.reshape((-1, len(off_comb[1])))]
                      for off_comb in layer]
                     for layer in offsets]
        # print("mutations init:\n", mutations)
        mutations = [[(off_comb[0], row_off_comb)
                      for off_comb in layer
                      for row_off_comb in off_comb[1]]
                     for layer in mutations]
        # print("mutations 2:\n", mutations)
        mutations = list(product(*mutations))
        # print("mutations final:\n", mutations)

        # input("goon?")

    return mutations
