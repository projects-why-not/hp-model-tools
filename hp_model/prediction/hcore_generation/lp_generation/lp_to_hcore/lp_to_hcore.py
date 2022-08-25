from .generation import *
from .mutations import get_core_mutations


def lp_to_cores(lp_variable_dict):
    def var_to_num(a):
        if a.count("_") > 1:
            return tuple([int(v.split("_")[-1])
                          for v in a[a.find("_") + 2:-1].split(",")])
        else:
            return int(a.split("_")[-1])

    def var_to_w_key(a):
        a_parsed = var_to_num(a)
        if type(a_parsed) == tuple and len(a_parsed) == 3:
            return ((a_parsed[0],
                     a_parsed[1]),
                    (a_parsed[0] - 1,
                     a_parsed[1] - 1 + 2 * a_parsed[2]))
        elif type(a_parsed) == tuple and len(a_parsed) == 2:
            return ((a_parsed[0],
                     a_parsed[1]),
                    (a_parsed[0],
                     a_parsed[1] - 2))
        else:
            return (a_parsed,
                    a_parsed - 1)

    n_row_pts = {var_to_num(k): int(lp_variable_dict[k].value())
                 for k in sorted([k for k in lp_variable_dict
                                  if "x_" in k],
                                 key=lambda a: var_to_num(a))}
    k_bw_rows = {var_to_w_key(k): int(lp_variable_dict[k].value())
                 for k in sorted([k for k in lp_variable_dict
                                  if "w_" in k or "q_" in k],
                                 key=lambda a: var_to_num(a))}
    is_2d = type(list(k_bw_rows.keys())[0][0]) != tuple

    mutations = get_core_mutations(n_row_pts,
                                   k_bw_rows)

    cores = []
    for mut in mutations:
        cores += [generate_core(n_row_pts, mut, is_2d)]

    return cores
