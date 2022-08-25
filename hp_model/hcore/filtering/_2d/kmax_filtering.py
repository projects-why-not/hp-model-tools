

def get_Kmax(hp_sequence):
    # corollary 1 filtering
    hp_seq_HP_dist_indices = {}
    for i in range(len(hp_sequence)):
        if hp_sequence[i] == "P":
            continue
        d_to_prev_p = hp_sequence[:i + 1][::-1].find("P")
        d_to_next_p = hp_sequence[i:].find("P")
        d_to_prev_p = d_to_prev_p if d_to_prev_p > 0 else len(hp_sequence)
        d_to_next_p = d_to_next_p if d_to_next_p > 0 else len(hp_sequence)
        d_to_nearest_p = min(d_to_next_p,
                             d_to_prev_p)
        if d_to_nearest_p not in hp_seq_HP_dist_indices:
            hp_seq_HP_dist_indices[d_to_nearest_p] = []
        hp_seq_HP_dist_indices[d_to_nearest_p] = hp_seq_HP_dist_indices[d_to_nearest_p] + [i]

    min_RH = len(hp_seq_HP_dist_indices[1])
    N = hp_sequence.count("H")
    K_max = 3 * N - min_RH - 3
    return K_max
