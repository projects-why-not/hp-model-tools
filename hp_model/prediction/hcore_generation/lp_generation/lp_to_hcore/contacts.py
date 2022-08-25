
def get_n_contacts_bw_layers(k, m, x0, lattice_kind="triangular"):
    """
    :param k: number of contacts in 1st row
    :param m: number of contacts in 2nd row
    :param x0: offset of 2nd layer from 1st
    :param lattice_kind: type of lattice. Can be either 'triangular' or 'quadratic'
    :return: number of contacts
    """
    if lattice_kind not in ["quadratic", "triangular"]:
        raise ValueError("lattice_kind must be either 'triangular' or 'quadratic'!")

    def get_layer_overlay(k, m, x0):
        if x0 <= -m or x0 >= k + 1:
            return 0
        lb, rb = max(0, x0), min(k, x0 + m)
        return rb - lb

    if x0 <= -m or x0 >= k + 1:
        return 0
    overlay_n = get_layer_overlay(k, m, x0)
    w_kmx0 = 0
    if overlay_n == 0:
        w_kmx0 = 0 # + (x0 + m == k)

    if lattice_kind == "quadratic":
        return overlay_n

    w_kmx0 = 2 * overlay_n # - (x0 <= 0) # + (x0 + m > k and x0 <= k)

    #     w_kmx0 += (k >= x0 and k < x0 + m)
    w_kmx0 -= (x0 <= 0)
    w_kmx0 += (x0 + m > k and x0 <= k)
    return w_kmx0
