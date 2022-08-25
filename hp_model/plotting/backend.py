from ..math import np


class Plotter:
    @classmethod
    def plot_hcore(cls, ax, hcore, dot_s=60, contact_width=1):
        vertices = hcore.vertices
        ax.scatter(*vertices.T, s=dot_s, c="black", zorder=2)
        # print(hcore.lattice.vectors)
        for i in range(len(vertices) - 1):
            for j in range(i + 1, len(vertices)):
                # print("\t", (vertices[i] - vertices[j])..tolist())

                eps = 1e-5
                node_dist = vertices[i] - vertices[j]
                vect_deltas = np.array(hcore.lattice.vectors) - np.full((len(hcore.lattice.vectors),
                                                                         len(node_dist)),
                                                                        node_dist)
                vect_deltas = np.abs(vect_deltas).sum(axis=1)
                if np.where(vect_deltas < eps)[0].shape[0] > 0:
                    # if (vertices[i] - vertices[j]).astype(int).tolist() in hcore.lattice.vectors:
                    # print(f"\t\t{vertices[i]} and {vertices[j]} are in contact")
                    ax.plot(*np.vstack(([vertices[i]],
                                        [vertices[j]])).T,
                            c="red",
                            linewidth=contact_width,
                            zorder=1)
                # else:
                #     print(f"{vertices[i] - vertices[j]} not in lattice!")
        # ax.axis("equal")
        # ax.axis("off")

    @classmethod
    def plot_fit(cls, ax, hcore, hp_sequence, coord_dict, dot_s=60, contact_width=1, chain_width=2):
        cls.plot_hcore(ax, hcore, dot_s, contact_width)
        coord_dict = {int(k.split("_")[1]): v for k, v in coord_dict.items()}
        coordinates = np.zeros((len(hp_sequence), len(hcore.basis[0])))
        coordinates[list(coord_dict.keys())] = list(coord_dict.values())
        ax.plot(*coordinates.T,
                c="black",
                linewidth=chain_width,
                zorder=3)
        ax.scatter(*coordinates[np.array(list(hp_sequence)) != "H"].T,
                   facecolors="white",
                   edgecolors="black",
                   s=100,
                   zorder=4)
