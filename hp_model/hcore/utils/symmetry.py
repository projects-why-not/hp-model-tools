from ...math import np, RotationMath, ContactMath
from ..hcore import HCore


class SymmetryInspection:
    __angles = [-np.pi / 2,
                np.pi / 2,
                np.pi,
                0]

    def __init__(self):
        pass

    @classmethod
    def _analyze_rotation(cls, core1, core2, angle, rot_axis=np.array([1,0])):
        def analyze_parallel_shift(c1, c2):
            # if ContactMath.compute_contact_num(c1, lattice.vectors) != ContactMath.compute_contact_num(c2, lattice.vectors):
            #     return False
            if c1.n_contacts != c2.n_contacts:
                return False

            ordrd_rot_core2 = RotationMath.order_vertices(c2.vertices.tolist())
            ordrd_core1 = RotationMath.order_vertices(c1.vertices.tolist())
            #         print("ordrds:\n", ordrd_core1, "\n", ordrd_rot_core2, "\n\n")
            offsets = ordrd_rot_core2 - ordrd_core1
            # print(f"OFFSETS FOR {core1}, {core2}: {offsets}")

            # discovering offset
            if (np.unique(offsets[:, 0]).shape[0] == 1) and (np.unique(offsets[:, 1]).shape[0] == 1):
                return True
            return False

        def analyze_mirroring(c1, c2, rel_ax_ind):
            # m_core1 = np.array(c1).astype(int)
            # m_core2 = np.array(c2).astype(int)
            # if rel_to_ox:
            #     corner_agg_fs = ((np.max, np.max),
            #                      (np.max, np.min))
            # else:
            #     corner_agg_fs = ((np.min, np.max),
            #                      (np.max, np.max))
            #
            # p0_1 = m_core1[m_core1[:, 0] == corner_agg_fs[0][0](m_core1[:, 0])]
            # p0_1 = p0_1[p0_1[:, 1] == corner_agg_fs[0][1](p0_1[:, 1])]
            # p0_1 = p0_1[0]
            # p0_2 = m_core2[m_core2[:, 0] == corner_agg_fs[1][0](m_core2[:, 0])]
            # p0_2 = p0_2[p0_2[:, 1] == corner_agg_fs[1][1](p0_2[:, 1])]
            # p0_2 = p0_2[0]
            #
            # m_core1 -= p0_1  # m_core1[0]
            # m_core2 -= p0_2  # m_core2[0]
            # m_core2[:, int(rel_to_ox)] *= -1
            #
            # # FIXME: needed or not?
            # # if ContactMath.compute_contact_num(m_core1, lattice.vectors) != ContactMath.compute_contact_num(m_core2, lattice.vectors):
            # #     return False
            #
            # m_core1_set = set(["_".join(v.astype(str).tolist()) for v in m_core1])
            # m_core2_set = set(["_".join(v.astype(str).tolist()) for v in m_core2])
            # # print("Ox", m_core1_set, m_core2_set)
            # if m_core1_set == m_core2_set:
            #     return True
            # return False
            axis = np.zeros(c1.vertices[0].shape[0])
            axis[rel_ax_ind] = 1
            mirror_matr = RotationMath.mirroring_matrix(axis)
            mirrored_c2 = HCore(np.array([mirror_matr @ node for node in c2.vertices]),
                                c2.lattice)
            return analyze_parallel_shift(c1, mirrored_c2)

        lattice = core2.lattice
        # core1_n_contacts = core1.n_contacts

        # core1 = np.array(core1.vertices)
        core2_nodes = np.array(core2.vertices)
        # core1[:,0] -= min(core1[:,0])
        # core1[:, 1] -= min(core1[:, 1])
        core2_nodes[:, 0] -= min(core2_nodes[:, 0])
        core2_nodes[:, 1] -= min(core2_nodes[:, 1])

        rot_core2_nodes = RotationMath.rotate_points(core2_nodes, angle, rot_axis)
        rot_core2 = HCore(rot_core2_nodes, lattice)
        if analyze_parallel_shift(core1, rot_core2):
            return True

        # discovering mirroring
        # if analyze_mirroring(core1, rot_core2, True) or analyze_mirroring(core1, rot_core2, False):
        #     return True
        for ax_ind in range(core1.vertices[0].shape[0]):
            if analyze_mirroring(core1, rot_core2, ax_ind):
                return True

        return False

    @classmethod
    def are_equal(cls, hcore1, hcore2, rotation_angles=None):
        # print(rotation_angles)
        # input("c?")
        if rotation_angles is None:
            return cls._analyze_rotation(hcore1, hcore2, 0)
        # rotation + offset + x-, y-symmetry
        if "3D" in hcore1.lattice.kind:
            for axis in range(3):
                ax_vec = np.zeros(3)
                ax_vec[axis] = 1
                for ang in rotation_angles:
                    if cls._analyze_rotation(hcore1, hcore2, ang, ax_vec):
                        return True
        else:
            for ang in rotation_angles:
                if cls._analyze_rotation(hcore1, hcore2, ang):
                    return True

        return False
