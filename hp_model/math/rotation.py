import numpy as np


class RotationMath:
    def __init__(self):
        pass

    @classmethod
    def get_axis_rotation_matrix(cls, vector, angle):
        if np.linalg.norm(vector) != 1:
            vector = vector / np.linalg.norm(vector)

        cf = np.cos(angle)
        sf = np.sin(angle)
        vf = 1 - cf
        x, y, z = vector

        return np.array([[x * x * vf + cf, x * y * vf - z * sf, x * z * vf + y * sf],
                         [x * y * vf + z * sf, y * y * vf + cf, y * z * vf - x * sf],
                         [x * z * vf - y * sf, y * x * vf + x * sf, z * z * vf + cf]])

    @classmethod
    def rotation_matrix(cls, angle, rot_axis=np.array([1, 0])):
        if rot_axis.shape[0] not in [2, 3]:
            raise Exception("Wrong n_dim value!")
        if rot_axis.shape[0] == 2:
            return np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]]) # .astype(int)
        return cls.get_axis_rotation_matrix(rot_axis, angle)

    @classmethod
    def mirroring_matrix(cls, axis_u):
        axis_u = np.array(axis_u) / np.linalg.norm(axis_u)
        axis_u = axis_u.reshape((np.reshape(axis_u, -1).shape[0], 1))
        return np.eye(axis_u.shape[0]) - 2 * axis_u @ axis_u.T

    @classmethod
    def rotate_points(cls, points, angle, rot_axis=np.array([1,0])):
        m = cls.rotation_matrix(angle, rot_axis)
        return np.vstack([np.round(m @ v) for v in points]).astype(int)

    @classmethod
    def order_vertices(cls, core):
        return np.array(sorted(core))
