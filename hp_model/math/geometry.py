import numpy as _np


class GeometryMath:
    @classmethod
    def project_point_on_line(cls, point, linePointAndVector):
        r0 = _np.array(point)
        r1 = _np.array(linePointAndVector[0])
        s1 = _np.array(linePointAndVector[1])
        return r1 - ((r1 - r0).dot(s1.T))/s1.dot(s1.T) * s1

    @classmethod
    def get_axis_rotation_matrix(cls, vector, angle):
        if _np.linalg.norm(vector) != 1:
            vector = vector / _np.linalg.norm(vector)

        cf = _np.cos(angle)
        sf = _np.sin(angle)
        vf = 1 - cf
        x, y, z = vector

        return _np.array([[x*x*vf + cf, x*y*vf - z*sf, x*z*vf + y*sf],
                            [x*y*vf + z*sf, y*y*vf + cf, y*z*vf - x*sf],
                            [x*z*vf - y*sf, y*x*vf + x*sf, z*z*vf + cf]])

    @classmethod
    def get_plane_by_3_points(cls, p0, p1, p2, verbose=False):
        cols12 = _np.vstack((p1 - p0, p2 - p0)).T
        if verbose:
            print('cols12:', cols12)
            print(cols12[[1,2]])
            print(cols12[[0,2]])
            print(cols12[[0,1]])

        A = _np.linalg.det(cols12[[1,2]])
        B = _np.linalg.det(cols12[[0,2]])
        C = _np.linalg.det(cols12[[0,1]])

        koeffs = _np.array([A,-B,C])
        koeffs = _np.append(koeffs, -koeffs.dot(p0))

        return koeffs

    @classmethod
    def get_plane_by_vector_and_point(cls, vec, p):
        plane = _np.array(vec)
        plane = _np.append(plane, -vec.dot(p))
        return plane

    @classmethod
    def project_point_on_plane(cls, p, plane):
        a,b,c = plane[:3]
        x0, y0, z0 = p
        A = _np.array([[b,-a,0],
                       [c,0,-a],
                       [a,b,c]])
        b = _np.array([b*x0 - a*y0,
                       c*x0 - a*z0,
                       -plane[-1]])
        return _np.linalg.solve(A, b)
