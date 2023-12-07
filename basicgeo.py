import math
from collections import namedtuple


# cheap verson of normalization
def octahedron_angle(x, y, z):
    s = abs(x) + abs(y) + abs(z)
    if z >= 0:
        return x/s, y/s
    if x >= 0:
        if y >= 0:
            return 1 - y/s, 1 - x/s
        return 1 + y/s, -1 + x/s
    if y >= 0:
        return -1 + y/s, 1 + x/s
    return -1 - y/s, -1 - x/s


class P2(namedtuple('P2', ['u', 'v'])):
    __slots__ = ()

    def __new__(self, u, v):
        return super(P2, self).__new__(self, float(u), float(v))

    def __repr__(self):
        return "P2(%s, %s)" % (self.u, self.v)

    def __add__(self, a):
        return P2(self.u + a.u, self.v + a.v)

    def __sub__(self, a):
        return P2(self.u - a.u, self.v - a.v)

    def __mul__(self, a):
        return P2(self.u*a, self.v*a)

    def __neg__(self):
        return P2(-self.u, -self.v)

    def __rmul__(self, a):
        raise TypeError

    def lensq(self):
        return self.u*self.u + self.v*self.v

    def len(self):
        if self.u == 0.0:
            return abs(self.v)
        if self.v == 0.0:
            return abs(self.u)
        return math.sqrt(self.u*self.u + self.v*self.v)

    def len_l_z(self):
        return math.sqrt(self.x*self.x + self.y*self.y)

    def arg(self):
        return math.degrees(math.atan2(self.v, self.u))

    def assertlen1(self):
        assert abs(self.len() - 1.0) < 0.0001
        return True
        
    @staticmethod
    def dot(a, b):
        return a.u*b.u + a.v*b.v

    @staticmethod
    def dot_l_z(a, b):
        return a.u*b.x + a.v*b.y

    @staticmethod
    def z_norm(v):
        ln = v.len()
        if ln == 0.0:  
            ln = 1.0
        return P2(v.u/ln, v.v/ln)
            
    @staticmethod
    def a_perp(v):
        return P2(-v.v, v.u)

    @staticmethod
    def c_perp(v):
        return P2(v.v, -v.u)

    @staticmethod
    def convert_l_z(p):
        return P2(p.x, p.y)


class P3(namedtuple('P3', ['x', 'y', 'z'])):
    __slots__ = ()

    def __new__(self, x, y, z):
        return super(P3, self).__new__(self, float(x), float(y), float(z))

    def __repr__(self):
        return "P3(%s, %s, %s)" % (self.x, self.y, self.z)

    def __add__(self, a):
        return P3(self.x + a.x, self.y + a.y, self.z + a.z)

    def __sub__(self, a):
        return P3(self.x - a.x, self.y - a.y, self.z - a.z)

    def __mul__(self, a):
        return P3(self.x*a, self.y*a, self.z*a)

    def __neg__(self):
        return P3(-self.x, -self.y, -self.z)

    def __rmul__(self, a):
        raise TypeError

    def lensq(self):
        return self.x*self.x + self.y*self.y + self.z*self.z

    def len(self):
        return math.sqrt(self.lensq())

    def len_l_z(self):
        return math.sqrt(self.x*self.x + self.y*self.y)

    def assertlen1(self):
        assert abs(self.len() - 1.0) < 0.0001
        return True
        
    @staticmethod
    def dot(a, b):
        return a.x*b.x + a.y*b.y + a.z*b.z

    @staticmethod
    def cross(a, b):
        return P3(a.y*b.z - b.y*a.z, -a.x*b.z + b.x*a.z, a.x*b.y - b.x*a.y)

    @staticmethod
    def diff(a, b, bfore):
        if bfore:
            return b - a
        return a - b
        
    @staticmethod
    def z_norm(v):
        ln = v.len()
        if ln == 0.0:  
            ln = 1.0
        return P3(v.x/ln, v.y/ln, v.z/ln)
            
    @staticmethod
    def convert_g_z(p, z):
        return P3(p.u, p.v, z)

    @staticmethod
    def convert_c_z(p, z):
        return P3(p.x, p.y, z)


# should this also be a named tuple?
class I1:
    def __init__(self, lo, hi):
        self.lo = lo
        self.hi = hi

    def absorb(self, v):
        assert self.lo <= self.hi
        if v < self.lo:
            self.lo = v
        elif v > self.hi:
            self.hi = v

    def __repr__(self):
        return "I1(%s, %s)" % (self.lo, self.hi)

    def d_is_valid(self):
        assert self.lo <= self.hi
        return self.lo <= self.hi

    def leng(self):
        assert self.d_is_valid()
        return self.hi - self.lo

    def along(self, lam):
        assert self.d_is_valid()
        return self.lo * (1 - lam) + self.hi * lam

    def contains(self, v):
        assert self.d_is_valid()
        return self.lo <= v <= self.hi

    def contains_strict(self, v):
        assert self.d_is_valid()
        return self.lo < v < self.hi

    def inflate(self, v):
        assert self.d_is_valid()
        self.lo -= v
        self.hi += v
        assert self.d_is_valid()

    @staticmethod
    def absorb_list(vs):
        v0 = next(vs)
        res = I1(v0, v0)
        for v in vs:
            res.absorb(v)
        return res


class Partition1:
    def d_along(self, lam):  # shouldn't be here
        return self.lo * (1 - lam) + self.hi * lam

    def __init__(self, lo, hi, nparts):
        self.lo = lo
        self.hi = hi
        self.nparts = nparts
        self.vs = [self.d_along(i * 1.0 / nparts) for i in range(0, nparts + 1)]
        assert (lo, hi) == (self.vs[0], self.vs[-1])
        assert len(self.vs) == nparts + 1
        
    def get_part(self, v):
        assert self.lo <= v <= self.hi, ("Getpart", v, "not between", self.lo, self.hi)
        i = int(self.nparts * (v - self.lo) / (self.hi - self.lo))
        # print(i, v)
        if i != 0 and v < self.vs[i]:
            i -= 1
        if i == self.nparts:
            i -= 1
        elif i < self.nparts - 1 and v >= self.vs[i + 1]:
            i += 1
        assert self.vs[i] <= v <= self.vs[i+1], (v, "not between", self.vs[i], self.vs[i+1])
        assert i == self.nparts - 1 or v < self.vs[i+1]
        return i

    def get_part_range(self, vlo, vhi):
        if vhi < self.lo or vlo > self.hi:
            return 0, -1
            
        ilo = int(self.nparts * (vlo - self.lo) / (self.hi - self.lo))
        if ilo > 0:
            if vlo <= self.vs[ilo]:
                ilo -= 1
        else:
            ilo = 0
        
        ihi = int(self.nparts * (vhi - self.lo) / (self.hi - self.lo))
        if ihi < self.nparts - 1:
            if vhi > self.vs[ihi + 1]:
                ihi += 1
        else:
            ihi = self.nparts - 1
            
        assert 0 <= ilo <= ihi <= self.nparts - 1
        return ilo, ihi
        

class Quat:
    def __init__(self, q0, q1, q2, q3):
        self.q0 = q0
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3
        self.iqsq = 1/((self.q0**2 + self.q1**2 + self.q2**2 + self.q3**2) or 1)

    def __mul__(self, a):
        return Quat(a.q0*self.q0 - a.q1*self.q1 - a.q2*self.q2 - a.q3*self.q3, 
                    a.q0*self.q1 + a.q1*self.q0 + a.q2*self.q3 - a.q3*self.q2, 
                    a.q0*self.q2 - a.q1*self.q3 + a.q2*self.q0 + a.q3*self.q1, 
                    a.q0*self.q3 + a.q1*self.q2 - a.q2*self.q1 + a.q3*self.q0) 
    
    def vec_dots(self):
        r00 = self.q0*self.q0*2 * self.iqsq
        r11 = self.q1*self.q1*2 * self.iqsq
        r22 = self.q2*self.q2*2 * self.iqsq
        r33 = self.q3*self.q3*2 * self.iqsq
        r01 = self.q0*self.q1*2 * self.iqsq
        r02 = self.q0*self.q2*2 * self.iqsq
        r03 = self.q0*self.q3*2 * self.iqsq
        r12 = self.q1*self.q2*2 * self.iqsq
        r13 = self.q1*self.q3*2 * self.iqsq
        r23 = self.q2*self.q3*2 * self.iqsq
        return (P3(r00 - 1 + r11,   r12 + r03,         r13 - r02),
                P3(r12 - r03,       r00 - 1 + r22,     r23 + r01),
                P3(r13 + r02,       r23 - r01,         r00 - 1 + r33))

    def vec_dots_0(self):
        r00 = self.q0*self.q0*2 * self.iqsq
        r11 = self.q1*self.q1*2 * self.iqsq
        r02 = self.q0*self.q2*2 * self.iqsq
        r03 = self.q0*self.q3*2 * self.iqsq
        r12 = self.q1*self.q2*2 * self.iqsq
        r13 = self.q1*self.q3*2 * self.iqsq
        return P3(r00 - 1 + r11,   r12 + r03,         r13 - r02)
        
    def vec_dots_1(self):
        r00 = self.q0*self.q0*2 * self.iqsq
        r22 = self.q2*self.q2*2 * self.iqsq
        r01 = self.q0*self.q1*2 * self.iqsq
        r03 = self.q0*self.q3*2 * self.iqsq
        r12 = self.q1*self.q2*2 * self.iqsq
        r23 = self.q2*self.q3*2 * self.iqsq
        return P3(r12 - r03,       r00 - 1 + r22,     r23 + r01)
        
    def vec_dots_2(self):
        r00 = self.q0*self.q0*2 * self.iqsq
        r33 = self.q3*self.q3*2 * self.iqsq
        r01 = self.q0*self.q1*2 * self.iqsq
        r02 = self.q0*self.q2*2 * self.iqsq
        r13 = self.q1*self.q3*2 * self.iqsq
        r23 = self.q2*self.q3*2 * self.iqsq
        return P3(r13 + r02,       r23 - r01,         r00 - 1 + r33)
                                            
    def vec_dots_t(self):   # transposed
        r00 = self.q0*self.q0*2 * self.iqsq
        r11 = self.q1*self.q1*2 * self.iqsq
        r22 = self.q2*self.q2*2 * self.iqsq
        r33 = self.q3*self.q3*2 * self.iqsq
        r01 = self.q0*self.q1*2 * self.iqsq
        r02 = self.q0*self.q2*2 * self.iqsq
        r03 = self.q0*self.q3*2 * self.iqsq
        r12 = self.q1*self.q2*2 * self.iqsq
        r13 = self.q1*self.q3*2 * self.iqsq
        r23 = self.q2*self.q3*2 * self.iqsq
        return (P3(r00 - 1 + r11,   r12 - r03,         r13 + r02),
                P3(r12 + r03,       r00 - 1 + r22,     r23 - r01),
                P3(r13 - r02,       r23 + r01,         r00 - 1 + r33))

    def __repr__(self):
        return "Quat(%f,%f,%f,%f)" % (self.q0, self.q1, self.q2, self.q3)


def along_acc(lam, a, b):
    if a == b:
        return a
    return a * (1 - lam) + b * lam


def along(lam, a, b):
    return a * (1 - lam) + b * lam


sendactivity = None


def set_sendactivity(lsendactivity):
    global sendactivity
    sendactivity = lsendactivity

    
def dplotrect(xlo, xhi, ylo, yhi, materialnumber=0):
    sendactivity("contours", contours=[[(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi), (xlo, ylo)]],
                 materialnumber=materialnumber)
