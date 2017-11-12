#!/usr/bin/python3

"""
File name      : hu_moments.py
Description    : This program computes Hu's moments invariant on scale and rotation
                 of a discrete image f(x,y).
Reference      : Ming-Kuei Hu, "Visual pattern recognition by moment invariants,"
                 in IRE Transactions on Information Theory, vol. 8, no. 2, pp. 179-187,
                 February 1962.
                 doi: 10.1109/TIT.1962.1057692
                 URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1057692&isnumber=22787
Date created   : 12/november/2017
Created by     : Adailton Júnior
E-mail         : adailtonjn68@gmail.com
Github project : https://github.com/adailtonjn68/hu_moments_in_python/
Python version : 3.5 (tested in version 3.5 but as there're no big changes in python
                 standard, it can be used with other versions).
Comments       : This implementation uses only python standard functions, but it
                 can be used also with numpy arrays.
"""


def m_pq(f, p, q):
    """
    Two-dimensional (p+q)th order moment of image f(x,y)
    where p,q = 0, 1, 2, ...
    """
    m = 0
    # Loop in f(x,y)
    for x in range(0, len(f)):
        for y in range(0, len(f[0])):
            # +1 is used because if it wasn't, the first row and column would
            # be ignored
            m += ((x+1)**p)*((y+1)**q)*f[x][y]
    return m


def centroid(f):
    """
    Computes the centroid of image f(x,y)
    """
    m_00 = m_pq(f, 0, 0)
    return [m_pq(f, 1, 0)/m_00, m_pq(f, 0 ,1)/m_00]


def u_pq(f, p, q):
    """
    Centroid moment invariant to rotation.
    This function is equivalent to the m_pq but translating the centre of image
    f(x,y) to the centroid.
    """
    u = 0
    centre = centroid(f)
    for x in range(0, len(f)):
        for y in range(0, len(f[0])):
            u += ((x-centre[0]+1)**p)*((y-centre[1]+1)**q)*f[x][y]
    return u


def hu(f):
    """
    This function computes Hu's seven invariant moments.
    """
    u_00 = u_pq(f, 0, 0)

    # Scale invariance is obtained by normalization.
    # The normalized central moment is given below
    eta = lambda f, p, q: u_pq(f, p, q)/(u_00**((p+q+2)/2))

    # normalized central moments used to compute Hu's seven moments invariat
    eta_20 = eta(f, 2, 0)
    eta_02 = eta(f, 0, 2)
    eta_11 = eta(f, 1, 1)
    eta_12 = eta(f, 1, 2)
    eta_21 = eta(f, 2, 1)
    eta_30 = eta(f, 3, 0)
    eta_03 = eta(f, 0, 3)

    # Hu's moments are computed below
    phi_1 = eta_20 + eta_02
    phi_2 = 4*eta_11 + (eta_20-eta_02)**2
    phi_3 = (eta_30 - 3*eta_12)**2 + (3*eta_21 - eta_03)**2
    phi_4 = (eta_30 + eta_12)**2 + (eta_21 + eta_03)**2
    phi_5 = (eta_30 - 3*eta_12)*(eta_30 + eta_12)*((eta_30+eta_12)**2 - 3*(eta_21+eta_03)**2) + (3*eta_21 - eta_03)*(eta_21 + eta_03)*(3*(eta_30 + eta_12) - (eta_21 + eta_03)**2)
    phi_6 = (eta_20 - eta_02)*((eta_30 + eta_12)**2 - (eta_21 + eta_03)**2) + 4*eta_11*(eta_30 + eta_12)*(eta_21 + eta_03)
    phi_7 = (3*eta_21 - eta_03)*(eta_30 + eta_12)*((eta_30 + eta_12)**2 - 3*(eta_21 + eta_03)**2) - (eta_30 - 3*eta_12)*(eta_21 + eta_03)*(3*(eta_30 + eta_12)**2 - (eta_21 + eta_03)**2)

    return [phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7]
