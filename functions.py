import numpy as np

def trajectory_time(params, t0=0):
    q0, qf, dq_m, ddq_m = params

    dq = qf-q0
    # triangular check
    c = np.sqrt(dq*ddq_m)

    if c <= dq_m:                 # Triangular
        t1 = np.sqrt(dq/ddq_m)
        T = t1
        tf = 2*t1
    else:                         # Trapezoidal
        t1 = dq_m/ddq_m
        T = dq/dq_m
        tf = T+t1

    return t0, t1, T, tf


def plan_trajectory(q_params, t_params):
    t0, t1, T, tf = t_params
    q0, qf, dq_m, ddq_m = q_params
    t = np.linspace(t0, tf, int(3E3))
    q = []
    v = []
    a = []

    for i in t:

        if i <= t1:
            qi = q0 + (0.5*ddq_m*(i-t0)**2)
            q02 = qi
            vi = ddq_m*i
            v02 = vi
            ai = ddq_m

        elif i > t1 and i <= T:
            vi = dq_m
            qi = q02 + v02*(i-t1)
            ai = 0

        elif i > T:
            vi = ddq_m*(tf-i)
            qi = qf - (0.5*ddq_m*(i-tf)**2)
            ai = -ddq_m

        q.append(qi)
        v.append(vi)
        a.append(ai)

    return t, q, v, a

def quintic_polynom(poly_params, t0, tf, known_members = 0):
    t = np.linspace(t0, tf, int(3E3))
    q = []
    v = []
    a = []

    for i in t:

        qi = (poly_params[0] * i ** 5) + (poly_params[1] * i ** 4) + (poly_params[2] * i ** 3) + (poly_params[3] * i ** 2) + (poly_params[4] * i) + poly_params[5]
        vi = (5 * poly_params[0] * i ** 4) + (4 * poly_params[1] * i ** 3) + (3 * poly_params[2] * i ** 2) + (2 * poly_params[3] * i) + poly_params[4] 
        ai = (20 * poly_params[0] * i ** 3) + (12 * poly_params[1] * i ** 2) + (6 * poly_params[2] * i)  + 2 * poly_params[3] 
        
        q.append(qi)
        v.append(vi)
        a.append(ai)

    return t, q, v, a

def cubic_polynom(poly_params, t0, tf, known_members = 0):
    t = np.linspace(t0, tf, int(3E3))
    q = []
    v = []
    a = []

    for i in t:

        qi = (poly_params[0] * i ** 3) + (poly_params[1] * i ** 2) + (poly_params[2] * i) + poly_params[3]
        vi = (3 * poly_params[0] * i ** 2) + (2 * poly_params[1] * i) + poly_params[2]
        ai = (6 * poly_params[0] * i ) + 2 * poly_params[1]  
        
        q.append(qi)
        v.append(vi)
        a.append(ai)

    return t, q, v, a
