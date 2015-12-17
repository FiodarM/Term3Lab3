__author__ = 'fiodar'

import numpy as np
from numpy import array, append, apply_along_axis

def odesolve(f, y0, sequence, method='runge_kutta', tolerance=1e-2):

    sequence = append([], sequence)
    t = t0 = sequence[0]
    y = array(y0)
    ytable = array([y0])
    dt = (sequence[-1] - t0) * tolerance
    steps = 0

    if method == 'euler':

        next_y = lambda dx: y + f(t, y) * dx
        rest = lambda dx: apply_along_axis(abs, 0, array(next_y(dx) - next_y(dx/2.)))

    else:

        k1 = lambda: f(t, y)
        k2 = lambda dx: f(t+dx/2., y+k1()*dx/2.)
        k3 = lambda dx: f(t+dx/2., y+k2(dx)*dx/2.)
        k4 = lambda dx: f(t+dx, y+k3(dx)*dx)
        next_y = lambda dx: y + dx/6.*(k1() + 2*k2(dx) + 2*k3(dx) + k4(dx))
        rest = lambda dx: apply_along_axis(abs, 0, array(next_y(dx) - next_y(dx/2.)))/15.

    for end in sequence[1:]:

        while abs(t - t0) < abs(end - t0):

            if np.equal(rest(dt), rest(dt/2.)).all():
                dt = end - t

            elif np.greater_equal(rest(dt), tolerance * apply_along_axis(abs, 0, next_y(dt/2.))).any():
                while np.greater_equal(rest(dt), tolerance * apply_along_axis(abs, 0, next_y(dt/2.))).any()\
                        and abs(dt) > 1e-6:

                    dt *= 0.5

            else:
                while np.less(rest(dt), tolerance * apply_along_axis(abs, 0, next_y(dt/2.))).all():
                    value = np.less(rest(dt), tolerance * apply_along_axis(abs, 0, next_y(dt/2.)))
                    dt *= 1.2
                    if abs(t + dt - t0) > abs(end - t0):
                        dt = end - t
                        break

            y = next_y(dt)
            t += dt
            steps += 1

        t0 = t

        ytable = append(ytable, [y], axis=0)

    print method, steps

    return ytable




