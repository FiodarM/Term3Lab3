__author__ = 'fiodar'

import pylab as pl
from scipy import integrate
from lab3 import *


def f(t, y):

    a = lambda t, x, v: pl.exp(-x)
    x, v = y

    return pl.array([v, a(t, x, v)])

f_sp = lambda y, t: f(t, y)

T = pl.linspace(0., 5.0, 30)

Xeu = odesolve(f, [0., 0.], T, method='euler', tolerance=1e-2).transpose()
Xrk = odesolve(f, [0., 0.], T, method='runge_kutta', tolerance=1e-2).transpose()
Xsp = integrate.odeint(f_sp, [0., 0.], T).transpose()

style = 'r-'
style_sp = 'k--'
width = 2

pl.figure('coordinate (time)', tight_layout=1)

pl.xlabel('$$t$$')
pl.ylabel('$$x$$')
pl.plot(T, Xrk[0], style, linewidth=width)
pl.plot(T, Xsp[0], style_sp, linewidth=width)
pl.grid()

pl.figure('velocity (coordinate)',  tight_layout=1)

pl.xlabel('x')
pl.ylabel('v')
pl.plot(Xrk[0], Xrk[1], style, linewidth=width)
pl.plot(Xsp[0], Xsp[1], style_sp)
pl.axhline(pl.sqrt(2), 0, Xrk[0][-1], color='g', ls='--', linewidth=width)
pl.text(0.1, 1.45, '$$v_stat$$', color='g')
pl.grid()


pl.show()
