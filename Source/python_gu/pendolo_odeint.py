import numpy as np

# =============================================================================
# def pend(y, t, b, c):
#     theta, omega = y
#     dydt = [omega, -b*omega - c*np.sin(theta)]
#     return dydt
# =============================================================================
dydt=np.zeros(2)
def pend(y, t, b, c):
    
    theta, omega = y
    #dydt = [omega, -b*omega - c*np.sin(theta)]
    dydt[0]=omega
    dydt[1]=-b*omega - c*np.sin(theta)
    return dydt


b = 0.25
c = 5.0

y0 = [np.pi - 0.1, 0.0]

t = np.linspace(0, 10, 101)

from scipy.integrate import odeint
sol = odeint(pend, y0, t, args=(b, c))

import matplotlib.pyplot as plt
plt.figure(dpi=200)
plt.plot(t, sol[:, 0], 'b', label=r'$\theta $(t)')
plt.plot(t, sol[:, 1], 'g', label=r'$\omega$(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()