import numpy as np
import CoolProp.CoolProp as CP

'caso in  cui si conservano le portate'

lib="HEOS"
fluidname= "CO2"
fluid   = CP.AbstractState(lib, fluidname )



dydt=np.zeros(2)

def serbatoio(Y, t, m_in, m_out_l, m_out_v, h_in, V):
    
    y, P = Y
        
    fluid.update(CP.PQ_INPUTS, P, 0)
    H_l=fluid.hmass()
    rho_l=fluid.rhomass()     
    de_rho_l = fluid.first_saturation_deriv(CP.iDmass, CP.iP)  #derivata parziale nella pressione
    de_H_l = fluid.first_saturation_deriv(CP.iHmass, CP.iP)  #derivata parziale nella pressione
    
    fluid.update(CP.PQ_INPUTS, P, 1)
    H_v=fluid.hmass()
    rho_v=fluid.rhomass()     
    de_rho_v= fluid.first_saturation_deriv(CP.iDmass, CP.iP)  #derivata parziale nella pressione
    de_H_v = fluid.first_saturation_deriv(CP.iHmass, CP.iP)  #derivata parziale nella pressione
    
    
    A= rho_l*H_l - rho_v*H_v
    B= y*(H_l*de_rho_l + rho_l*de_H_l)  +  (1-y)*(H_v*de_rho_v + rho_v*de_H_v) - 1
    C= rho_l - rho_v
    D= y*de_rho_l  +  (1-y)*de_rho_v 
    M=  (m_in - m_out_l - m_out_v)/V
    E=  (m_in*h_in - m_out_l*H_l - m_out_v*H_v)/V
    
    
    
    dydt[0]= (E-M*B/D)/(A-C*B/D)  #y
    dydt[1]= (M-dydt[0]*C)/D      #p
    
    
    fluid.update(CP.HmassP_INPUTS, h_in, P)
    q=fluid.Q()
    print('titolo = ',q)
    yyy= (1-q)*rho_v/(rho_l*q + rho_v*(1-q))
    print('y = ',yyy)
    
    return dydt

Q_start=0.2  #Q iniziale
P_start=40e5
#fluid.update(CP.PQ_INPUTS, P_start, Q_start)

'condizioni al cintorno in ingresso'
Q_in=0.5
fluid.update(CP.PQ_INPUTS, P_start, Q_in)

h_in=fluid.hmass()
m_in=        1
m_out_l=     m_in*(1-Q_start)
m_out_v=     m_in*Q_start
V=           10 #m^3

b = 0.25
c = 5.0

y0 = [0.5, P_start]

t = np.linspace(0, 50000, 101)

from scipy.integrate import odeint
sol = odeint(serbatoio, y0, t, args=(m_in, m_out_l, m_out_v, h_in, V))
y, P = sol[-1,:]

import matplotlib.pyplot as plt
plt.figure(dpi=200)
plt.plot(t, sol[:, 0], 'b', label=r'$y $(t)')
#plt.plot(t, sol[:, 1], 'g', label=r'$\omega$(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

plt.figure(dpi=200)
#plt.plot(t, sol[:, 0], 'b', label=r'$\theta $(t)')
plt.plot(t, sol[:, 1], 'g', label=r'$p$(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()









'riverifica:'
fluid.update(CP.HmassP_INPUTS, h_in, P)
print('titolo finale  = ',fluid.Q())
print('titolo iniziale= ',Q_start)


'verifica delta rho'
fluid.update(CP.PQ_INPUTS, P_start, 0)
rho_l=fluid.rhomass()     

fluid.update(CP.PQ_INPUTS, P_start, 1)
rho_v=fluid.rhomass()     

print('delta rho iniziale= ',rho_l-rho_v)

fluid.update(CP.PQ_INPUTS, P, 0)
rho_l=fluid.rhomass()     

fluid.update(CP.PQ_INPUTS, P, 1)
rho_v=fluid.rhomass()     

print('delta rho finale= ',rho_l-rho_v)

