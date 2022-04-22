import numpy as np, sympy as sp
# import mplpub
from simupy.systems.symbolic import DynamicalSystem, dynamicsymbols
from simupy.block_diagram import BlockDiagram
from simupy.array import Array, r_
import matplotlib.pyplot as plt

def shape_figure():

    pass


if __name__ == "__main__":

    x, v, u = dynamicsymbols('x v u')
    l, m = sp.symbols('l m')

    parameters = {

        l: 1,
        m: 1

    }

    inertia = DynamicalSystem(

        state_equation=r_[v, u / (m * l ** 2)],
        state=r_[x, v],
        input_= u,
        constants_values=parameters

    )

    g = sp.symbols('g')
    parameters[g] = 9.81

    gravity = DynamicalSystem(

        output_equation=-g * m * l * sp.sin(x),
        input_= x,
        constants_values = parameters

    )

    BD = BlockDiagram(inertia, gravity)
    BD.connect(gravity, inertia)
    BD.connect(inertia, gravity, outputs=[0])
    BD.simulate(8)

    print(BD)
