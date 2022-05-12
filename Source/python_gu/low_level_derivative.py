import CoolProp

HEOS = CoolProp.AbstractState("HEOS", "Water")

HEOS.update(CoolProp.PT_INPUTS, 101325, 300)

HEOS.cpmass()

HEOS.first_partial_deriv(CoolProp.iHmass, CoolProp.iT, CoolProp.iP)


# See how much faster this is?
CoolProp.CoolProp.PropsSI('d(Hmass)/d(T)|P', 'P', 101325, 'T', 300, 'Water')


HEOS.keyed_output(CoolProp.iDmass)
HEOS.rhomass()

HEOS.first_partial_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP)
HEOS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iHmass)
