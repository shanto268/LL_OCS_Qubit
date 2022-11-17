from pyEPR.calcs import Convert
import numpy as np
import matplotlib.pyplot as plt

def get_Ec_from_Ej(Ej, target_omega):
    a = np.sqrt(8*Ej)
    w = target_omega
    Ec_minus = 0.5*(a**2 - 2*w - np.sqrt(a**4 - 4*a**2 *w))
    #Ec_plus = 0.5 *(a**2 - 2*w + np.sqrt(a**4 - 4*a**2 *w))
    return Ec_minus #GHz

def get_Cj_from_Lj(Ljs, target_omega):
    Ejs = Convert.Ej_from_Lj(Ljs, 'nH', "GHz")
    Ecs = get_Ec_from_Ej(Ejs, target_omega=3) #GHz
    Cjs = Convert.Cs_from_Ec(Ecs, units_in='GHz', units_out='fF')
    return Cjs
