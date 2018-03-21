###################################################################
#         Code for cylindrical AXIALLY-MAGNETIZED SOLENOID        #
#                        ALL UNITS ARE SI                         #
###################################################################

from matplotlib import pyplot as plt
import math as m
import numpy as np
import psf_import as psf
from scipy.constants import elementary_charge as ee
from scipy.constants import m_e as mo
from scipy.constants import c
from scipy.interpolate import UnivariateSpline as UnivariateSpline
from generate_OUTSF7 import FullOUTSF
import pandas as pd

###################################################################
#                      PARTNUM                                    #
###################################################################

partNum = "9963-65252"

###################################################################
#                PROGRAMME CONFIG                                 #
###################################################################

psf_map = 1
map_exists = 0
plot_field = 0
params_from_database = 1
magnet_file = "G:\Programmes\LANL\Solenoids\hkcm_magnets.xlsx"

###################################################################
#                MAGNET PARAMS FOR FIELD MAP GENERATION           #
###################################################################

# outer diameter in m
D = 20e-3
# inner diameter in m
d = 8e-3
# height in m
h = 30.0e-3
# Remanent field in Tesla
BR = 1.17

###################################################################
#                MAGNET PARAMETERS IN METERS                      #
###################################################################
if params_from_database:
    df = pd.read_excel(magnet_file, index_col=0)
    [h, d, D] = [df.loc[partNum, 'h']*1e-3,
                   df.loc[partNum, 'd']*1e-3,
                   df.loc[partNum, 'D']*1e-3]
else:
    Ri = 4.5e-3
    Ro = 23.5e-3
    L  = 6e-3

###################################################################
#                ELECTRON PARAMS                                  #
###################################################################

# Beam parameters # energy in MeV, R is twice the RMS size of beam
e_energy = 2
R = 2e-3

def gamma(energy):
    return energy*ee*1e6/(mo*(c*c))
def beta(energy):
    return 1/gamma(energy)*m.sqrt(abs(gamma(energy)*gamma(energy) - 1))
def vel(energy):
    return beta(energy)*c
def momentum(energy):
    return gamma(energy)*mo*beta(energy)*c

print "Electron energy = " + str(e_energy) + " MeV"
print "Gamma = " + str(gamma(e_energy))

###################################################################
#                MAGNETIC FIELD AND FIELD INTEGRALS               #
###################################################################


def Bz(z_, L_, Ri_, Ro_):
    p = z_ + L_ / 2
    n = z_ - L_ / 2
    t1 = p / m.sqrt(Ri_ * Ri_ + p * p)
    t2 = -n / m.sqrt(Ri_ * Ri_ + n * n)
    t3 = -p / m.sqrt(Ro_ * Ro_ + p * p)
    t4 = n / m.sqrt(Ro_ * Ro_ + n * n)
    return BR / 2 * (t1 + t2 + t3 + t4)

if psf_map == 0:
    # z-axis in m
    z = 1e-3 * np.array(range(-1000, 1000, 1))
    BzNumSq = [Bz(x, L, Ri, Ro) * Bz(x, L, Ri, Ro) for x in z]
    BzNum = [Bz(x, L, Ri, Ro) for x in z]

else:
    if map_exists == 0 :
        FullOUTSF(D, d, h, BR, partNum)
        [z, BzNum] = psf.BImport(partNum)
        BzNumSq = np.multiply(BzNum, BzNum)
    else:
        [z, BzNum] = psf.BImport(partNum)
        BzNumSq = np.multiply(BzNum, BzNum)


# Interpolation function for Bz, BzSquared and second derivative
BzInterp = UnivariateSpline(z, BzNum, s=0, k=5)

# Interpolation function for BzSquared
BzSqInterp = UnivariateSpline(z, BzNumSq)
BzInterp_2d = BzInterp.derivative(n=2)

# Array from 2nd derivative of spline
BzNum_2d = BzInterp(z)

# Product of B" and B^2
Bz2d_Bz = np.multiply(BzNum_2d, BzNum)
Bz2d_Bz_interp = UnivariateSpline(z, Bz2d_Bz)

# Field Integrals
integral_F2 = BzSqInterp.integral(min(z), max(z))
integral_F3 = Bz2d_Bz_interp.integral(min(z), max(z))


###################################################################
#               OPTIONAL PLOT OF B FIELD AT ZERO                  #
###################################################################
if plot_field == 1:
    plt.plot(z, BzNum)
    plt.plot(z, BzInterp(z))
    plt.show()


###################################################################
#               FOCAL LENGTH AND EMITTANCE GROWTH                 #
###################################################################

def f():
    return pow(ee*ee/(4*gamma(e_energy)*gamma(e_energy)*mo*mo*vel(e_energy)*vel(e_energy))*integral_F2, -1)

print "Magnetic field at center of solenoid = " + str(max(BzNum)) + " T"
print "Focal Length = " + str(f()*1e3) + " mm"

# Emittance growth from spherical aberration in mm.mrad (hence 1e6 factor)
def emittance(R_, p_):
    return 1e6/(mo*c)*ee*ee*m.pow(R_, 4)/(m.sqrt(2)*48*momentum(e_energy))*integral_F3

print "Spherical Aberration induced Emittance growth = " + str(emittance(R, momentum(e_energy))) + " mm.mrad"

