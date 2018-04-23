###################################################################
#         Code for cylindrical AXIALLY-MAGNETIZED SOLENOID        #
#                        ALL UNITS ARE SI                         #
#    RETURNS WRITTEN FILES ADDRESS                                #
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


def psf_generate(params_from_database, magnet_file, PSF_W_DIR, partNum, map_exists, scanParam,
                 fixedParam, e_energy, plot_field):

    subFolName = ''

    if fixedParam == 'Energy':
        subFolName = 'MeV'
    if fixedParam == 'Lmap':
        subFolName = 'mm'

    w_dir = '{0}\\{1}{2}'.format(PSF_W_DIR,
                                      fixedParam,
                                      subFolName)

    ###################################################################
    #      MAGNET PARAMETERS for FIELD MAP GENERATION, in meters      #
    ###################################################################
    if params_from_database:
        df = pd.read_excel(magnet_file, index_col=0)
        [h, d, D, BR] = [df.loc[partNum, 'h']*1e-3,
                         df.loc[partNum, 'd']*1e-3,
                         df.loc[partNum, 'D']*1e-3,
                         df.loc[partNum, 'BR']]

    else: # Enter parameters by hand here if required
        # outer diameter in m
        D = 20e-3
        # inner diameter in m
        d = 8e-3
        # height in m
        h = 30.0e-3
        # Remanent field in Tesla
        BR = 1.17
    ###################################################################
    #                ELECTRON PARAMS                                  #
    ###################################################################


    def gamma(energy):
        return energy*ee*1e6/(mo*(c*c))
    def beta(energy):
        return 1/gamma(energy)*m.sqrt(abs(gamma(energy)*gamma(energy) - 1))
    def vel(energy):
        return beta(energy)*c
    def momentum(energy):
        return gamma(energy)*mo*beta(energy)*c



    ###################################################################
    #                MAGNETIC FIELD AND FIELD INTEGRALS               #
    ###################################################################

    if map_exists == 0:
        out_sf7_file = FullOUTSF(D, d, h, BR, partNum, w_dir, e_energy)
        [z, BzNum] = psf.BImport(partNum, w_dir, e_energy)

        BzNumSq = np.multiply(BzNum, BzNum)
    else:
        [z, BzNum] = psf.BImport(partNum, w_dir, e_energy)
        BzNumSq = np.multiply(BzNum, BzNum)
        out_sf7_file = '{0}\\{1}MeV\\{2}\\OUTSF7.TXT'.format(w_dir, e_energy, partNum)

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

    # bmap_offset is the distance between the edge of the field map and the entrance of the magnet.
    # The current size of the field map in z is 120mm, hence 60-h/2
    # The value of the edge of the field map is defined in intro_text.txt
    bmap_offset = (0.06 - h/2.)

    return [out_sf7_file, w_dir, bmap_offset]


