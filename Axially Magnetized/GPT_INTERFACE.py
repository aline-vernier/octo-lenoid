###################################################################
#         Code for cylindrical AXIALLY-MAGNETIZED SOLENOID        #
#                        ALL UNITS ARE SI                         #
#     Poisson superfish maps must exist for this script to run    #
#     SCAN OVER ENERGY AND DISTANCE OF SOLENOID TO SOURCE         #
###################################################################

from matplotlib import pyplot as plt
import math as m
import numpy as np
import psf_import as psf
from scipy.constants import elementary_charge as ee
from scipy.constants import m_e as mo
from scipy.constants import c
from scipy.interpolate import UnivariateSpline as UnivariateSpline
import os
import pandas
import subprocess

###################################################################
#                PROGRAMME CONFIG                                 #
###################################################################
show_plot = 0
# Root location of GPT simulations
GPT_ROOT = "G:\GPT\Salle Noire"
SOL_TYPE = "\Axially Magnetized Solenoid"

# Number of steps for source-solenoid distance
num_l_map = 10

# Scans to be performed
multi_run = 1
energy_scan = 0
l_map_scan = 1

###################################################################
#                MAGNET PARAMS FOR FIELD MAP GENERATION           #
###################################################################

magnet_file = "G:\Programmes\LANL\Solenoids\hkcm_magnets.xlsx"
partNum = "9963-56987"
df = pandas.read_excel(magnet_file, index_col=0)
h = df.loc[partNum, 'h']

w_dir = GPT_ROOT + SOL_TYPE + "\\" + partNum

###################################################################
#                ELECTRON BEAM PARAMS                             #
###################################################################

# Beam parameters # energy in MeV
e_energy_start = 1.0
e_energy_stop = 3.0
num_e_steps = 3

# RMS angular divergence in rad
sigdiv = 60e-3
# Length of electron packet in m
length = 10e-6
# Transverse size of electron packet in m
sigr = 10e-6
# Number of particles
npart = 1000


def gamma(energy):
    return energy * ee * 1e6 / (mo * (c * c))
def beta(energy):
    return 1 / gamma(energy) * m.sqrt(abs(gamma(energy) * gamma(energy) - 1))
def vel(energy):
    return beta(energy) * c
def momentum(energy):
    return gamma(energy) * mo * beta(energy) * c

###################################################################
#                MAGNETIC FIELD AND FIELD INTEGRALS               #
###################################################################

[z, BzNum] = psf.BImport("\\" + partNum)
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
if show_plot == 1:
    plt.plot(z, BzNum)
    plt.plot(z, BzInterp(z))
    plt.show()


###################################################################
#               FOCAL LENGTH AND EMITTANCE GROWTH                 #
###################################################################

def f(e_energy):
    return pow(
        ee * ee / (4 * gamma(e_energy) * gamma(e_energy) * mo * mo * vel(e_energy) * vel(e_energy)) * integral_F2, -1)


###################################################################
#               GENERATE FOCAL LENGTH FILE                        #
###################################################################

# Create folder if it doesn't exist
if not os.path.exists(w_dir):
    os.makedirs(w_dir)

l_map_array = []
text = ""
for e_energy in np.linspace(e_energy_start, e_energy_stop, num=num_e_steps):
    f_length = f(e_energy)
    l_map = f_length + min(z)
    l_map_array.append(l_map)
    text = text + str(e_energy) + "\t" + str(f_length) + "\t" + str(l_map) + "\n"

focLength_file = w_dir + "\\focalLengths.dat"

with open(focLength_file, 'w') as out:
    out.write(text)

###################################################################
#               GENERATE GPT INPUT FILE                           #
#          SCAN PARAMS = E_Scan and L_map optional                #
###################################################################

if energy_scan == 0:
    energy = str(1000 * e_energy_start)
else:
    energy = "EScan"

if l_map_scan == 0:
    Lmap = str(min(l_map_array) + 5e-3)
else:
    Lmap = "Lmap"

text = \
    "# Magnet Number " + str(partNum) + '\n\n' \
    + "# Define beam parameters \n" \
    + "E0 = " + energy + "; # keV\n" \
    + "gamma=E0/511+1; \n" \
    + "Gbetaz=sqrt(gamma^2-1); \n" \
    + "vz=c*Gbetaz/gamma; \n" \
    + "sigdiv = " + str(sigdiv) + "; \n" \
    + "sigGbetar=sigdiv*Gbetaz;\n\n" \
    + "sigr = " + str(sigr) + "; \n" \
    + "len = " + str(length) + "; \n\n" \
    + "# Start initial beam \n npart = " + str(npart) + ";\n\n" \
    + "setparticles(\"beam\",npart,me,qe,0.0) ;\n\n" \
    + "# transverse \n" \
    + "setrxydist(\"beam\",\"g\",0,sigr,0,3) ;\n" \
    + "setphidist(\"beam\",\"u\",0,2*pi) ; \n" \
    + "setGBrxydist(\"beam\",\"g\",0,sigGbetar,0,3);" \
    + "setGBphidist(\"beam\",\"u\",0,2*pi);\n\n" \
    + "# longitudinal\n" \
    + "setzdist(\"beam\",\"u\",0,len);\n" \
    + "setGdist(\"beam\",\"u\",gamma,0) ;\n\n" \
    + "# Positions of various elements:  iris + the solenoid + detectors\n" \
    + "# Lmap = Scanned over \n" \
    + "# Liris = Scanned over \n" \
    + "Ldet = 0.3;\n" \
    + "Rpinhole = 1000e-6;	# diameter of entrance pinhole \n\n" \
    + "rmax(\"wcs\",\"z\"," + Lmap + ",Rpinhole,0.5e-3);		# thickness of lead sheet\n " \
    + "map2D_B(\"wcs\",\"z\"," + Lmap + ",\"fieldmap.gdf\",\"R\",\"Z\",\"Br\",\"Bz\",1.0);\n\n" \
    + "# Specify output times \n" \
    + "dtmax=1e-3/vz;\n" \
    + "Tdet=Ldet/vz;\n" \
    + "Tstep=Tdet/100;\n" \
    + "#screen(\"wcs\",\"I\",Ldet);\n" \
    + "snapshot(0,Tdet,Tstep) ;\n\n"

GPT_input_file = w_dir + "\\SalleNoire_beam.in"
with open(GPT_input_file, 'w') as out:
    out.write(text)

###################################################################
#               GENERATE GPT MULTIPLE RUN FILE .MR                #
#               IF MULTI_RUN == 1                                 #
###################################################################
if multi_run == 1:
    l_map_min = min(l_map_array)
    l_map_max = 2 * max(l_map_array)
    l_map_steps = (l_map_max - l_map_min) / float(num_l_map)
    e_energy_step = (e_energy_stop - e_energy_start) / float(num_e_steps)

    EScan_text = "EScan " + str(e_energy_start * 1e3) + "; " \
                 + str(e_energy_stop * 1e3) + " " \
                 + str(e_energy_step * 1e3) + "\n"
    Lmap_text = "Lmap " + str(l_map_min) + " " \
                + str(l_map_max) + " " \
                + str(l_map_steps) + "\n"
    gdfa_cmd = ""
    if energy_scan == 1 and l_map_scan == 1:
        text = EScan_text + Lmap_text
        gdfa_cmd = " EScan Lmap"
    if energy_scan == 1 and l_map_scan == 0:
        text = EScan_text
        gdfa_cmd = " EScan"
    if energy_scan == 0 and l_map_scan == 1:
        text = Lmap_text
        gdfa_cmd = " Lmap"

    GPT_mr_file = w_dir + "\\SalleNoire_beam.mr"
    with open(GPT_mr_file, 'w') as out:
        out.write(text)

    ###################################################################
    #               GENERATE BATCH FILE                               #
    ###################################################################

    text = "fish2gdf -o fieldmap.gdf " + "\"" + w_dir + "\OUTSF7.TXT\"\n"\
           + "gdf2a -o fieldmap.txt fieldmap.gdf\n"\
           + "mr -v -o results_SalleNoire_beam.gdf SalleNoire_beam.mr gpt SalleNoire_beam.in\n"\
           + "gdfa -o std_SalleNoire_beam.gdf results_SalleNoire_beam.gdf " + gdfa_cmd +" stdx stdy stdz avgz\n" \
           + "gdf2a -o std_SalleNoire_beam.txt std_SalleNoire_beam.gdf\n"\


    GPT_bat_file = w_dir + "\\SalleNoire_beam.bat"
    with open(GPT_bat_file, 'w') as out:
        out.write(text)

    ###################################################################
    #               RUN BATCH FILE IF MULTIPLE RUN                    #
    ###################################################################

    subprocess.call(['G:\Programmes\GPT\GPTwin\GPTwin.exe', GPT_bat_file])

else:

    ###################################################################
    #               GENERATE BATCH FILE                               #
    ###################################################################

    text = "fish2gdf -o fieldmap.gdf " + w_dir + "\OUTSF7.TXT\n" \
           + "gdf2a -o fieldmap.txt fieldmap.gdf\n" \
           + "gpt -v -o results_SalleNoire_beam.gdf SalleNoire_beam.in\n"

    GPT_bat_file = w_dir + "\\SalleNoire_beam.bat"
    with open(GPT_bat_file, 'w') as out:
        out.write(text)

    ###################################################################
    #               RUN BATCH FILE                                    #
    ###################################################################

    subprocess.call(['G:\Programmes\GPT\GPTwin\GPTwin.exe', GPT_bat_file])
    ###################################################################
    #               END OF SCRIPT                                     #
    ###################################################################
