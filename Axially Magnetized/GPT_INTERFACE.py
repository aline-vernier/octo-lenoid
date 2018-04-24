###################################################################
#         Code for cylindrical AXIALLY-MAGNETIZED SOLENOID        #
#                        ALL UNITS ARE SI                         #
#     Poisson superfish maps must exist for this script to run    #
#     SCAN OVER ENERGY AND DISTANCE OF SOLENOID TO SOURCE         #
#     NOT YET DEBUGGED FOR SINGLE RUN                             #
#                                                                 #
#     EXPLORATION MODE ALLOWS WIDE SCAN OF LMAP                   #
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


def gpt_run(MAGNET_FILE, OUTSF7_LOC, psf_w_dir, partNum, bmap_offset, GPT_ROOT, scanParam,
                 fixedParam, set_energy, show_plot, explore_mode):
    zero = 0
    ###################################################################
    #                PROGRAMME CONFIG                                 #
    ###################################################################

    # Root location of GPT simulations
    subFolName = ''
    if fixedParam == 'Energy':
        subFolName = 'MeV'
    if fixedParam == 'Lmap':
        subFolName = 'mm'

    w_dir = '{0}\\{1}{2}\\{3}MeV\\{4}'.format(GPT_ROOT, fixedParam, subFolName, set_energy, partNum)
    e_val = set_energy

    # Number of steps for source-solenoid distance
    num_l_map = 5

    # Scans to be performed
    multi_run = 1
    energy_scan = 0
    l_map_scan = 1
    l_map_steps = 1
    l_map_max = 0
    l_map_min = 0

    # Screen or snapshot
    snapshot = 1
    time_steps = 202
    if snapshot == 1:
        screen_or_snapshot = "snapshot(0,Tdet,Tstep) ;\n"
    else:
        screen_or_snapshot = "screen(\"wcs\",\"I\",Ldet);\n"

    # Data grouping if multirun
    group_by_time = 1

    # Focal length exploration mode : set to 1 if wide scan is required
    explore = 0

    ###################################################################
    #                MAGNET PARAMS FOR FIELD MAP GENERATION           #
    ###################################################################

    magnet_file = MAGNET_FILE
    df = pandas.read_excel(magnet_file, index_col=0)
    [h, rSol] = [df.loc[partNum, 'h'], df.loc[partNum, 'D']/2]

    ###################################################################
    #                ELECTRON BEAM PARAMS                             #
    ###################################################################

    # Beam parameters # energy in MeV
    e_energy_start = set_energy - 0.5
    e_energy_stop = set_energy + 0.5
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

    [z, BzNum] = psf.BImport(partNum, psf_w_dir, e_val)
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
    #               FOCAL LENGTH                                      #
    ###################################################################

    def f(set_energy):
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
    #               GENERATE Lmap array manually if required          #
    ###################################################################
    if explore == 0:
        # l_map_array = np.linspace(max(l_map_array), min(l_map_array) + 0.08, 10)
        l_map_array = np.linspace(0.0395, 0.04, 5)
    ###################################################################
    #               GENERATE GPT INPUT FILE                           #
    #          SCAN PARAMS = E_Scan and L_map optional                #
    ###################################################################

    if fixedParam == "Energy":
        energy = str(1000 * set_energy)
    else:
        energy = "EScan"

    if l_map_scan == 0:
        Lmap = str(min(l_map_array) + 5e-3)
        Lpinhole = str(min(l_map_array) + bmap_offset -1e-3)
    else:
        Lmap = "Lmap"
        Lpinhole = "Lmap + " + str(bmap_offset - 1e-3)
        print 'bmap_offset = {0}'.format(bmap_offset)

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
        + "# Start initial beam \n " \
        + "npart = " + str(npart) + ";\n\n" \
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
        + "Ldet = 0.5;\n" \
        + "Rpinhole = 1000e-6;	# diameter of entrance pinhole \n\n" \
        + "rmax(\"wcs\",\"z\"," + Lpinhole + ",Rpinhole, 0.5e-3);		# thickness of lead sheet\n" \
        + "rmax(\"wcs\",\"z\"," + str(zero) + ", " + str(rSol*1e-3) + ");		# thickness of lead sheet\n" \
        + "map2D_B(\"wcs\",\"z\"," + Lmap + ",\"" + w_dir + "\\fieldmap.gdf\",\"R\",\"Z\",\"Br\",\"Bz\",1.0);\n\n" \
        + "# Specify output times \n" \
        + "dtmax=1e-3/vz;\n" \
        + "Tdet=Ldet/vz;\n" \
        + "Tstep=Tdet/" + str(time_steps) +";\n" \
        + screen_or_snapshot

    GPT_input_file = w_dir + "\\SalleNoire_beam.in"
    with open(GPT_input_file, 'w') as out:
        out.write(text)

    ###################################################################
    #               GENERATE GPT MULTIPLE RUN FILE .MR                #
    #               IF MULTI_RUN == 1                                 #
    ###################################################################
    if multi_run == 1:
        l_map_min = min(l_map_array)
        l_map_max = 2*max(l_map_array)
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
        if group_by_time == 1:
            gdfa_cmd = " time"
        GPT_mr_file = w_dir + "\\SalleNoire_beam.mr"
        with open(GPT_mr_file, 'w') as out:
            out.write(text)

        ###################################################################
        #               GENERATE BATCH FILE for MULTIRUN                  #
        ###################################################################

        bat_text = "fish2gdf -o \"{0}\\fieldmap.gdf\" \"{1}\"\ngdf2a -o \"{0}\\fieldmap.txt\"" \
                   " \"{0}\\fieldmap.gdf\"\n" \
                   "mr -v -o \"{0}\\results_SalleNoire_beam.gdf\" \"{0}\\SalleNoire_beam.mr\" " \
                   "gpt \"{0}\\SalleNoire_beam.in\"\n" \
                   "gdfa -o \"{0}\\std_SalleNoire_beam.gdf\" " \
                   "\"{0}\\results_SalleNoire_beam.gdf\" {2} stdx numpar stdz avgz nemirrms\n" \
                   "gdf2a -o \"{0}\\std_SalleNoire_beam.txt\" " \
                   "\"{0}\\std_SalleNoire_beam.gdf\"\n\n".format(w_dir, OUTSF7_LOC, gdfa_cmd)
    else:

        ###################################################################
        #               GENERATE BATCH FILE FOR SIMPLE RUN                #
        ###################################################################

        gdfa_cmd = ""

        bat_text = "fish2gdf -o \"{0}\\fieldmap.gdf\" \"{1}\"\ngdf2a -o \"{0}\\fieldmap.txt\"" \
                   " \"{0}\\fieldmap.gdf\"\n" \
                   "gpt -v -o \"{0}\\results_SalleNoire_beam.gdf\" \"{0}\\SalleNoire_beam.in\"\n" \
                   "gdfa -o \"{0}\\std_SalleNoire_beam.gdf\" " \
                   "\"{0}\\results_SalleNoire_beam.gdf\" {2} stdx numpar stdz avgz nemirrms \n"\
            .format(w_dir, OUTSF7_LOC, gdfa_cmd)

    ###################################################################
    #               DEFINE AND PREPEND HEADER                         #
    ###################################################################

    n_lmap = m.modf((l_map_max - l_map_min)/l_map_steps)[1] + 1
    header = 'multirun {0} \nScanned_Over_Energy {1} \nEScan_Steps  {2} \nScanned_Over_Lmap {3} \nn_L_map {4}\n'.\
        format(str(multi_run), str(energy_scan), str(num_e_steps), str(l_map_scan), str(n_lmap)) \
        + 'group_by {0} \n'.format(str(gdfa_cmd)) \
        + 'time_steps {0} \n'.format(str(time_steps))\
        + 'electron_energy {0} \n'.format(str(energy))\
        + 'l_map {0} \n'.format(str(Lmap))\
        + 'l_pinhole {0} \n\n'.format(str(Lpinhole))

    GPT_head_file = w_dir + "\\std_SalleNoire_beam_h.txt"
    with open(GPT_head_file, 'w') as out:
        out.write(header)

    bat_text += 'type std_SalleNoire_beam.txt >> std_SalleNoire_beam_h.txt \n'

    GPT_bat_file = w_dir + "\\SalleNoire_beam.bat"
    with open(GPT_bat_file, 'w') as out:
        out.write(bat_text)

    ###################################################################
    #               RUN                                               #
    ###################################################################

    subprocess.call(['G:\Programmes\GPT\GPTwin\GPTwin.exe', GPT_bat_file])

    ###################################################################
    #               END OF SCRIPT                                     #
    ###################################################################
    return [GPT_head_file, w_dir]