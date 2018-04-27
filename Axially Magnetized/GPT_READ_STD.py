###################################################################
#         Code that reads sdt file generated from GPT             #
#                        ALL UNITS ARE SI                         #
#                                                                 #
#                                                                 #
###################################################################

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pristinifier as ps
import pandas as pd
from matplotlib.lines import Line2D
import matplotlib.patches as patches

def read_std(GPT_HEAD_FILE, w_dir, partNum, energy, bmap_offset):
    ###################################################################
    #                      PARTNUM                                    #
    ###################################################################

#    partNum = "9963-67475"

    ###################################################################
    #                LOAD MAGNET DATA                                 #
    ###################################################################

    magnet_file = "G:\Programmes\LANL\Solenoids\hkcm_magnets.xlsx"
    df = pd.read_excel(magnet_file, index_col=0)
    [h, d, D] = [df.loc[partNum, 'h'],
                 df.loc[partNum, 'd'],
                 df.loc[partNum, 'D']]
    ###################################################################
    #                IMPORT FILE GENERATED FROM GDFA                  #
    ###################################################################

#    w_dir = "G:\GPT\Salle_Noire\Axially_Magnetized_Solenoid" + "\\" + partNum
#    w_file = "\\std_SalleNoire_beam_h.txt"

    file_data = []
    with open(GPT_HEAD_FILE) as my_file:
        for line in my_file:
            file_data.append(line.split())


    multirun = bool(file_data[0][1])
    energy_scan = bool(file_data[1][1])
    n_e_steps = float(file_data[2][1])
    l_map_scan = bool(file_data[3][1])
    n_lmap_steps = float(file_data[4][1])
    group_by = file_data[5][1]
    n_t_steps = int(file_data[6][1])
    energy = file_data[7][1]
    Lmap = file_data[8][1]
    Lpinhole = float(file_data[9][3])

    file_data = file_data[10:-1]
    data = []
    data_temp = []
    l_map_vals = []
    l_map_vals_float = []

    pos_mag_vals = []
    l_pinhole_vals = []
    labels = []
    if file_data[n_t_steps+3] != []:
        n_t_steps += 1
    if multirun:
        idx_2_max = n_t_steps + 3
        for idx_1 in range(len(file_data) / (n_t_steps + 3)):

            data_temp = []
            l_map_vals.append(
                file_data[idx_1*idx_2_max + 1 + idx_1][0] + ' = '
                + str(float(file_data[idx_1*idx_2_max + 1][1])) + ' cm')

            l_map_vals_float.append(float(file_data[idx_1*idx_2_max + 1][1]))
            Lmapval  = float(file_data[idx_1*idx_2_max + 1][1])
            Lmagnet_t = Lmapval*1e3 + 60 - h/2.

            pos_mag_vals.append(
                'magpos = '
                + str(1e-1*Lmagnet_t) + ' cm')

            print 'Lmap = {2}, Lmagnet = {1}, Lpinhole = {0}'.format(Lpinhole*1e3 + Lmapval*1e3, Lmapval*1e3 + 60 -h/2., Lmapval*1e3)
            for idx_2 in range(4, idx_2_max):
                this_element = map(float,
                                   file_data[idx_1*idx_2_max + idx_2])
                data_temp.append(this_element)

            data.append(data_temp)

    data = np.array(data)


    ###################################################################
    #                PLOT DATA FROM FILE                              #
    ###################################################################

    params = ps.plot_params()
    plt.rcParams.update(params)

    [tableau20, tableau20Edge] = ps.rgb_array()

    fig, ax = plt.subplots(3, sharex=True)
    mpl.interactive(False)
    # stdx vs avgz
    for idx in range(len(l_map_vals)):
        stdx = [1e3*element for element in np.transpose(data[idx])[1]]
        charge = np.transpose(data[idx])[2]
        avgz = np.transpose(data[idx])[4]
        nemirrms = [1e6*element for element in np.transpose(data[idx])[5]]

        tidx = idx%len(tableau20)
        ax[0].plot(avgz, stdx, color=tableau20[idx], marker='o', markeredgecolor=tableau20Edge[idx])
        ax[0].legend(pos_mag_vals, loc='upper right')
        ax[1].plot(avgz, nemirrms, color=tableau20[idx], marker='o', markeredgecolor=tableau20Edge[idx])
        ax[2].semilogy(avgz, charge, color=tableau20[idx], marker='o', markeredgecolor=tableau20Edge[idx])

    for idx in range(len(l_map_vals)):

        Lp = l_map_vals_float[idx] + Lpinhole
        line = Line2D([Lp, Lp], [0, 5], color=tableau20[idx])
        line2 = Line2D([Lp++1e-3*h, Lp++1e-3*h], [0, 5], color=tableau20[idx])
        ax[0].add_line(line)
        ax[0].add_line(line2)

#    ax[0].set_ylim(0, 2.5)
    ax[1].set_ylim(0, 5)

    ax[2].set_xlabel('z (m)')

    ax[0].set_ylabel('stdx (mm)')
    ax[1].set_ylabel('emittance (mm.mrad)')
    ax[2].set_ylabel('num particles')

    fig.show()
    fig.savefig(w_dir + '\\Lmap_stdx_and_nemirstd_vs_avgz.jpg', bbox_inches='tight')
