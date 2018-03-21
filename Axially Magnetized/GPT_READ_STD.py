###################################################################
#         Code that reads sdt file generated from GPT             #
#                        ALL UNITS ARE SI                         #
#                                                                 #
#                                                                 #
###################################################################


from matplotlib import pyplot as plt
import numpy as np
import pristinifier as ps
import pandas as pd

###################################################################
#                      PARTNUM                                    #
###################################################################

partNum = "9963-65252"

###################################################################
#                IMPORT FILE GENERATED FROM GDFA                  #
###################################################################

w_dir = "G:\GPT\Salle_Noire\Axially_Magnetized_Solenoid" + "\\" + partNum
w_file = "\\std_SalleNoire_beam_h.txt"

file_data = []
with open(w_dir + w_file) as my_file:
    for line in my_file:
        file_data.append(line.split())

# print file_data[0:7]

multirun = bool(file_data[0][1])
energy_scan = bool(file_data[1][1])
n_e_steps = float(file_data[2][1])
l_map_scan = bool(file_data[3][1])
n_lmap_steps = float(file_data[4][1])
group_by = file_data[5][1]
n_t_steps = int(file_data[6][1])

file_data = file_data[8::]
data = []
data_temp = []
l_map_vals = []
labels = []

if multirun:
    idx_2_max = n_t_steps + 3
    labels = file_data[1]
    for idx_1 in range(len(file_data) % (n_t_steps + 3)):
        data_temp = []
        l_map_vals.append(
            file_data[idx_1 * (idx_2_max + 1)][0] + ' = '
        + str(1e2*float(file_data[idx_1*(idx_2_max + 1)][1])) + ' cm')
        for idx_2 in range(2, idx_2_max):
            this_element = map(float,
                               file_data[idx_1*(idx_2_max + 1) + idx_2])
            data_temp.append(this_element)
        data.append(data_temp)

data = np.array(data)
###################################################################
#                LOAD MAGNET DATA                                 #
###################################################################

magnet_file = "G:\Programmes\LANL\Solenoids\hkcm_magnets.xlsx"
df = pd.read_excel(magnet_file, index_col=0)
[h, d, D] = [df.loc[partNum, 'h'],
             df.loc[partNum, 'd'],
             df.loc[partNum, 'D']]
###################################################################
#                PLOT DATA FROM FILE                              #
###################################################################

params = ps.plot_params()
plt.rcParams.update(params)

[tableau20, tableau20Edge] = ps.rgb_array()

# stdx vs avgz
for idx in range(len(l_map_vals)):
    stdx = np.transpose(data[idx])[1]
    avgz = np.transpose(data[idx])[4]
    plt.plot(avgz, stdx, color=tableau20[idx], marker='o', markeredgecolor=tableau20Edge[idx])
plt.legend(l_map_vals, loc='upper left')
plt.title('magnet \# {0}, h = {1} mm, d = {2} mm, D = {3} mm'.format(partNum, h, d, D), fontsize = 15)

plt.savefig(w_dir+'\\Lmap_stdx_vs_avgz.pdf', bbox_inches='tight')
plt.show()



