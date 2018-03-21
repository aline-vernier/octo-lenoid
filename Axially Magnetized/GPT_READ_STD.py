###################################################################
#         Code that reads sdt file generated from GPT             #
#                        ALL UNITS ARE SI                         #
#                                                                 #
#                                                                 #
###################################################################

from matplotlib import pyplot as plt
import math as m
import numpy as np
import pandas as pd
from scipy.constants import elementary_charge as ee
from scipy.constants import m_e as mo
from scipy.constants import c
from scipy.interpolate import UnivariateSpline as UnivariateSpline
import re

###################################################################
#                IMPORT FILE GENERATED FROM GDFA                  #
###################################################################

partNum = "\9963-65252"
w_dir = "G:\GPT\Salle_Noire\Axially_Magnetized_Solenoid" + partNum
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
        l_map_vals.append(float(file_data[idx_1*(idx_2_max + 1)][1]))
        for idx_2 in range(2, idx_2_max):
            this_element = map(float,
                               file_data[idx_1*(idx_2_max + 1) + idx_2])
            data_temp.append(this_element)
            data.append(data_temp)

data = np.array(data)
print labels
