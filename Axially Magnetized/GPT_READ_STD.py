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
        this_line = line.split()
        file_data.append(this_line)

header = file_data[0:7]
data = file_data[8::]

print header

