from matplotlib import pyplot as plt
import numpy as np

def rgb_array():
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    tableau20Edge = list(tableau20)

    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)


    for i in range(len(tableau20Edge)):
        r, g, b = tableau20Edge[i]
        r = r * 1/2
        g = g * 1/2
        b = b * 1/2
        tableau20Edge[i] = (r / 255., g / 255., b / 255.)

    return [tableau20, tableau20Edge]




def plot_params():
    fig_width_pt = 4*246.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0 / 72.27  # Convert pt to inches
    golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = 1*fig_width * 1  # height in inches
    fig_size = [fig_width, fig_height]

    params = {'font.size' : 20,
              'font.family' : 'Helvetica',
              'font.monospace' : 'Computer Modern',
              'axes.labelsize' : 15,
              'backend' : 'ps',
              'legend.fontsize': 10,
              'xtick.labelsize' : 16,
              'ytick.labelsize' : 16,
              'text.usetex': True,
              'figure.figsize' : fig_size}
    return params

# plt.rcParams.update(params)

# fig, ax = plt.subplots(1, sharex=True)

# ax[0].plot(_x_array, _y_array, color=tableau20[0], markeredgecolor=tableau20Edge[0])
# ax[0].plot(xPlot, ref, color=tableau20[4], markeredgecolor=tableau20Edge[6])

# ax[0].errorbar(xRaw - RefX0, refRaw,fmt='o', yerr=refErr, color='black', markerfacecolor='none', markersize=5)
# ax[0].errorbar(xRaw - RefX0, atomRaw, fmt='d', yerr=atomErr, color='black', markerfacecolor='none', markersize=5)
# ax[0].plot([RefX0 - RefX0, RefX0 - RefX0], [0, 140], '--', color=tableau20[4])
# ax[0].plot([AtomX0 - RefX0, AtomX0 - RefX0], [0, 140], '--', color=tableau20[2])
