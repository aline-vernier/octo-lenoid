import numpy as np


def BImport(partnum):
    partnum = '\\' + partnum
    w_dir = "G:\Programmes\LANL\Solenoids" + partnum
    w_file = "\OUTSF7.TXT"
    data = []

    print w_dir + w_file
    with open(w_dir + w_file) as my_file:
        for line in my_file:
            data.append(
                line.split()
            )

    BField_data = data[33::]

    minBounds = map(float, data[27][2][1:-1].split(','))
    maxBounds = map(float, data[28][2][1:-1].split(','))
    r_z_inc = [float(x) for x in [data[29][-1], data[29][-2]]]

    rBounds = [minBounds[0], maxBounds[0]]
    zBounds = [minBounds[1], maxBounds[1]]

    data = []
    for element in BField_data:
        data.append([float(x) for x in element])
    data = np.array(data).transpose()

    _r = 1e-2*data[0]
    _z = 1e-2*data[1]
    _Br = 1e-4*data[2]
    _Bz = 1e-4*data[3]
    _B = 1e-4*data[4]

    r0 = np.where(min(abs(_r)) == abs(_r))[0]
    B0 = np.array([abs(_Bz[element]) for element in r0])
    z0 = np.array([_z[element] for element in r0])

    return [z0, B0]
