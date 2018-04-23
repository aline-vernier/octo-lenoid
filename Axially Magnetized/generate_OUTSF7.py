import os
import subprocess


def FullOUTSF(_D, _d, _h, _BR, _partNum, _w_dir, _e_val):

    intro_text = ""
    with open('intro_text.txt') as input_file:
        for line in input_file:
            intro_text += line

###################################################################
#                MAGNET PARAMS IN CM                              #
###################################################################

    D = 100*_D
    d = 100*_d
    h = 100*_h
    BR = 1e4*_BR

###################################################################
#                BOUNDARIES OF MAGNET IN PSF                      #
###################################################################

    p1 = [d/2., -h/2.]
    p2 = [D/2., -h/2.]
    p3 = [D/2., h/2.]
    p4 = [d/2., h/2.]
    points = [p1, p2, p3, p4, p1]

    magnet = "&reg mat=3, mtid=6, mshape=1 &  ! Magnet \n"
    for point in points:
        magnet = magnet + "&po x=" + str(point[0]) + ",y=" + str(point[1]) + "&\n"

    magnet += "&mt mtid=6\n"
    magnet += "HCEPT=" + str(-BR) + ", BCEPT=" + str(BR) + ", aeasy=90& ! From HKCM data"

###################################################################
#               GENERATE FILE                                     #
###################################################################

    text = intro_text + magnet

    _w_dir = '{0}\\{1}MeV\\{2}\\'.format(_w_dir, _e_val, _partNum)
    _w_name = _partNum

    am_file = _w_dir + _w_name + ".am"
    t35_file = _w_dir + _w_name + ".t35"

    if not os.path.exists(_w_dir):
        os.makedirs(_w_dir)

    with open(am_file, 'w') as out:
        out.write(text)

###################################################################
#               CALL AUTOMESH, PANDIRA AND  SF7                   #
###################################################################

    subprocess.call(
        ['G:\Programmes\LANL\AUTOMESH.EXE', am_file])
    subprocess.call(
        ['G:\Programmes\LANL\PANDIRA.EXE', t35_file])
    subprocess.call(
        ['G:\Programmes\LANL\SF7.EXE', t35_file])

    return _w_dir + 'OUTSF7.TXT'
