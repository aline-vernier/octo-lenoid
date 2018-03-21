import os
import subprocess


def FullOUTSF(_D, _d, _h, _BR, _partNum ):
    _partNum = "\\" + _partNum
    intro_text = ""
    with open('intro_text.txt') as input_file:
        for line in input_file:
            intro_text += line


###################################################################
#                MAGNET PARAMS IN CM                              #
###################################################################

#    D = 4
#    d = 1.5
#    h = 3.0
#    BR = 12800

    D = 100*_D
    d = 100*_d
    h = 100*_h
    BR = 1e4*_BR

    print D
    print d
    print h

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
    partNum = _partNum
    w_dir = "G:\Programmes\LANL\Solenoids" + partNum
    w_name = "\AxiallyMagnetizedSolenoid"

    am_file = w_dir + w_name + ".am"
    t35_file = w_dir + w_name + ".t35"

    if not os.path.exists(w_dir):
        os.makedirs(w_dir)

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

    return
