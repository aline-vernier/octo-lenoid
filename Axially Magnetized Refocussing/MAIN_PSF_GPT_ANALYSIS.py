import PSF_INTERFACE as PSF
import GPT_INTERFACE as GPT
import GPT_READ_STD as GPT_READ
import pandas as pd

PSF_W_DIR = 'G:\Programmes\LANL\Solenoids'
GPT_W_DIR = 'G:\GPT\Salle_Noire\Axially_Magnetized_Solenoid\FullStudy'

scan_param = "Lmap"
fixed_param = "Energy"
energy = 3.0 #in MeV
l_map = None

plot_field = 0
map_exists = 0
explore_mode = 0

# params_from_database = 1
MAGNET_FILE = "G:\Programmes\LANL\Solenoids\hkcm_magnets.xlsx"

partNumLoop = False
params_from_database = True
genField_and_GPT = True

df = pd.read_excel(MAGNET_FILE, index_col=0)

if partNumLoop:
    for partNum in df.index.values:

        [OUTSF7_LOC, psf_w_dir, bmap_offset] = PSF.psf_generate(params_from_database, MAGNET_FILE, PSF_W_DIR, partNum,
                             map_exists, scan_param, fixed_param, energy, plot_field)

        [GPT_HEAD_FILE, FILE_LOC] = GPT.gpt_run(MAGNET_FILE, OUTSF7_LOC, psf_w_dir, partNum, bmap_offset, GPT_W_DIR,
                                                scan_param, fixed_param, energy, plot_field, explore_mode)

        GPT_READ.read_std(GPT_HEAD_FILE, FILE_LOC, partNum, energy, bmap_offset)

else:
    partNum = '0000-00004'
    if genField_and_GPT:
        [OUTSF7_LOC, psf_w_dir, bmap_offset] = PSF.psf_generate(params_from_database, MAGNET_FILE, PSF_W_DIR, partNum,
                                                                map_exists, scan_param, fixed_param, energy, plot_field)

        [GPT_HEAD_FILE, FILE_LOC] = GPT.gpt_run(MAGNET_FILE, OUTSF7_LOC, psf_w_dir, partNum, bmap_offset, GPT_W_DIR,
                                                scan_param, fixed_param, energy, plot_field, explore_mode)

        GPT_READ.read_std(GPT_HEAD_FILE, FILE_LOC, partNum, energy, bmap_offset)
    else:

        GPT_HEAD_FILE = "G:\GPT\Salle_Noire\Axially_Magnetized_Solenoid\FullStudy\EnergyMeV\\3.0MeV\9963-67475\std_SalleNoire_beam_h.txt"
        FILE_LOC = "G:\GPT\Salle_Noire\Axially_Magnetized_Solenoid\FullStudy\EnergyMeV\\3.0MeV\9963-67475"
        bmap_offset = 0.04

        GPT_READ.read_std(GPT_HEAD_FILE, FILE_LOC, partNum, energy, bmap_offset)


