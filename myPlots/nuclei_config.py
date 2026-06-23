"""
nuclei_config.py

Plain-data configuration: which files exist, and how each nucleus / sim
model should be styled. This replaces the inline `f_Data_He =
ROOT.TFile.Open(...)` + `Series(...)` chains that used to live at the top of
makePlots.py -- adding a 6th nucleus or a 6th sim model is now a single
entry here, not a copy-pasted block in the driver script.

This module only holds data -- it does not open any ROOT files itself
(driver scripts call graph_io.open_file(entry['file']) when they actually
need one), so importing this module can never fail because a file is
temporarily missing.
"""

DATA_DIR = '/work/clas12/users/awild/RGM/CLAS_note/Boston2024/Main/Dec2024/Sys_Err/'
SIM_DIR = '/work/clas12/users/nwright/RGM/Q2/Sys_Err/rootFiles/'

# The reference nucleus everything else in the A-dependence plots is
# normalized to.
REFERENCE_KEY = 'He'

# One entry per nucleus used in the cross-nucleus (A-dependence) plots.
# 'offset' is a small x-shift so markers for different nuclei don't sit
# exactly on top of each other in the overlay plots; 'marker' is the
# matplotlib marker style.
NUCLEI = [
    {'key': 'He',   'label': r'${}^{4}He$',   'file': DATA_DIR + 'Data_He.root',
     'color': 'black',           'marker': None, 'offset': 0.0},
    {'key': 'C',    'label': r'${}^{12}C$',   'file': DATA_DIR + 'Data_C.root',
     'color': 'darkred',         'marker': '<',  'offset': -0.01},
    {'key': 'Ca40', 'label': r'${}^{40}Ca$',  'file': DATA_DIR + 'Data_40Ca.root',
     'color': 'darkgreen',       'marker': 'v',  'offset': -0.005},
    {'key': 'Ca48', 'label': r'${}^{48}Ca$',  'file': DATA_DIR + 'Data_48Ca.root',
     'color': 'darkblue',        'marker': '^',  'offset': 0.005},
    {'key': 'Sn',   'label': r'${}^{120}Sn$', 'file': DATA_DIR + 'Data_Sn.root',
     'color': 'darkgoldenrod',   'marker': '>',  'offset': 0.01},
]

# One entry per simulation model used in the He data-vs-sim plots.
SIM_MODELS = [
    {'key': 'AV18',  'label': 'AV18',         'file': SIM_DIR + 'Sim_He_AV18.root',
     'color': 'red'},
    {'key': 'AV4',   'label': "AV4'",         'file': SIM_DIR + 'Sim_He_AV4.root',
     'color': 'pink'},
    {'key': 'N2LO10', 'label': 'N2LO (1.0fm)', 'file': SIM_DIR + 'Sim_He_N2LO10.root',
     'color': 'blue'},
    {'key': 'N2LO12', 'label': 'N2LO (1.2fm)', 'file': SIM_DIR + 'Sim_He_N2LO12.root',
     'color': 'green'},
    {'key': 'NV',    'label': 'Norfolk',      'file': SIM_DIR + 'Sim_He_NV.root',
     'color': 'orange'},
]


def find(entries, key):
    for e in entries:
        if e['key'] == key:
            return e
    raise KeyError("No entry with key %r in this config list" % key)


def reference_nucleus():
    return find(NUCLEI, REFERENCE_KEY)
