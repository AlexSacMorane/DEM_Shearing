# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import shutil
from datetime import datetime
from pathlib import Path
import pickle
import matplotlib.pyplot as plt

#Own function and class
import Create_IC
import Create_IC_Polygonal
import Report
import User

#-------------------------------------------------------------------------------
#IC simulation
#-------------------------------------------------------------------------------

if Path('Debug').exists():
    shutil.rmtree('Debug')
os.mkdir('Debug')
os.mkdir('Debug/Configuration')
os.mkdir('Debug/Configuration/Init')

simulation_report = Report.Report('Debug/Report.txt',datetime.now())

#get data
dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations = User.All_parameters()

#find name
i_run = 1
folderpath = Path('ICs/'+dict_algorithm['template_simulation_name']+str(i_run)+'_dict_ic')
while folderpath.exists():
    i_run = i_run + 1
    folderpath = Path('ICs/'+dict_algorithm['template_simulation_name']+str(i_run)+'_dict_ic')
dict_algorithm['name_folder'] = dict_algorithm['template_simulation_name']+str(i_run)

#initial ic
Create_IC.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#-------------------------------------------------------------------------------
#Shear simulation
#-------------------------------------------------------------------------------

#define periodic bc

#define packs
i_bottom = 0
i_top = 0
for grain in dict_ic['L_g_tempo'] :
    if grain.is_group(0, 2*dict_geometry['R_mean'], 'Bottom') :
        i_bottom = i_bottom + 1
    elif grain.is_group(dict_sample['y_box_max']-2*dict_geometry['R_mean'], dict_sample['y_box_max'], 'Top') :
        i_top = i_top + 1
simulation_report.write_and_print(str(i_bottom)+' grains in Bottom group\n'+str(i_top)+' grains in Top group\n', str(i_bottom)+' grains in Bottom group\n'+str(i_top)+' grains in Top group')
