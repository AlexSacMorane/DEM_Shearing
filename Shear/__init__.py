# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains ??.nt functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math
import matplotlib.pyplot as plt

#Own
import Create_IC
import Create_IC.Grain_ic
import Create_IC.Contact_gg_ic
import Create_IC.Contact_gimage_ic

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def DEM_vertical_load(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Loading the granular system with vertical load.

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simultion report (a report)
        Output :
            Nothing, but initial condition dictionnary is updated
    """
    i_DEM_0 = dict_ic['i_DEM_IC']
    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gimage'] = []
    dict_ic['L_contact_ij_gimage'] = []
    dict_ic['id_contact'] = 0
    dict_ic['L_g_image'] = []
    dict_ic['L_i_image'] = []

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    Ymax_stop = 0
    for grain in dict_ic['L_g_tempo']:
        Force_stop = Force_stop + 0.5*grain.mass*dict_sollicitation['gravity']
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(dict_ic['Ecin_ratio_IC']*grain.radius/dict_ic['dt_DEM_IC'])**2

    while DEM_loop_statut :

        dict_ic['i_DEM_IC'] = dict_ic['i_DEM_IC'] + 1

        #create image
        for grain in dict_ic['L_g_tempo']:
            #left wall
            if (grain.center[0] - dict_sample['x_box_min']) < dict_algorithm['d_to_image'] :
                if grain.id in dict_ic['L_i_image'] : #image exists
                    image = dict_ic['L_g_image'][dict_ic['L_i_image'].index(grain.id)]
                    if image.position == 'right' :
                        image.position = 'left'
                else : #image does not exist
                    dict_ic['L_g_image'].append(Create_IC.Grain_ic.Grain_Image(grain, 'left'))
                    dict_ic['L_i_image'].append(grain.id)
            #right wall
            elif (dict_sample['x_box_max'] - grain.center[0]) < dict_algorithm['d_to_image'] :
                if grain.id in dict_ic['L_i_image'] : #image exists
                    image = dict_ic['L_g_image'][dict_ic['L_i_image'].index(grain.id)]
                    if image.position == 'left' :
                        image.position = 'right'
                else : #image does not exist
                    dict_ic['L_g_image'].append(Create_IC.Grain_ic.Grain_Image(grain, 'right'))
                    dict_ic['L_i_image'].append(grain.id)
            #center
            else :
                if grain.id in dict_ic['L_i_image'] : #image exists
                    i_toremove = dict_ic['L_i_image'].index(grain.id)
                    dict_ic['L_g_image'].pop(i_toremove)
                    dict_ic['L_i_image'].pop(i_toremove)
        #translate image
        for image in dict_ic['L_g_image']:
            if image.position == 'left' :
                image.translation(np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0]))
            elif image.position == 'right' :
                image.translation(np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0]))

        #Contact detection
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_ic['i_update_neighborhoods_com']  == 0:
            Create_IC.Contact_gg_ic.Update_Neighborhoods(dict_ic)
            Create_IC.Contact_gimage_ic.Update_Neighborhoods(dict_ic)
        Create_IC.Contact_gg_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)
        Create_IC.Contact_gimage_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in  dict_ic['L_contact']+dict_ic['L_contact_gimage']:
            contact.normal()
            contact.tangential(dict_ic['dt_DEM_IC'])

        #Move grains (only Current)
        for grain in dict_ic['L_g_tempo']:
            if grain.group == 'Current' :
                grain.euler_semi_implicite(dict_ic['dt_DEM_IC'],10*dict_ic['Ecin_ratio_IC'])

        #periodic condition
        for grain in dict_ic['L_g_tempo']:
            #left wall
            if grain.center[0] < dict_sample['x_box_min'] :
                grain.center = grain.center.copy() + np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0])
                for i in range(len(grain.l_border)):
                    grain.l_border[i] = grain.l_border[i].copy() + np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0])
                    grain.l_border_x[i] = grain.l_border_x[i].copy() + dict_sample['x_box_max'] - dict_sample['x_box_min']
                #contact gimage needed to be convert into gg
                Create_IC.convert_gimage_into_gg(grain, dict_ic, dict_material)
                #contact gg needed to be convert into gimage
                Create_IC.convert_gg_into_gimage(grain, dict_ic, dict_material)

            #right wall
            elif grain.center[0] > dict_sample['x_box_max'] :
                grain.center = grain.center.copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                for i in range(len(grain.l_border)):
                    grain.l_border[i] = grain.l_border[i].copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                    grain.l_border_x[i] = grain.l_border_x[i].copy() + dict_sample['x_box_min'] - dict_sample['x_box_max']
                #contact gimage needed to be convert into gg
                Create_IC.convert_gimage_into_gg(grain, dict_ic, dict_material)
                #contact gg needed to be convert into gimage
                Create_IC.convert_gg_into_gimage(grain, dict_ic, dict_material)

        #Control the top group to have the pressure target
        dy_top, Fv = Control_Top_NR(dict_sollicitation['Vertical_Confinement_Force'],dict_ic['L_contact_gg']+dict_ic['L_contact_gimage'],dict_ic['L_g_tempo'])
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy_top
        for grain in dict_ic['L_g_tempo'] :
            if grain.group == 'Top':
                grain.move_as_a_group(np.array([0, dy_top]), dict_ic['dt_DEM_IC'])

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(dict_sample['y_box_max'])

        if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] ==0:
            if dict_sollicitation['gravity'] > 0 :
                print('i_DEM',dict_ic['i_DEM_IC'],'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            else :
                print('i_DEM',dict_ic['i_DEM_IC'],'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            if dict_ic['Debug_DEM'] :
                Plot_Config_Loaded(dict_ic,dict_ic['i_DEM_IC'])

        #Check stop conditions for DEM
        if dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_IC'] + i_DEM_0:
             DEM_loop_statut = False
        if dict_sollicitation['gravity'] > 0:
            if Ecin < Ecin_stop and F < Force_stop and (0.95*dict_sollicitation['Vertical_Confinement_Force']<Fv and Fv<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                  DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_IC']*0.1 + i_DEM_0 and (0.95*dict_sollicitation['Vertical_Confinement_Force']<Fv and Fv<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False

    #Update dict
    dict_ic['L_L_g_tempo'].append(dict_ic['L_g_tempo'].copy())

#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    """
    Compute total kinetic energy.

        Input :
            a list of temporary grains (a list)
        Output :
            the total kinetic energy (a float)
    """
    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def F_total(L_g):
    """
    Compute total force applied on grains in the sample.

        Input :
            a list of temporary grains (a list)
        Output :
            the total force applied (a float)
    """
    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy])
    return F

#-------------------------------------------------------------------------------

def Control_Top_NR(Force_target,L_contact_gg,L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a confinement value (a float)
            a list of contact grain - grain and grain - image (a list)
            a list of temporary grain (a list)
        Output :
            the displacement of the top group (a float)
            a force applied on the top group before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gg:
        #compute force applied, save contact overlap and spring
        if (contact.g1.group == 'Current' and contact.g2.group == 'Top') :
            F = F - contact.F_2_1_n * contact.pc_normal[1] #- because F_2_1_n < 0
            overlap_L.append(contact.overlap)
            k_L.append(contact.k*contact.pc_normal[1])
        elif (contact.g1.group == 'Top' and contact.g2.group == 'Current') :
            F = F + contact.F_2_1_n * contact.pc_normal[1] #+ because F_2_1_n < 0 and pc_normal from top to current
            overlap_L.append(contact.overlap)
            k_L.append(-contact.k*contact.pc_normal[1]) #because pc_normal from top to current

    if overlap_L != []:
        i_NR = 0
        dy_top = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_ymax_f(dy_top,overlap_L,k_L,Force_target) and error_on_ymax_f(dy_top,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy_top = dy_top - error_on_ymax_f(dy_top,overlap_L,k_L,Force_target)/error_on_ymax_df(dy_top,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy_top,overlap_L,k_L,Force_target) and error_on_ymax_f(dy_top,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False

    else :
        #if there is no contact with the upper wall, the wall is reset
        raise ValueError('Not available for the moment')
        y_max = Reset_y_max(L_g,Force_target)

    return dy_top, F

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
            a confinement force (a float)
        Output :
            the difference between the force applied and the confinement (a float)
    """
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    """
    Compute the derivative function df to control the upper wall (error_on_ymax_f()).

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_y_max(L_g,Force):
    """
    The Top group is located as a single contact verify the target value.

        Input :
            the list of temporary grains (a list)
            the confinement force (a float)
        Output :
            the upper wall position (a float)
    """
    y_max = None
    id_grain_max = None
    for id_grain in range(len(L_g)):
        grain = L_g[id_grain]
        y_max_grain = grain.center[1] + grain.radius

        if y_max != None and y_max_grain > y_max:
            y_max = y_max_grain
            id_grain_max = id_grain
        elif y_max == None:
            y_max = y_max_grain
            id_grain_max = id_grain

    factor = 5
    k = factor*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].radius)
    y_max = y_max - (Force/k)**(2/3)

    return y_max

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(dict_ic,i):
    """
    Plot loaded configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    L_color_group = ['k','r','b']
    L_group = ['Current', 'Bottom', 'Top']
    plt.figure(1,figsize=(16,9))
    for grain in dict_ic['L_g_tempo']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x,grain.l_border_y,L_color_group[i_group])
    for grain in dict_ic['L_g_image']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group]:
                plt.plot(grain.l_border_x,grain.l_border_y,'-.', color = L_color_group[i_group])
    plt.axis('equal')
    plt.savefig('Debug/Configuration/Shear/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------
