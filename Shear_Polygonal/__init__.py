# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#Own
import Create_IC_Polygonal
from Create_IC_Polygonal.Grain_ic_polygonal import Grain_Tempo_Polygonal, Grain_Image_Polygonal
import Create_IC_Polygonal.Contact_gg_ic_polygonal
import Create_IC_Polygonal.Contact_gimage_ic_polygonal

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def DEM_shear_load(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Loading the granular system with vertical load and shear.

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
    Shear_strain = 0
    #compute the sample height
    min_value = min(dict_ic['L_g_tempo'][0].l_border_y)
    max_value = max(dict_ic['L_g_tempo'][0].l_border_y)
    for grain in dict_ic['L_g_tempo']:
        if min(grain.l_border_y) < min_value :
            min_value = min(grain.l_border_y)
        if max(grain.l_border_y) > max_value :
            max_value = max(grain.l_border_y)
    Sample_height = max_value - min_value
    #track total displacement of grains
    for grain in dict_ic['L_g_tempo'] :
        grain.track_u = True
        grain.total_ux = 0
        grain.total_uy = 0
    DEM_loop_statut = True
    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gimage'] = []
    dict_ic['L_contact_ij_gimage'] = []
    dict_ic['id_contact'] = 0
    dict_ic['L_g_image'] = []
    dict_ic['L_i_image'] = []

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
                    dict_ic['L_g_image'].append(Grain_Image_Polygonal(grain, 'left'))
                    dict_ic['L_i_image'].append(grain.id)
            #right wall
            elif (dict_sample['x_box_max'] - grain.center[0]) < dict_algorithm['d_to_image'] :
                if grain.id in dict_ic['L_i_image'] : #image exists
                    image = dict_ic['L_g_image'][dict_ic['L_i_image'].index(grain.id)]
                    if image.position == 'left' :
                        image.position = 'right'
                else : #image does not exist
                    dict_ic['L_g_image'].append(Grain_Image_Polygonal(grain, 'right'))
                    dict_ic['L_i_image'].append(grain.id)
            #center
            else :
                if grain.id in dict_ic['L_i_image'] : #image exists
                    i_toremove = dict_ic['L_i_image'].index(grain.id)
                    dict_ic['L_g_image'].pop(i_toremove)
                    dict_ic['L_i_image'].pop(i_toremove)
                    L_i_toremove = []
                    for ij_gimage in dict_ic['L_contact_ij_gimage'] :
                        if grain.id == ij_gimage[1] :
                            L_i_toremove.append(dict_ic['L_contact_ij_gimage'].index(ij_gimage))
                    L_i_toremove.reverse()
                    for i_toremove in L_i_toremove:
                        dict_ic['L_contact_gimage'].pop(i_toremove)
                        dict_ic['L_contact_ij_gimage'].pop(i_toremove)
        #translate image
        for image in dict_ic['L_g_image']:
            if image.position == 'left' :
                image.translation(np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0]))
            elif image.position == 'right' :
                image.translation(np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0]))

        #Contact detection
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_ic['i_update_neighborhoods_com']  == 0:
            Create_IC_Polygonal.Contact_gg_ic_polygonal.Update_Neighborhoods(dict_ic)
            Create_IC_Polygonal.Contact_gimage_ic_polygonal.Update_Neighborhoods(dict_ic)
        Create_IC_Polygonal.Contact_gg_ic_polygonal.Grains_contact_Neighborhoods(dict_ic,dict_material)
        Create_IC_Polygonal.Contact_gimage_ic_polygonal.Grains_contact_Neighborhoods(dict_ic,dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in  dict_ic['L_contact']+dict_ic['L_contact_gimage']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_ic['dt_DEM_IC'])

        #Delete contacts gg and gimage with no overlap
        L_i_toremove = []
        for i_contact in range(len(dict_ic['L_contact'])):
            if dict_ic['L_contact'][i_contact].overlap_normal < 0:
                L_i_toremove.append(i_contact)
        L_i_toremove.reverse()
        for i_toremove in L_i_toremove:
            dict_ic['L_contact'].pop(i_toremove)
            dict_ic['L_contact_ij'].pop(i_toremove)
        L_i_toremove = []
        for i_contact in range(len(dict_ic['L_contact_gimage'])):
            if dict_ic['L_contact_gimage'][i_contact].overlap_normal < 0:
                L_i_toremove.append(i_contact)
        L_i_toremove.reverse()
        for i_toremove in L_i_toremove:
            dict_ic['L_contact_gimage'].pop(i_toremove)
            dict_ic['L_contact_ij_gimage'].pop(i_toremove)

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
                Create_IC_Polygonal.convert_gimage_into_gg(grain, dict_ic, dict_material)
                #contact gg needed to be convert into gimage
                Create_IC_Polygonal.convert_gg_into_gimage(grain, dict_ic, dict_material)
            #right wall
            elif grain.center[0] > dict_sample['x_box_max'] :
                grain.center = grain.center.copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                for i in range(len(grain.l_border)):
                    grain.l_border[i] = grain.l_border[i].copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                    grain.l_border_x[i] = grain.l_border_x[i].copy() + dict_sample['x_box_min'] - dict_sample['x_box_max']
                #contact gimage needed to be convert into gg
                Create_IC_Polygonal.convert_gimage_into_gg(grain, dict_ic, dict_material)
                #contact gg needed to be convert into gimage
                Create_IC_Polygonal.convert_gg_into_gimage(grain, dict_ic, dict_material)

        #Control the top group to have the pressure target
        dy_top, Fv = Control_Top_NR(dict_sollicitation['Vertical_Confinement_Force'],dict_ic['L_contact']+dict_ic['L_contact_gimage'],dict_ic['L_g_tempo'])
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy_top
        for grain in dict_ic['L_g_tempo'] :
            if grain.group == 'Top':
                grain.move_as_a_group(np.array([dict_sollicitation['Shear_velocity']*dict_ic['dt_DEM_IC'], dy_top]), dict_ic['dt_DEM_IC'])
        #Update shear strain
        Shear_strain = Shear_strain + dict_sollicitation['Shear_velocity']*dict_ic['dt_DEM_IC'] / Sample_height

        if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] == 0:
            print('i_DEM',dict_ic['i_DEM_IC'],'and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'% and Shear',int(100*Shear_strain/dict_sollicitation['Shear_strain_target']),'%')
            if dict_ic['Debug_DEM'] :
                Plot_Config_Loaded(dict_ic,dict_ic['i_DEM_IC'])

        #Check stop conditions for DEM
        if Shear_strain >= dict_sollicitation['Shear_strain_target'] :
             DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False

    #plot total displacement field
    Plot_total_U(dict_ic)
    
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
        #only the normal component is considered, no tangential one, no damping
        if (contact.g1.group == 'Current' and contact.g2.group == 'Top') :
            F = F - contact.F_2_1_n * contact.pc_normal[1] #- because F_2_1_n < 0
            overlap_L.append(contact.overlap_normal)
            k_L.append(contact.k*contact.pc_normal[1])
        elif (contact.g1.group == 'Top' and contact.g2.group == 'Current') :
            F = F + contact.F_2_1_n * contact.pc_normal[1] #+ because F_2_1_n < 0 and pc_normal from top to current
            overlap_L.append(contact.overlap_normal)
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
        dy_top = Reset_y_max(L_g)

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

def Reset_y_max(L_g):
    """
    The Top group is going down by a fraction of the mean radius.

    The confinement force is not verified.

        Input :
            the list of temporary grains (a list)
        Output :
            the upper wall position (a float)
    """
    R_mean = 0
    n_grain = 0
    for grain in L_g :
        if grain.group == 'Top' :
            R_mean = R_mean + grain.radius
            n_grain = n_grain + 1
    dy_top = - R_mean / n_grain / 20
    return dy_top

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(dict_ic,x_min,x_max,y_min,y_max,i):
    """
    Plot loaded configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    for grain in dict_ic['L_g_tempo']:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
    for grain in dict_ic['L_g_image']:
        plt.plot(grain.l_border_x,grain.l_border_y,'-.k')
        plt.plot(grain.center[0],grain.center[1],'xk')
    plt.plot([x_min,x_max],[y_min,y_min],'k')
    plt.plot([x_min,x_max],[y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/Configuration/Init/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_total_U(dict_ic):
    """
    Plot the map of the total displacement tracked.

        Input :
            a list of temporary grain (a list)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    L_color_group = ['k','r','b']
    L_group = ['Current', 'Bottom', 'Top']
    plt.figure(1,figsize=(16,9))
    x_L = []
    y_L = []
    u_L = []
    v_L = []
    for grain in dict_ic['L_g_tempo']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x,grain.l_border_y,L_color_group[i_group])
        if grain.group == 'Current' :
            x_L.append(grain.center[0])
            y_L.append(grain.center[1])
            u_L.append(grain.total_ux)
            v_L.append(grain.total_uy)
    plt.quiver(x_L,y_L,u_L,v_L)
    plt.axis('equal')
    plt.savefig('Debug/Configuration/Total_U_Sheared.png')
    plt.close(1)

#-------------------------------------------------------------------------------
