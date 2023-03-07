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
import Create_IC.Grain_ic
import Create_IC.Contact_gg_ic
import Create_IC.Contact_gimage_ic
import Create_IC.Contact_gw_ic

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Create an initial condition

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    #define the y_max for the grains generation
    radius_mean = 0
    for i in range(len(dict_geometry['L_R'])):
        radius_mean = radius_mean + dict_geometry['L_R'][i]*dict_geometry['L_percentage_R'][i]
    dy_creation = dict_geometry['N_grain']/dict_ic['n_generation']*dict_ic['factor_ymax_box']*(2*radius_mean)**2/(dict_sample['x_box_max']-dict_sample['x_box_min'])
    dict_sample['dy_creation'] = dy_creation

    #plan the grains generation
    L_n_grain_radius_ite = []
    L_n_grain_radius_final = []
    L_n_grain_radius_done = []
    for percentage in dict_geometry['L_percentage_R']:
        L_n_grain_radius_ite.append(round(dict_geometry['N_grain']*percentage/dict_ic['n_generation'],0))
        L_n_grain_radius_final.append(round(dict_geometry['N_grain']*percentage,0))
        L_n_grain_radius_done.append(0)
    dict_ic['L_n_grain_radius_ite'] = L_n_grain_radius_ite
    dict_ic['L_n_grain_radius_final'] = L_n_grain_radius_final
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    dict_ic['i_DEM_IC']  = 0
    dict_ic['L_L_g_tempo'] = []
    dict_sample['y_box_min_ic'] = dict_sample['y_box_min']
    dict_ic['last_id'] = 0

    #---------------------------------------------------------------------------

    for i_generation in range(1,dict_ic['n_generation']+1) :

        print(f'Generation {i_generation} of grains')

        #add elements in dicts
        dict_ic['L_g_tempo'] = []
        dict_ic['i_generation'] = i_generation

        #create not overlaping grains
        Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, simulation_report)

        #DEM to find the steady-state configuration after loading
        #find the maximum y (center+radius)
        y_max = dict_sample['y_box_min_ic']
        for grain in dict_ic['L_g_tempo']:
            if grain.center[1] + grain.radius > y_max:
                y_max = grain.center[1] + grain.radius

        #add element in dict
        dict_sample['y_box_max'] = y_max

        DEM_loading(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

        #update element in dict
        dict_sample['y_box_min_ic'] = dict_sample['y_box_max']

    #---------------------------------------------------------------------------

    print('Combine generations of grains')

    dict_ic['i_generation'] = dict_ic['n_generation']+1

    dict_ic['L_g_tempo'] = []
    for L_g_tempo in dict_ic['L_L_g_tempo']:
        for g_tempo in L_g_tempo:
            dict_ic['L_g_tempo'].append(g_tempo)

    DEM_loading(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    simulation_report.write_and_print(str(len(dict_ic['L_g_tempo']))+' / '+str(dict_geometry['N_grain'])+' grains have been created\n','\n'+str(len(dict_ic['L_g_tempo']))+' / '+str(dict_geometry['N_grain'])+' grains have been created\n')

#-------------------------------------------------------------------------------

def DEM_loading(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Loading the granular system.

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simultion report (a report)
        Output :
            Nothing, but initial condition dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    if dict_ic['i_generation'] == dict_ic['n_generation']+1 :
        i_update_neighborhoods = dict_ic['i_update_neighborhoods_com']
        y_min = dict_sample['y_box_min']
    else :
        i_update_neighborhoods = dict_ic['i_update_neighborhoods_gen']
        y_min = dict_sample['y_box_min_ic']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    i_DEM_0 = dict_ic['i_DEM_IC']
    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gimage'] = []
    dict_ic['L_contact_ij_gimage'] = []
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []
    dict_ic['id_contact'] = 0
    dict_ic['L_g_image'] = []
    dict_ic['L_i_image'] = []

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    Fv_tracker = []
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
                    dict_ic['L_g_image'].append(Grain_ic.Grain_Image(grain, 'left'))
                    dict_ic['L_i_image'].append(grain.id)
            #right wall
            elif (dict_sample['x_box_max'] - grain.center[0]) < dict_algorithm['d_to_image'] :
                if grain.id in dict_ic['L_i_image'] : #image exists
                    image = dict_ic['L_g_image'][dict_ic['L_i_image'].index(grain.id)]
                    if image.position == 'left' :
                        image.position = 'right'
                else : #image does not exist
                    dict_ic['L_g_image'].append(Grain_ic.Grain_Image(grain, 'right'))
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
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % i_update_neighborhoods  == 0:
            Contact_gg_ic.Update_Neighborhoods(dict_ic)
            Contact_gimage_ic.Update_Neighborhoods(dict_ic)
        Contact_gg_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)
        Contact_gimage_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)

        # Detection of contacts between grain and walls
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % i_update_neighborhoods  == 0:
            wall_neighborhood = Contact_gw_ic.Update_wall_Neighborhoods(dict_ic['L_g_tempo'],dict_ic['factor_neighborhood_IC'],dict_sample['x_box_min'],dict_sample['x_box_max'],y_min,dict_sample['y_box_max'])
        Contact_gw_ic.Grains_Polyhedral_Wall_contact_Neighborhood(wall_neighborhood,dict_sample['x_box_min'],dict_sample['x_box_max'],y_min,dict_sample['y_box_max'], dict_ic, dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in  dict_ic['L_contact']+dict_ic['L_contact_gimage']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_ic['dt_DEM_IC'])

        #Move grains
        for grain in dict_ic['L_g_tempo']:
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
                convert_gimage_into_gg(grain, dict_ic, dict_material)
                #contact gg needed to be convert into gimage
                convert_gg_into_gimage(grain, dict_ic, dict_material)
            #right wall
            elif grain.center[0] > dict_sample['x_box_max'] :
                grain.center = grain.center.copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                for i in range(len(grain.l_border)):
                    grain.l_border[i] = grain.l_border[i].copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                    grain.l_border_x[i] = grain.l_border_x[i].copy() + dict_sample['x_box_min'] - dict_sample['x_box_max']
                #contact gimage needed to be convert into gg
                convert_gimage_into_gg(grain, dict_ic, dict_material)
                #contact gg needed to be convert into gimage
                convert_gg_into_gimage(grain, dict_ic, dict_material)

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(dict_ic['L_g_tempo'])):
            if dict_ic['L_g_tempo'][id_grain].center[1] < y_min :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] > dict_sample['y_box_max'] :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box')
            dict_ic['L_g_tempo'].pop(id_grain)

        #Control the y_max to have the pressure target
        dict_sample['y_box_max'], Fv = Control_y_max_NR(dict_sample['y_box_max'],dict_sollicitation['Vertical_Confinement_Force'],dict_ic['L_contact_gw'],dict_ic['L_g_tempo'])

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(dict_sample['y_box_max'])
        Fv_tracker.append(Fv)

        if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] ==0:
            if dict_sollicitation['gravity'] > 0 :
                print('i_DEM',dict_ic['i_DEM_IC'],'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            else :
                print('i_DEM',dict_ic['i_DEM_IC'],'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            if dict_ic['Debug_DEM'] :
                Plot_Config_Loaded(dict_ic,dict_sample['x_box_min'],dict_sample['x_box_max'],y_min,dict_sample['y_box_max'],dict_ic['i_DEM_IC'])

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

    #plot trackers
    if dict_ic['Debug_DEM'] and dict_ic['i_generation'] == dict_ic['n_generation']+1 :
        fig, ((ax1, ax2)) = plt.subplots(1,2, figsize=(16,9),num=1)

        ax1.set_title('Total kinetic energy (e-12 J)')
        ax1.plot(Ecin_tracker)
        ax1.plot([0, len(Ecin_tracker)-1],[Ecin_stop, Ecin_stop],'r')

        ax2.set_title('About the upper plate')
        ax2.plot(Ymax_tracker, color = 'blue')
        ax2.set_ylabel('Coordinate y (µm)', color = 'blue')
        ax2.tick_params(axis ='y', labelcolor = 'blue')
        ax2a = ax2.twinx()
        ax2a.plot(range(50,len(Fv_tracker)),Fv_tracker[50:], color = 'orange')
        ax2a.plot([50, len(Fv_tracker)-1],[dict_sollicitation['Vertical_Confinement_Force'], dict_sollicitation['Vertical_Confinement_Force']], color = 'red')
        ax2a.set_ylabel('Force applied (µN)', color = 'orange')
        ax2a.tick_params(axis ='y', labelcolor = 'orange')

        fig.savefig('Debug/Init_disks_trackers.png')
        plt.close(1)

#-------------------------------------------------------------------------------

def Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, simulation_report):
    """
    Generate the grains.

    A position is tried. This position must not create overlap with already creared temporary grain. If there is no overlap, a temporary grai nis created.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a sample dictionnary (a dict)
            a material dictionnary (a dict)
            a generation id (a int)
            a simulation report (a report)
        Output :
            Nothing, but initial configuration dictionnary is updated
    """
    #Parameters for the method
    n_not_created = 0
    L_n_grain_radius = []
    if dict_ic['i_generation'] < dict_ic['n_generation']+1:
        for i in range(len(dict_ic['L_n_grain_radius_ite'])):
            L_n_grain_radius.append(dict_ic['L_n_grain_radius_ite'][i]*dict_ic['i_generation'])
    else :
        L_n_grain_radius = dict_ic['L_n_grain_radius_final']

    for i in range(len(dict_geometry['L_R'])):
        radius = dict_geometry['L_R'][i]
        n_grain = L_n_grain_radius[i]
        n_grain_done = dict_ic['L_n_grain_radius_done'][i]
        last_id_grain_created = dict_ic['last_id']
        for id_grain in range(last_id_grain_created, int(last_id_grain_created + n_grain - n_grain_done)):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < dict_ic['N_test_max']:
                i_test = i_test + 1
                center = np.array([random.uniform(dict_sample['x_box_min']+1.1*radius,dict_sample['x_box_max']-1.1*radius),random.uniform(dict_sample['y_box_min_ic']+1.1*radius,dict_sample['y_box_min_ic'] + dict_sample['dy_creation'])])
                g_tempo = Grain_ic.Grain_Tempo(id_grain-n_not_created, center, radius, dict_material)
                grain_created = True
                for grain in dict_ic['L_g_tempo']:
                    if Contact_gg_ic.Grains_contact_f(g_tempo,grain):
                        grain_created = False
            if i_test == dict_ic['N_test_max'] and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write_and_print('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n','Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries')
            else :
                dict_ic['L_g_tempo'].append(g_tempo)
                dict_ic['L_n_grain_radius_done'][i] = dict_ic['L_n_grain_radius_done'][i] + 1
                dict_ic['last_id'] = dict_ic['last_id'] + 1

#-------------------------------------------------------------------------------

def convert_gimage_into_gg(grain, dict_ic, dict_material):
    """
    Convert a contact grain-image in a contact grain-grain.

        Input :
            a grain (a grain_tempo)
            an initial dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial dictionnary is updated
    """
    L_i_contact_to_delete = []
    for ij_gimage in dict_ic['L_contact_ij_gimage'] :
        i_contact_gimage = dict_ic['L_contact_ij_gimage'].index(ij_gimage)
        if ij_gimage[0] == grain.id or ij_gimage[1] == grain.id:
            if ij_gimage[0] > ij_gimage[1] :
                ij_gg = (ij_gimage[1], ij_gimage[0])
            else :
                ij_gg = ij_gimage
            i_grain = 0
            grain_i = dict_ic['L_g_tempo'][i_grain]
            while not grain_i.id == ij_gg[0] :
                i_grain = i_grain + 1
                grain_i = dict_ic['L_g_tempo'][i_grain]
            j_grain = i_grain + 1
            grain_j = dict_ic['L_g_tempo'][j_grain]
            while not grain_j.id == ij_gg[1] :
                j_grain = j_grain + 1
                grain_j = dict_ic['L_g_tempo'][j_grain]
            if ij_gg not in dict_ic['L_contact_ij'] :
                #creation of contact
                dict_ic['L_contact_ij'].append(ij_gg)
                dict_ic['L_contact'].append(Create_IC.Contact_gg_ic.Contact_Tempo(dict_ic['id_contact'], grain_i, grain_j, dict_material))
                dict_ic['id_contact'] = dict_ic['id_contact'] + 1
                #transmit data
                dict_ic['L_contact'][-1].convert_gimage_in_gg(dict_ic['L_contact_gimage'][i_contact_gimage])
                #update neighborhood
                grain_i.neighbourood.append(grain_j)
            L_i_contact_to_delete.append(i_contact_gimage)
    #delete previous contact gimage
    L_i_contact_to_delete.reverse()
    for i_contact_to_delete in L_i_contact_to_delete :
        dict_ic['L_contact_gimage'].pop(i_contact_to_delete)
        dict_ic['L_contact_ij_gimage'].pop(i_contact_to_delete)

#-------------------------------------------------------------------------------

def convert_gg_into_gimage(grain, dict_ic, dict_material):
    """
    Convert a contact grain-grain in a contact grain-gimage.

        Input :
            a grain (a grain_tempo)
            an initial dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial dictionnary is updated
    """
    L_i_contact_to_delete = []
    for ij_gg in dict_ic['L_contact_ij'] :
        i_contact_ij = dict_ic['L_contact_ij'].index(ij_gg)
        if ij_gg[0] == grain.id or ij_gg[1] == grain.id:
            ij_gimage = ij_gg
            if ij_gimage[1] in dict_ic['L_i_image'] :
                #contact gimage 1
                i_grain = 0
                grain = dict_ic['L_g_tempo'][i_grain]
                while not grain.id == ij_gimage[0] :
                    i_grain = i_grain + 1
                    grain = dict_ic['L_g_tempo'][i_grain]
                i_image = 0
                image = dict_ic['L_g_image'][i_image]
                while not image.id == ij_gimage[1] :
                    i_image = i_image + 1
                    image = dict_ic['L_g_image'][i_image]
                #creation of contact
                dict_ic['L_contact_ij_gimage'].append(ij_gimage)
                dict_ic['L_contact_gimage'].append(Create_IC.Contact_gimage_ic.Contact_Image(dict_ic['id_contact'], grain, image, dict_material))
                dict_ic['id_contact'] = dict_ic['id_contact'] + 1
                #transmit data
                dict_ic['L_contact_gimage'][-1].convert_gimage_in_gg(dict_ic['L_contact'][i_contact_ij])
                #update neighborhood
                grain.neighbourood_image.append(image)

            #contact gimage 2
            ij_gimage = (ij_gg[1], ij_gg[0])
            if ij_gimage[1] in dict_ic['L_i_image'] :
                #contact gimage 1
                i_grain = 0
                grain = dict_ic['L_g_tempo'][i_grain]
                while not grain.id == ij_gimage[0] :
                    i_grain = i_grain + 1
                    grain = dict_ic['L_g_tempo'][i_grain]
                i_image = 0
                image = dict_ic['L_g_image'][i_image]
                while not image.id == ij_gimage[1] :
                    i_image = i_image + 1
                    image = dict_ic['L_g_image'][i_image]
                #creation of contact
                dict_ic['L_contact_ij_gimage'].append(ij_gimage)
                dict_ic['L_contact_gimage'].append(Create_IC.Contact_gimage_ic.Contact_Image(dict_ic['id_contact'], grain, image, dict_material))
                dict_ic['id_contact'] = dict_ic['id_contact'] + 1
                #transmit data
                dict_ic['L_contact_gimage'][-1].convert_gimage_in_gg(dict_ic['L_contact'][i_contact_ij])
                #update neighborhood
                grain.neighbourood_image.append(image)

            L_i_contact_to_delete.append(i_contact_ij)
    #delete previous contact gimage
    L_i_contact_to_delete.reverse()
    for i_contact_to_delete in L_i_contact_to_delete :
        dict_ic['L_contact'].pop(i_contact_to_delete)
        dict_ic['L_contact_ij'].pop(i_contact_to_delete)

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

def Control_y_max_NR(y_max,Force_target,L_contact_gw,L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a coordinate of the upper wall (a float)
            a confinement value (a float)
            a list of contact grain - wall (a list)
            a list of temporary grain (a list)
        Output :
            the coordinate of the upper wall (a float)
            a force applied on the upper wall before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dy = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy = dy - error_on_ymax_f(dy,overlap_L,k_L,Force_target)/error_on_ymax_df(dy,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        y_max = y_max + dy

    else :
        #if there is no contact with the upper wall, the wall is reset
        y_max = Reset_y_max(L_g,Force_target)

    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            #reactualisation
            contact.limit = y_max

    return y_max, F

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
    The upper wall is located as a single contact verify the target value.

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
    plt.savefig('Debug/Configuration/Init_disks/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------
