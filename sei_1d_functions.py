# -*- coding: utf-8 -*-

import numpy as np
import cantera as ct

"""----------Residual functions for IDA solver----------"""


"""------------Spatially Heterogeneous Model-----------"""
def residual(t, SV, SV_dot):
    from sei_1d_init import objs, params, voltage_lookup, SVptr

    # Initialize residual equal to all zeros:
    res = SV_dot

    # Read out cantera objects:
    WE = objs['WE']
    sei = objs['SEI']
    elyte = objs['elyte']
    sei_elyte = objs['SEI_elyte']
    sei_conductor = objs['conductor']
    WE_sei = objs['WE_SEI']

    phi_WE = np.interp(t,voltage_lookup['time'],voltage_lookup['voltage'])
    
    i_sei = np.zeros(params['Ny']+1,)
    i_Far = np.zeros(params['Ny'],)

    # Set anode and anode/electrolyte interface potentials
    WE.electric_potential = phi_WE
    
    # Start in the volume adjacent to the working electrode:
    j=0

    # SEI electric potential:
    phi_sei_loc = SV[SVptr['phi sei'][j]]
    sei.electric_potential = phi_sei_loc
    sei_conductor.electric_potential = phi_sei_loc

    # SEI volume fraction:
    eps_sei_loc = SV[SVptr['eps sei'][j]]

    # The current in the SEI entering this voluem is that produced by
    #   charge-transfer reactions at the WE-SEI interface:
    i_sei[j] = eps_sei_loc*WE_sei.get_net_production_rates(WE)*ct.faraday

    # Surface area equals the geometric area, for this interface:
    rxn_scale = 1

    # Loop through the remaining volumes (except for the very last one):
    for j in range(params['Ny']-1):

        # Read out local SEI composition and set Cantera object:
        Ck_sei_loc = SV[SVptr['Ck sei'][j]]
        Xk_sei_loc = Ck_sei_loc/sum(Ck_sei_loc)
        sei.X = Xk_sei_loc

        # SEI electric potential:
        phi_sei_loc =  SV[SVptr['phi sei'][j]]
        sei.electric_potential = phi_sei_loc
        sei_conductor.electric_potential = phi_sei_loc

        # Electrolyte electric potential assumed to be zero:
        elyte.electric_potential = 0.

        # sei-electrolyte area per unit volume.  This is scaled by 
        #   (1 - eps_sei)*eps_sei so that available area goes to zero
        #   as the sei volume fraction approaches either zero or one:
        sei_APV = 4.*eps_sei_loc*(1.0-eps_sei_loc)/params['d_sei']

        # Production rates from chemical reactions at sei-electrolyte interface:
        Rates_sei_elyte = sei_elyte.get_net_production_rates(sei)*sei_APV

        # Calculate residual for chemical molar concentrations:
        dSVdt_ck_sei = rxn_scale*Rates_sei_elyte
        res[SVptr['Ck sei'][j]] = SV_dot[SVptr['Ck sei'][j]] - dSVdt_ck_sei

        # Calculate residual for sei volume fraction:
        dSVdt_eps_sei = np.dot(dSVdt_ck_sei, params['vol_k sei'])
        res[SVptr['eps sei'][j]] = SV_dot[SVptr['eps sei'][j]] - dSVdt_eps_sei

        # Calculate faradaic current density due to charge transfer at SEI-elyte
        #   interface, in A/m2 total.
        Rates_sei_conductor = sei_elyte.get_net_production_rates(sei_conductor)
        i_Far[j] = Rates_sei_conductor*sei_APV/params['dyInv']

        # sei electric potential at next volume:
        phi_sei_next = SV[SVptr['phi sei'][j+1]]

        # Sei volume fraction at interface between volumes:
        eps_sei_next = SV[SVptr['eps sei'][j+1]]
        eps_sei_int = 0.5*(eps_sei_loc + eps_sei_next)
    
        # Current = (Conductivity)*(volume fraction)*(-grad(Phi))
        dPhi = phi_sei_loc - phi_sei_next
        sigma_sei = np.dot(params['sigma sei'],sei.X)
        i_sei[j+1] = eps_sei_int*sigma_sei*dPhi

        # Rates for next node are scaled by sei volume fraction in this node:
        rxn_scale = 0.5*(eps_sei_loc+eps_sei_next)
        eps_sei_loc = eps_sei_next

    # Repeat calculations for final node, where the boundary condition is
    #   that i_sei = 0 at the interface with the electrolyte:
    j = int(params['Ny']-1)
    Ck_sei_loc = SV[SVptr['Ck sei'][j]]
    Xk_sei_loc = Ck_sei_loc/sum(Ck_sei_loc)
    sei.X = Xk_sei_loc

    elyte.electric_potential = 0.

    phi_sei_loc =  SV[SVptr['phi sei'][j]]
    sei.electric_potential = phi_sei_loc
    sei_conductor.electric_potential = phi_sei_loc

    # SEI surface Area Per unit Volume (APV)
    sei_APV = 4.*eps_sei_loc*(1-eps_sei_loc)**2/params['d_sei']
    Rates_sei_elyte = sei_elyte.get_net_production_rates(sei)*sei_APV
    dSVdt_ck_sei = rxn_scale*Rates_sei_elyte
    res[SVptr['Ck sei'][j]] = SV_dot[SVptr['Ck sei'][j]] - dSVdt_ck_sei

    # Calculate faradaic current density due to charge transfer at SEI-elyte
    #   interface, in A/m2 total.
    Rates_sei_conductor = sei_elyte.get_net_production_rates(sei_conductor)
    i_Far[j] = Rates_sei_conductor*sei_APV/params['dyInv']

    dSVdt_eps_sei = np.dot(dSVdt_ck_sei, params['vol_k sei'])
    res[SVptr['eps sei'][j]] = SV_dot[SVptr['eps sei'][j]] - dSVdt_eps_sei

    i_dl = i_Far - i_sei[:-1] + i_sei[1:]
    dSVdt_phi_dl = -i_dl/params['C_dl WE_sei']
    res[SVptr['phi sei']] = SV_dot[SVptr['phi sei']] - dSVdt_phi_dl

    return res


    

"""------------Spatially Homogeneous Model-----------"""
def residual_homogeneous(t, SV, SV_dot):
    from sei_1d_init import objs, params, voltage_lookup, SVptr

    # Initialize residual equal to all zeros:
    res = SV_dot

    # Read out cantera objects:
    WE = objs['WE']
    sei = objs['SEI']
    elyte = objs['elyte']
    sei_elyte = objs['SEI_elyte']
    sei_conductor = objs['conductor']
    WE_sei = objs['WE_SEI']

    phi_WE = np.interp(t,voltage_lookup['time'],voltage_lookup['voltage'])
    
    # Set anode and anode/electrolyte interface potentials
    WE.electric_potential = phi_WE

    # SEI electric potential:
    phi_sei_WE = phi_WE + SV[SVptr['WE_dl']]
    sei.electric_potential = phi_sei_WE
    sei_conductor.electric_potential = phi_sei_WE

    # The current in the SEI entering this voluem is that produced by
    #   charge-transfer reactions at the WE-SEI interface:
    i_WE = WE_sei.get_net_production_rates(WE)*ct.faraday


    phi_sei_elyte = SV[SVptr['Elyte_dl']]
    sei.electric_potential = phi_sei_elyte
    sei_conductor.electric_potential = phi_sei_elyte
    
    Xk_sei = SV[SVptr['Xk_sei']]
    sei.X = Xk_sei

    

    i_elyte = sei_elyte.get_net_production_rates(sei_conductor)*ct.faraday
    sdot_sei = sei_elyte.get_net_production_rates(sei)

    # Repeat calculations for final node, where the boundary condition is
    #   that i_sei = 0 at the interface with the electrolyte:
    j = int(params['Ny']-1)
    Ck_sei_loc = SV[SVptr['Ck sei'][j]]
    Xk_sei_loc = Ck_sei_loc/sum(Ck_sei_loc)
    sei.X = Xk_sei_loc

    elyte.electric_potential = 0.

    phi_sei_loc =  SV[SVptr['phi sei'][j]]
    sei.electric_potential = phi_sei_loc
    sei_conductor.electric_potential = phi_sei_loc

    # SEI surface Area Per unit Volume (APV)
    sei_APV = 4.*eps_sei_loc*(1-eps_sei_loc)**2/params['d_sei']
    Rates_sei_elyte = sei_elyte.get_net_production_rates(sei)*sei_APV
    dSVdt_ck_sei = rxn_scale*Rates_sei_elyte
    res[SVptr['Ck sei'][j]] = SV_dot[SVptr['Ck sei'][j]] - dSVdt_ck_sei

    # Calculate faradaic current density due to charge transfer at SEI-elyte
    #   interface, in A/m2 total.
    Rates_sei_conductor = sei_elyte.get_net_production_rates(sei_conductor)
    i_Far[j] = Rates_sei_conductor*sei_APV/params['dyInv']

    dSVdt_eps_sei = np.dot(dSVdt_ck_sei, params['vol_k sei'])
    res[SVptr['eps sei'][j]] = SV_dot[SVptr['eps sei'][j]] - dSVdt_eps_sei

    i_dl = i_Far - i_sei[:-1] + i_sei[1:]
    dSVdt_phi_dl = -i_dl/params['C_dl WE_sei']
    res[SVptr['phi sei']] = SV_dot[SVptr['phi sei']] - dSVdt_phi_dl

    return res