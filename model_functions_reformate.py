#! /usr/local/bin/python

#########################################################################
## import necessary modules and define global variables
import numpy
import pylab
import scipy.stats
from thermodynamic_variables import Sol_prod,equil_cont
from clim_model import clim_fun_CO2_only,clim_fun_lowCH4,clim_fun_highCH4 
#########################################################################

# The function 'Forward_Model' takes input parameters and returns time-evolution of the key carbon cycle variables from 0 to 4 Ga
# (and mass imbalance). The input parameters to Forward_Model are defined in main_code_parallel.py.
def Forward_Model(W,F_outgass,n,climp,tdep_weath,mod_sea,alt_frac,Mp_frac,lfrac,carb_exp,sed_thick,F_carbw,fpel,CWF,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas):

    # define global variables:
    global pK1,pK2,H_CO2,salt,Pore_mod,T_surface0,Ca_p,Ca_o,options_array
    global landf_vary,lum_vary,salt0,k_w,k_c,k_o,k_r,Mp,ppCO2_o,Mo,ko1,ko2,ALK_o,ALK_p
    # load options defined in main_code_parallel.py:
    options_array=numpy.load('options_array.npy')
    
    # Initial (modern) ocean pH and atmospheric pCO2 (bar) 
    pH_o=8.2 
    ppCO2_o=.000280 
   
    T=clim_fun_CO2_only(ppCO2_o,1.0) # find initial (modern) temperature
    # Calculate equilibrium constants for carbon chemistry:
    [pK1,pK2,H_CO2]=equil_cont(T)
    
    Mo=1.35e21 # Mass of ocean (kg)
    W=Mo/W # convert mixing time of ocean (yr) to mass flux (kg/yr)
    Mp=Mp_frac*Mo # Mass of water in pore space         
    
    ##############################################################
    ############ Compute other initial conditions ################
    # First, we use the assumed values for modern seafloor carbonate precipitation, outgassing, and seafloor dissolution to carbonate precipitation ratio
    # along with the carbon cycle steady state assumption to calculated the modern seafloor dissolution and continental silciate weathering
    partition=mod_sea/F_outgass
    F_diss0=alt_frac*partition*F_outgass #initial seafloor dissolution flux i.e. dissolution = alt_frac*precip_seafloor
    F_sil0=(1-partition)*F_outgass+(1-alt_frac)*partition*F_outgass #initial continental silicate weathering flux
    Fprecip_p=partition*F_outgass #initial pore space carbonate precipation flux (reduntant)
    
    #vary solar luminosity? y or n (don't modify these options - better to change elsewhere).
    lum_vary="y" 
    landf_vary="y"
    
    ### Initial conditions for modern atmosphere and ocean (equations S12 to S14)
    CO2_aq_o=H_CO2 * ppCO2_o #aqueous dissolved CO2
    hco3_o=(10**-pK1) * CO2_aq_o / (10**-pH_o) # aqueous bicarbonate
    co3_o=(10**-pK2) * hco3_o / (10**-pH_o) #aqueous carbonate
    DIC_o=co3_o+hco3_o+CO2_aq_o #dissolved inorganic carbon
    ALK_o=2*co3_o+hco3_o #carbonate alkalinity
    Ca_o=0.0100278 # initial (modern) calcium molality
    salt = ALK_o-2*Ca_o ## set salt to ensure ALK is consistent with Ca
    salt0=salt # initial (modern) salt
    
    T_surface0=clim_fun_CO2_only(ppCO2_o,1.0)
    interc=274.037-deep_grad*T_surface0
    buffer_T=deep_grad*T_surface0+interc #initial deep ocean temperature
    Pore_mod=9.0 #difference between pore space and deep ocean temperature for modern Earth
    
    ## Use these constants and steady state assumption to calculate remaining initial conditions:
    b1=2*F_diss0+W*ALK_o-2*W*DIC_o
    b2=2*F_sil0-W*ALK_o-2*F_outgass+2*W*DIC_o
    a1=W
    a2=-2*W
    a3=-W
    a4=2*W
    
    # Solving system of equations for pore space properties
    DIC_p=DIC_o-Fprecip_p/W
    ALK_p=(b1-a2*DIC_p)/a1
    Fprecip_o=(F_outgass+F_carbw)-(DIC_o-DIC_p)*W

    if DIC_p==0:
        print "WARNING: NO PHYSIAL STEADY STATE EXISTS! DIC is negative"
   
    #Solve quadratic ignoring H in ALK but including CO2aq
    [vv,xx]=pylab.roots([ALK_p/(10**-pK1*10**-pK2),(ALK_p-DIC_p)/(10**-pK2),ALK_p-2*DIC_p]) #equation S15
    H_p=numpy.max([vv,xx]) # take positive root
    pH_p=-pylab.log10(H_p)
    
    #Find remaining carbon chemistry from equilibrium conditions:
    CO2aq_p=DIC_p/((10**-pK1)*(10**-pK2)/H_p**2+10**-pK1/H_p+1)
    co3_p=ALK_p-DIC_p-H_p+CO2aq_p ## only true for high pH, okay for fitting I think (actually should include CO2aq to get exact).
    hco3_p=co3_p*10**-pH_p/(10**-pK2) #Modern bicarbonate alkalinity in pore space (not used)

    Ca_p = 0.5*(ALK_p - salt) #initial calcium molality in pore space
    omega_p=Ca_p *co3_p/Sol_prod(buffer_T+Pore_mod) # saturation state pore space
    omega_o=Ca_o *co3_o/Sol_prod(T_surface0) # saturation state ocean
    ## All initital conditions have been calculated ############################
    ############################################################################
    
    ############ Calculate proportionality constants ###########################
    # Given modern precipitation fluxes and saturation states, calculate proportionality constants
    k_c=Fprecip_p/(omega_p - 1)**n
    k_o=Fprecip_o/(omega_o - 1)**n

    ## Partition ocean carbonate sink into shelf precipitation and carbonate precipitation
    frac_pel=fpel # no pelagic carbonates in this version of the code (i.e. default is fpel=0)
    ko1=(1-frac_pel)*Fprecip_o/(omega_o - 1)**n 
    ko2=frac_pel*Fprecip_o/(omega_o**2.84)        
    
    ### remaining proportionality constants
    k_r=F_diss0/(2.88*10**-14*10**(-coef_for_diss*pH_p)*numpy.exp(-Ebas/(8.314*(buffer_T+Pore_mod)))) #for seafloor dissolution
    k_w=F_sil0 # for continental silicate weathering
        
    time=numpy.linspace(0,4e9,100) #define time array from 0 to 4 Ga in 100 time steps. Number of time steps can be increased for greater precision (takes longer)
    # run ode solver with initial conditions defined above, and inputs from main_code_parallel.py. The ode solver will output the dissolved inorganic carbon and
    # alkalinities of the ocean and pore space through time (all containted in 'out'). The diagnostic 'mes' is also returned
    [out,mes]=scipy.integrate.odeint(system_of_equations, [DIC_o+1.8e20/Mo*ppCO2_o,ALK_o,DIC_p,ALK_p], time, args=(W,F_outgass,n,climp,tdep_weath,mod_sea,alt_frac,Mp_frac,lfrac,carb_exp,sed_thick,F_carbw,CWF,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas),full_output=1)
 
    # Given alkalinities and dissolved inorganic carbon, calculate time-evolution for other carbon cycle variables
    # define empty arrays for variables of interest:

    pH_array_o=0*time #ocean pH
    CO2_array_o=0*time # atmospheric CO2
    pH_array_p=0*time # pore space pH
    CO2_array_p=0*time # pore space CO2
    Ca_array_o=0*time # ocean calcium molality
    Ca_array_p=0*time # pore space calcium molality
    CO3_array_o=0*time # ocean carbonate molality
    CO3_array_p=0*time # pore space carbonate molality
    HCO3_array_o=0*time # ocean bicarbonate molality
    HCO3_array_p=0*time # pore space bicarbonate molality
    omega_o=0*time # saturation state ocean
    omega_p=0*time # saturation state pore space
    Tsurf_array=0*time # Surface temperature
    Tdeep_array=0*time # Deep ocean temperature
    Fd_array=0*time # Seafloor dissolution
    Fs_array=0*time # Continental silciate weathering
    Precip_ocean_array=0*time # ocean carbonate precipitation
    Precip_pore_array=0*time # pore space carbonate precipitation
    volc_array=0*time # volcanic outgassing
    carbw_ar=0*time # Continental carbonate weathering
    co2_aq_o=0*time # ocean aqueous CO2
    co2_aq_p=0*time # pore space aqueous CO2
    spread_ar=0*time # spreading rate relative to modern
    T_pore_space_ar=0*time # pore space temperature
   
    for i in range(0,len(time)): 
        # using ode solution (alkalinities and DIC abundances) and input parameters, fill array pl1 with other carbon cycle evolution variables
        pl1=carbon_cycle([out[i,0],out[i,1],out[i,2],out[i,3]],time[i],W,F_outgass,n,climp,tdep_weath,mod_sea,alt_frac,Mp_frac,lfrac,carb_exp,sed_thick,F_carbw,CWF,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas) #ocean
        # Fill in arrays define above
        pH_array_o[i]=pl1[0]
        co2_aq_o[i]=pl1[3]
        CO2_array_o[i]=pl1[4]
        pH_array_p[i]=pl1[19]
        CO2_array_p[i]=pl1[20]
        Ca_array_o[i]=pl1[5]
        Ca_array_p[i]=pl1[21]
        CO3_array_o[i]=pl1[1]
        CO3_array_p[i]=pl1[16]
        HCO3_array_o[i]=pl1[2]
        HCO3_array_p[i]=pl1[17]
        co2_aq_p[i]=pl1[18]
        omega_o[i]=pl1[6]
        omega_p[i]=pl1[13]
        Tsurf_array[i]=pl1[7] 
        Tdeep_array[i]=pl1[8]
        Fd_array[i]=pl1[9]  
        Fs_array[i]=pl1[10]
        Precip_ocean_array[i]=pl1[11]
        Precip_pore_array[i]=pl1[12]
        carbw_ar[i]=pl1[15]
        volc_array[i]=pl1[14]
        spread_ar[i]=pl1[26]
        T_pore_space_ar[i]=pl1[27]
        
    max_el=numpy.size(time)-1 #number of elements in time array
    tstep=numpy.max(time)/max_el #time step (yrs) in time array

    # calculate mass imbalance
    flux_bal=numpy.sum(volc_array[2:]+carbw_ar[2:])*tstep-numpy.sum(Precip_ocean_array[2:]+Precip_pore_array[2:])*tstep # integral of carbon inputs minus outpus over all timesteps 
    res_change=1.8e20*(CO2_array_o[max_el]-CO2_array_o[2])+Mo*(HCO3_array_o[max_el]+CO3_array_o[max_el]+co2_aq_o[max_el]-co2_aq_o[2]-CO3_array_o[2]-HCO3_array_o[2])+Mp*(HCO3_array_p[max_el]+CO3_array_p[max_el]+co2_aq_p[max_el]-co2_aq_p[2]-CO3_array_p[2]-HCO3_array_p[2]) # final - initial carbon abundances (starting at second element to avoid transients)
    res_change2=(out[max_el,0]-out[2,0])*Mo+(out[max_el,2]-out[2,2])*Mp # alternative way of calculating final - initial carbon abundances
    
    imbalance=(flux_bal-res_change)/(out[max_el,0]*Mo) # difference between integrated flux balance and final-initial abundances (should be zero).

    # return selected carbon cycle outputs to main_code_parallel.py
    return [out[:,0],out[:,1],out[:,2],out[:,3],time,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array,spread_ar,volc_array,T_pore_space_ar],imbalance 
    
### The function system_of_equations takes as inputs the current state of the carbon cycle (alkalinities and DICs for the ocean and pore space) plus input parameters, 
### and calculates the time derivatives of the current state. It calls the function carbon_cycle to obtain steady state fluxes, which are necessary for calculating analytic derivatives.    
def system_of_equations(y,t0,W,F_outgass,n,climp,tdep_weath,mod_sea,alt_frac,Mp_frac,lfrac,carb_exp,sed_thick,F_carbw,CWF,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas):
   
    global pK1,pK2,H_CO2,salt,Pore_mod,T_surface0
    global landf_vary,lum_vary,salt0,k_w,k_c,k_o,k_r,Mp,ppCO2_o,Mo,ko1,ko2
    s=1.8e20/Mo #correction factor for mass balance (see manuscript) 
    
    # Call carbon cycle forward model to obtain current state of carbon cycle, namely ocean chemistry and carbon cycle fluxes, 
    # from parameters, time, and DIC and ALK for ocean and pore space.
    [EE_pH_o,EE_co3_o,EE_hco3_o, EE_co2aq_o,EE_ppCO2_o,EE_Ca_o,EE_omega_o,T_surface,buffer_T,F_diss,F_silic,Precip_ocean,Precip_pore,EE_omega_p,F_outg,carb_weath,EE_co3_p,EE_hco3_p, EE_co2aq_p,EE_pH_p,EE_ppCO2_p,EE_Ca_p,EE_H_o,EE_H_p,T_surface_diff,F_outg,spread_value,T_pore_space]=carbon_cycle([y[0],y[1],y[2],y[3]],t0,W,F_outgass,n,climp,tdep_weath,mod_sea,alt_frac,Mp_frac,lfrac,carb_exp,sed_thick,F_carbw,CWF,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas)
      
    DICo=y[0]-EE_ppCO2_o*s # Correct dissolved inorganic carbon content of ocean by subtracting atmospheric pCO2 from atmosphere-ocean reservoir (this ensures mass balance)
    
    # Calculate time derivatives from current state (equation S1).
    dy0_dt=F_outg/Mo+carb_weath/Mo-W*(DICo-y[2])/Mo-Precip_ocean/Mo # Time derivative of atmosphere-ocean carbon abundance
    dy1_dt=2*carb_weath/Mo-W*(y[1]-y[3])/Mo+2*F_silic/Mo-2*Precip_ocean/Mo # Time derivative of ocean alkalinity
    dy2_dt=W*(DICo-y[2])/Mp-Precip_pore/Mp # Time derivative of pore space carbon abundance
    dy3_dt=W*(y[1]-y[3])/Mp+2*F_diss/Mp-2*Precip_pore/Mp # Time derivative of pore space alkalinity
    return [dy0_dt,dy1_dt,dy2_dt,dy3_dt] # return derivatives           
    
## The function lf calculates land fraction relative to modern as a function of time, the Archean land fraction, and the timing of continental growth
## See equation S6 for details
def lf (t0,lfrac,growth_timing):
    LL=1/(1-lfrac)
    land_frac=numpy.max([0.0,1-1/(LL+numpy.exp(-10*(t0/1e9-growth_timing)))]) # max needed to prevent negative land fractions
    return land_frac

## The function carbon_cycle takes as inputs carbonate alkalinity and carbon abundance in the ocean and pore space, along with input parameters and time (in years), and calculates other steady state
## carbon cycle variables including complete ocean chemistry, surface temperature, and carbon cycle fluxes.
def carbon_cycle (y,t0,W,F_outgass,n,climp,tdep_weath,mod_sea,alt_frac,Mp_frac,lfrac,carb_exp,sed_thick,F_carbw,CWF,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas):

    ## Ca_cofactor is only important for the sensitivity tests where a wide range of Archean Ca abundances are imposed
    Ca_cofactor=new_add_Ca/1000.0*numpy.exp(t0/1e9-4) # perscribed ocean Ca molality evolution
    
    global pK1,pK2,H_CO2,salt,Pore_mod,T_surface0,ALK_o,ALK_p,Ca_p,Ca_o,options_array
    global landf_vary,lum_vary,salt0,k_w,k_c,k_o,k_r,Mp,ppCO2_o,Mo,ko1,ko2#,Mg_fun,Ca_Tfun,coef_for_diss
    
    #ocean variables
    s=1.8e20/Mo #correction factor for mass balance (see manuscript) 

    #### Temperature loop ####
    # This loop is necessary because the carbon chemistry equilibrium constants are temperatue dependent,
    # but surface temperature is controlled by atmospheric pCO2, which is in turn modulated by 
    # the carbon chemistry of the ocean. Therefore, several iterations are required to ensure the temperature
    # used to calculate equilibrium constants is the same as the surface climate. Sensitivity tests reveal
    # that three iterations are sufficient for accurate (within 1 degree) convergence.
    # Note that if Carbon_chem is set to 0 then the only a single iteration is executed and 
    # the carbon chemistry constants are assumed to be fixed. This is faster and only slightly less
    # accurate than the full calculation.

    ## initial guesses for surface temperature and pore space, to be used to determine equilibrium carbon chemistry constants:
    T_upd=273+18 # initial guess for surface temperature
    T_upd2=273+18 # initial guess for pore-space temperature
    
    T_surface=100.0 # Filler value (T_surface will be calcauled within the loop)
    
     
    # Check to see if Carbon_chem is 0 or 1
    if options_array[1]==0:
        max_temp_iter=1 # run single iteration
    else:
        max_temp_iter=3 # run 3 iterations
    
    jk=0   
    for jk in range(0,max_temp_iter):
        ## Calculate temperature-dependent carbon chemistry constants for the ocean:
        [pK1,pK2,H_CO2]=equil_cont(T_upd)

        # Calculate equilibrium carbon chemistry for the ocean:
        [c,d]=pylab.roots([y[1]/(10**-pK2*10**-pK1)*(1+s/H_CO2),(y[1]-y[0])/(10**-pK2),(y[1]-2*y[0])]) #equation S15
        EE_H_o=numpy.max([c ,d])## make sure get right root!
        EE_pH_o=-pylab.log10(EE_H_o) # equation S16
        EE_co3_o=y[1]/(2+EE_H_o/(10**-pK2))
        EE_hco3_o=y[1]-2*EE_co3_o
        EE_co2aq_o=( EE_hco3_o*EE_H_o/(10**-pK1) ) # equation S13
        EE_ppCO2_o = EE_co2aq_o /H_CO2 # equation S12
        EE_Ca_o = Ca_o+0.5*(y[1]-ALK_o)+Ca_cofactor #equation S17
        
        ## Calculate temperature-dependent carbon chemistry constants for the pore space:
        [pK1,pK2,H_CO2]=equil_cont(T_upd2)

        # Calculate equilibrium carbon chemistry for the pore space:
        [c,d]=pylab.roots([y[3]/(10**-pK2*10**-pK1),(y[3]-y[2])/(10**-pK2),(y[3]-2*y[2])]) # equation S15
        EE_H_p=numpy.max([c ,d])## make sure get right root!
        EE_pH_p=-pylab.log10(EE_H_p) # equation S16
        EE_co3_p=y[3]/(2+EE_H_p/(10**-pK2))
        EE_hco3_p=y[3]-2*EE_co3_p
        EE_co2aq_p=( EE_hco3_p*EE_H_p/(10**-pK1) ) #equation S14
        EE_ppCO2_p = EE_co2aq_p /H_CO2 #equation S12
        EE_Ca_p = Ca_p+0.5*(y[3]-ALK_p)+Ca_cofactor # equation S17
    
        # Calculate biological modification of continental weathering using equation S7:
        fb_EE=1/(1-CWF)
        biology_modifier=1-1/(fb_EE+numpy.exp(-10*(t0-0.6e9)/1e9))
        #biology_modifier=1.0

        # Calculate internal heatflow, outgassing, and spreading rate using equations S8-10:
        Q=(1-t0/(4.5e9))**-n_out
        #Q=1.0
        #Q=1-(1-n_out)*(t0/4e9)**2.0 # alternative for lowering outgassing
        F_outg=F_outgass*Q**mm
        #F_outg=F_outgass*Q # alternative for lowering outgassing
        spread= Q 

        if landf_vary=="y":
            land=lf(t0,lfrac,growth_timing)
        else:
            land=1.0

        if lum_vary=="y":
            L_over_Lo = 1/(1+0.4*(t0/4.6e9)) # Evolution of solar luminosity from Gough 1981.
            if options_array[2]==1: # Check to see if methane option is on. If yes, use methane climate models:
                ### Weighting functions for Phanerozoic, Proterozoic, and Archean methane climates (only used high methane tests)   
                ### Smooth weighting functions are required because sudden temperature jumps break the model.         
                w0 = 1-1/(1+numpy.exp(-15*(t0/1e9-0.541))) # 1 in the Phanerozoic, 0 elsewhere
                w1 = 1/(1+numpy.exp(-15*(t0/1e9-0.541)))-1/(1+numpy.exp(-15*(t0/1e9-2.5))) # 1 in the Proterozoic, 0 elsewhere
                w2 = 1/(1+numpy.exp(-15*(t0/1e9-2.5))) # 1 in the Archean, 0 elsewhere
                arch_temp=clim_fun_highCH4(EE_ppCO2_o,L_over_Lo) ## Archean temperature with high methane
                Protero_temp = clim_fun_lowCH4(EE_ppCO2_o,L_over_Lo) ## Proterozoic temperature with high methane
                
                ### Climate model for eleveated Proterozoic and Archean methane abundances:
                T_surface=(w0*clim_fun_CO2_only(EE_ppCO2_o,L_over_Lo)+w1*Protero_temp+w2*arch_temp)/(w0+w1+w2)

            else: # methane option off -> use CO2 only climate model (nominal model)
                T_surface=clim_fun_CO2_only(EE_ppCO2_o,L_over_Lo)
                ###T_surface=clim_fun_CO2_only(EE_ppCO2_o,L_over_Lo)+30.0 * (1/(1+numpy.exp(-15*(t0/1e9-2.5)))) ## for +30 K sensitivity test
        else:
            T_surface=clim_fun_CO2_only(EE_ppCO2_o,1.0) # For CO2-only, no evolution in solar luminosity

        ## Calculate deep ocean temperature
        T_surface_diff=T_surface-T_surface0
        interc=274.037-deep_grad*T_surface0 # intercept chosen to reproduce initial (modern) temperature
        buffer_T=deep_grad*T_surface+interc 
        buffer_Ti=numpy.max([deep_grad*T_surface+interc,271.15]) # ensures deep ocean is not frozen
        buffer_T=numpy.min([buffer_Ti,T_surface]) # This is deep ocean temperature (cannot be warmer than the surface)
    
        ## Calculate pore space temperature, 
        ## accounting for effect of sediment thickness and heatflow
        Pore_mod_dynamic=Pore_mod
        sed_depth=700-700*(1-sed_thick)*t0/4e9 # equation S5
        Pore_mod_dynamic=Q*sed_depth/77.7777777 ## Difference is 9 K for Q=1
        T_pore_space=buffer_T+Pore_mod_dynamic
        #T_pore_space=buffer_T+9.0 # sensitivity test constant temperature difference between pore space and deep ocean
        
        T_upd=T_surface #update surface temperature for next iteration
        T_upd2=T_pore_space #update pore space temperature for next iteration
 
 
    # Calculate saturation state of ocean and pore space (equation S18): 
    # combined activity coefficients for Ca and CO3. Assumed to be 1 except for Ca sensitiviy test where complexing is important
    act_o=1.0 #-7.89*EE_Ca_o**3+10.847*EE_Ca_o**2-5.2246*EE_Ca_o+1.0551 #from Aspen polynomial fit Fig. S20, or 1.0 for nominal case
    act_p=1.0 #-7.89*EE_Ca_p**3+10.847*EE_Ca_p**2-5.2246*EE_Ca_p+1.0551 #from Aspen polynomial fit Fig. S20, or 1.0 for nominal case
    EE_omega_o= act_o*EE_Ca_o *EE_co3_o/(Sol_prod(T_surface)) # saturation state ocean
    EE_omega_p=act_p*EE_Ca_p *EE_co3_p/(Sol_prod(T_pore_space)) # saturation state pore space                
    
    #Calculate carbonate sinks
    Carb_sink=ko1*land*(EE_omega_o - 1)**n+ko2*(EE_omega_o**2.84) # ocean carbonate from equation S19, noting that ko2 is 0 in this version of the model      
    Precip_ocean=Carb_sink ## Precipitaion of carbonate in ocean
    Precip_pore=k_c*(EE_omega_p -1)**n ## Precipitation of carbonate in pore space

    # Carbonate weathering flux 
    carb_weath=F_carbw*biology_modifier*(EE_ppCO2_o/ppCO2_o)**(carb_exp)*numpy.exp((T_surface_diff)/tdep_weath)*land # equation S2
            
    ### Dissolution of seafloor flux:
    F_diss=k_r*2.88*10**-14*10**(-coef_for_diss*EE_pH_p)*spread**beta*numpy.exp(-Ebas/(8.314*(T_pore_space))) # equation S3
    #F_diss=k_r*2.88*10**-14*10**(-coef_for_diss*8.57)*numpy.exp(-Ebas/(8.314*(283))) # sensitivity test holding dissolution artifically constant
    
    #Continetnal silicate weathering flux
    F_silic=k_w*biology_modifier*land*(EE_ppCO2_o/ppCO2_o)**climp*numpy.exp((T_surface_diff)/tdep_weath) # equation 1
    
    #return selected carbon cycle variables and fluxes   
    return [EE_pH_o,EE_co3_o,EE_hco3_o, EE_co2aq_o,EE_ppCO2_o,EE_Ca_o,EE_omega_o,T_surface,buffer_T,F_diss,F_silic,Precip_ocean,Precip_pore,EE_omega_p,F_outg,carb_weath,EE_co3_p,EE_hco3_p, EE_co2aq_p,EE_pH_p,EE_ppCO2_p,EE_Ca_p,EE_H_o,EE_H_p,T_surface_diff,F_outg,spread,T_pore_space] #T_pore_space
      

