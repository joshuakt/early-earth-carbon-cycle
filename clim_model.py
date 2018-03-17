## This script contains polynomial parameterizations of the three climate models used in this study.
## Each function takes as inputs atmospheric pCO2 (bar) and solar luminosity relative to modern, and returns surface temperature.
## The functions describing these polynomials fits are described in the supplementary material to the manuscript (Appendix E).

import numpy as np
import pylab

# This function is the CO2-only climate model
def clim_fun_CO2_only(CO2,S): 
    # polynomial coefficients from fit to 1-D climate model:
    coeff=np.array([  3145.89156036,   -894.19626463, -14231.92214014,  -3791.69941395,  18527.77037214, -33560.74243889,  26297.60786657,  -7674.7687673 ,  4461.16266511,  -1569.24844923,  11329.24842061, -21270.62144395,  -14022.32137627,  26860.83381806,   7684.55056014,   5722.63936369,  -178.26420159,   -396.14681854,   6399.09085414,    875.42303333,  -1605.36365431,   1304.94455772,  -1569.30078917,  -8012.66961608,  -3415.93123913])
    x=pylab.log10(CO2) # polynomial fit done in log space
    y=S
    inputs=[x*0+1, x, y, x**2, x**2*y, x**2*y**2, y**2, x*y**2, x*y,x**3,x**3*y**3,y**3,x**3*y**2,x**2*y**3,x**3*y,x*y**3,  x**4,x**4*y**4,y**4,x**4*y,x**4*y**2,x**4*y**3,y**4*x,y**4*x**2,y**4*x**3]
    return np.sum(coeff*inputs)

# 100 ppm methane representing the Proterozoic
def clim_fun_lowCH4(CO2,S): #CO2 in bar, S relative to S0
    # polynomial coefficients from fit to 1-D climate model:
    coeff=np.array([-2.803182960364228649e+01,-3.557582084216555813e+02,1.134404914511084826e+03,-1.969683575705347778e+02,8.939465316235038017e+02,-1.178979589385960026e+03,-1.253544234276430871e+03,-2.246677400515449335e+03,1.716529689916493226e+03,-2.795988531878487038e+01,6.842046486014059781e+01,5.012078688730417184e+02,-1.638805507113463875e+02,4.888049629422617954e+02,1.241243707536893623e+02,9.177501872534526228e+02])
    x=pylab.log10(CO2) # polynomial fit done in log space
    y=S
    inputs=[x*0+1, x, y, x**2, x**2*y, x**2*y**2, y**2, x*y**2, x*y,x**3,x**3*y**3,y**3,x**3*y**2,x**2*y**3,x**3*y,x*y**3]
    return np.sum(coeff*inputs)

# 1% methane representing the Archean
def clim_fun_highCH4(CO2,S): 
    # polynomial coefficients from fit to 1-D climate model:
    coeff=np.array([ -7.54993474e+00,  -1.56537857e+02,   1.09569449e+03, -4.23166796e+01,   3.33328419e+02,  -5.36250637e+02,-1.24044934e+03,  -1.54156844e+03,   1.04071399e+03,-1.46500650e+00,   2.45480651e+01,   5.13587132e+02, -4.87777690e+01,   2.51789370e+02,   2.63929259e+01, 6.86937462e+02])
    x=pylab.log10(CO2) # polynomial fit done in log space
    y=S
    inputs=[x*0+1, x, y, x**2, x**2*y, x**2*y**2, y**2, x*y**2, x*y,x**3,x**3*y**3,y**3,x**3*y**2,x**2*y**3,x**3*y,x*y**3]
    return np.sum(coeff*inputs)
