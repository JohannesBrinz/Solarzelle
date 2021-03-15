''' Versuch Solarzelle
 Johannes Brinz & Caterina Vanelli
 Datum: 12.03.2021
 '''

import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from scipy import optimize
from scipy import constants
import math
import matplotlib.font_manager as fm
import matplotlib.mlab as mlab
from scipy.stats import norm


#Datenimport
Kippunng = pd.read_csv("Vanelli Brinz 12-03-21/Kippung.dat", sep = "\t\t", header = 0, \
    names = ["Winkel[rad]", "I_K,o[yA]", "I_K,ano[mA]"], engine='python')    #40 < T [C] < 46
Leerlaufspannung = pd.read_csv("Vanelli Brinz 12-03-21/Leerlaufspannungen.dat", sep = "\t", header = 0, \
    names = ["T[C]", "U_L[mV]"])
a_an_11mV = pd.read_csv("Vanelli Brinz 12-03-21/a_an_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_an_du = pd.read_csv("Vanelli Brinz 12-03-21/a_an_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org2_11 = pd.read_csv("Vanelli Brinz 12-03-21/a_org2_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org2_du = pd.read_csv("Vanelli Brinz 12-03-21/a_org2_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org_11mV = pd.read_csv("Vanelli Brinz 12-03-21/a_org_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org_du = pd.read_csv("Vanelli Brinz 12-03-21/a_org_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_si_11mV = pd.read_csv("Vanelli Brinz 12-03-21/a_si_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_si_du = pd.read_csv("Vanelli Brinz 12-03-21/a_si_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
b_5 = pd.read_csv("Vanelli Brinz 12-03-21/b_5.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
b_20 = pd.read_csv("Vanelli Brinz 12-03-21/b_20.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
b_sonne = pd.read_csv("Vanelli Brinz 12-03-21/b_sonne.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c1_11mV = pd.read_csv("Vanelli Brinz 12-03-21/c1_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c3_v1_11 = pd.read_csv("Vanelli Brinz 12-03-21/c3_v1_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c3_v2_11 = pd.read_csv("Vanelli Brinz 12-03-21/c3_v2_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c3_v3_11 = pd.read_csv("Vanelli Brinz 12-03-21/c3_v3_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c4_11verb = pd.read_csv("Vanelli Brinz 12-03-21/c4_11verb.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
d_30 = pd.read_csv("Vanelli Brinz 12-03-21/d_30.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
r68_s2 = pd.read_csv("Vanelli Brinz 12-03-21/r68_s2.4.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
r68_s3 = pd.read_csv("Vanelli Brinz 12-03-21/r68_s3.6.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
r73_s2 = pd.read_csv("Vanelli Brinz 12-03-21/r73_s2,6.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])



#______A________
#Fit
T = 313
def U(I, I_s, n, R):
    return ((n*constants.k*T)/constants.e) * np.log((I+I_s)/I_s) + I*R



params, params_cov = optimize.curve_fit(U, a_si_du["current / A"][0:970],\
 a_si_du["voltage / V"][0:970])

#Plot
plt.errorbar(a_si_du["voltage / V"][0:970], a_si_du["current / A"][0:970], xerr = 0, yerr = 0,\
            linewidth = 2, color = "green", capsize=3)
plt.errorbar(U(np.linspace(-1, 1, 1000), params[0], params[1], params[2]), np.linspace(-1, 1, 1000), \
 lw=1, fmt = "--", c = "black")
#plt.errorbar(np.linspace(0, 0.5, 1000), np.linspace((1/2.718),(1/2.718), 1000))
plt.title('U-I Kennlinie Si dunkel', fontsize = 15)
#plt.text( 0, 0,"$L_c = $" + str(round(t_c * 3*10e8, 3)))
plt.xlabel('U [mV]' , fontsize = 13)
plt.ylabel('I [mA]', fontsize = 13)
plt.grid(True)
#plt.legend(['experimental data', "weighted modes", "identical modes"], fontsize = 13)
plt.savefig('Plots/Dunkel_Si.png', dpi=300)
plt.clf()
