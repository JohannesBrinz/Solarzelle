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
Kippunng = pd.read_csv("Daten/Kippung.dat", sep = "\t\t", header = 0, \
    names = ["Winkel[rad]", "I_K,o[yA]", "I_K,ano[mA]"], engine='python')    #40 < T [C] < 46
Leerlaufspannung = pd.read_csv("Daten/Leerlaufspannungen.dat", sep = "\t", header = 0, \
    names = ["T[C]", "U_L[mV]"])
a_an_11mV = pd.read_csv("Daten/a_an_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_an_du = pd.read_csv("Daten/a_an_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org2_11 = pd.read_csv("Daten/a_org2_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org2_du = pd.read_csv("Daten/a_org2_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org_11mV = pd.read_csv("Daten/a_org_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_org_du = pd.read_csv("Daten/a_org_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_si_11mV = pd.read_csv("Daten/a_si_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
a_si_du = pd.read_csv("Daten/a_si_du.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
b_5 = pd.read_csv("Daten/b_5.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
b_20 = pd.read_csv("Daten/b_20.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
b_sonne = pd.read_csv("Daten/b_sonne.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c1_11mV = pd.read_csv("Daten/c1_11mV.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c3_v1_11 = pd.read_csv("Daten/c3_v1_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c3_v2_11 = pd.read_csv("Daten/c3_v2_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c3_v3_11 = pd.read_csv("Daten/c3_v3_11.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
c4_11verb = pd.read_csv("Daten/c4_11verb.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
d_30 = pd.read_csv("Daten/d_30.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
r68_s2 = pd.read_csv("Daten/r68_s2.4.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
r68_s3 = pd.read_csv("Daten/r68_s3.6.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])
r73_s2 = pd.read_csv("Daten/r73_s2,6.dat", sep = "\t", header = 8, \
    names = ["voltage / V", "current / A", "time"])



#_______A________
#Fit
T = 313

def U_lin(I, R, x):
    return I*R + x

params_lin, params_lin_cov = optimize.curve_fit(U_lin, a_si_du["current / A"][800:970],\
 a_si_du["voltage / V"][800:970])       #Linerarer Fit zur Ermittlung des Serienwiderstands R_s

print("Serienwiderstand R_s:    ", params_lin[0] )

R = params_lin[0]
R = R-0.01          #-0.01 damit der ln nicht negariv wird

def U(I, I_s, n):
    return ((n*constants.k*T)/constants.e) * np.log((I+I_s)/I_s)    #Funktion für U' = U - RI

params, params_cov = optimize.curve_fit(U, a_si_du["voltage / V"][0:970]-R*a_si_du["current / A"][0:970],\
 a_si_du["current / A"][0:970]) #Der Fit funktioniert nicht


#Plot
plt.errorbar(a_si_du["voltage / V"][0:970], a_si_du["current / A"][0:970], xerr = 0, yerr = 0,\
            linewidth = 2, color = "green", capsize=3)
plt.errorbar(U_lin(np.linspace(0, 1, 1000), params_lin[0], params_lin[1]), np.linspace(0, 1, 1000), \
 lw=1, fmt = "--", c = "black")
plt.errorbar(U(np.linspace(0.01, 1, 1000), params[0], params[1]), np.linspace(0.01, 1, 1000), \
  lw=1, fmt = "--", c = "black")

plt.errorbar(a_si_du["voltage / V"][0:970]-R*a_si_du["current / A"][0:970], a_si_du["current / A"][0:970], xerr = 0, yerr = 0,\
            linewidth = 2, color = "green", capsize=3)      #Plot U'

plt.title('U-I Kennlinie Si dunkel', fontsize = 15)
plt.text( 0, 0.8, "$R_S = $" + str(round(params_lin[0], 3)) + "$\Omega$")   #R_s
plt.text( 0, 0.7, "$n = $" + str(round(params[0], 3)))   #n
plt.text( 0, 0.6, "$I_S = $" + str(round(params[1], 3)) + "A")   #I_s
plt.xlabel('U [mV]' , fontsize = 13)
plt.ylabel('I [mA]', fontsize = 13)
plt.grid(True)
#plt.legend(['experimental data', "weighted modes", "identical modes"], fontsize = 13)
plt.savefig('Plots/Dunkel_Si.png', dpi=300)
plt.clf()
