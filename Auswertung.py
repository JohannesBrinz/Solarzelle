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

U_err = 0
I_err = 0

T = 313

def U_lin(I, R, x):
    return I*R + x

params_lin, params_lin_cov = optimize.curve_fit(U_lin, a_si_du["current / A"][800:970],\
 a_si_du["voltage / V"][800:970])       #Linerarer Fit zur Ermittlung des Serienwiderstands R_s

print("Serienwiderstand R_s:    ", params_lin[0] )

R = params_lin[0]
R = R-0.13     #-0.13 damit der ln nicht negariv wird, variation

def U(I, n, I_s, x):
    return ((n*constants.k*T)/constants.e) * np.log(I/I_s) + x   #Funktion für U' = U - RI

params, params_cov = optimize.curve_fit(U, a_si_du["current / A"][700:900], a_si_du["voltage / V"][700:900]-R*a_si_du["current / A"][700:900])

#Plot Fit Bestimmung Si dunkel
plt.errorbar(a_si_du["voltage / V"][0:970], a_si_du["current / A"][0:970], xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "blue", capsize=3)  #Normaler Plot

plt.errorbar(U_lin(np.linspace(0, 1, 1000), params_lin[0], params_lin[1]), np.linspace(0, 1, 1000), \
 lw=2, fmt = "--", c = "black")   #Linerarer Fit

plt.errorbar(a_si_du["voltage / V"][700:970]-R*a_si_du["current / A"][700:970], a_si_du["current / A"][700:970], xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "green", capsize=3)   #Plot U'

plt.errorbar(U(np.linspace(0, 0.6, 1000), params[0], params[1], params[2]), np.linspace(0, 0.6, 1000), \
  lw=2, fmt = "-.", c = "black")  #Fit Log

plt.title('U-I Kennlinie Si dunkel', fontsize = 15)
plt.text( -0.25, 0.2, "$R_S = $" + str(round(params_lin[0], 3)) + "$\Omega$", fontsize = 14)   #R_s
plt.text( -0.25, 0.3, "$n = $" + str(round(params[0], 3)), fontsize = 14)   #n
plt.text( -0.25, 0.4, "$I_S = $" + str(round(params[1], 3)) + "A", fontsize = 14)   #I_s
plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
plt.grid(True)
plt.legend(['U-I Kennlinie', "Linearer Fit", "Plot  $U' = U - R_s I$", "logarithmischer Fit"], fontsize = 13)
plt.savefig('Plots/Dunkel_Si.png', dpi=300)
plt.clf()


#Plot Kennlinien dunkel
Ao1 = 0.06
Ao2 =  25
Aa = 26

fig, ax1 = plt.subplots()
color = 'black'
ax1.set_xlabel('U [V]' , fontsize = 13)
ax1.set_ylabel('I [$mA/cm^2$]', fontsize = 13, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()
color = 'black'

ax2.set_ylabel('I [$A/cm^2$]', fontsize = 13, color=color)  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor=color)

ax1.errorbar(a_org_du["voltage / V"], a_org_du["current / A"]*1e3/Ao1, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "blue", capsize=3, label = "Kennlinie organisch 1")

ax1.errorbar(a_org2_du["voltage / V"], a_org2_du["current / A"]*1e5/Ao2, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "green", capsize=3, label = "Kennlinie organisch 2 multipliziert mit 100")

ax2.errorbar(a_si_du["voltage / V"], a_si_du["current / A"]/Aa, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "red", capsize=3)

plt.title('U-I Kennlinien dunkel', fontsize = 15)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.subplots_adjust(top=0.85)

plt.grid(True)
ax1.legend(loc = 2)
ax2.legend(["Kennlinie anorganisch"], loc = 4)
plt.savefig('Plots/Kennlinie_Dunkel.png', dpi=300)
plt.clf()


#Kennlinie hell
fig, ax1 = plt.subplots()
color = 'black'
ax1.set_xlabel('U [V]' , fontsize = 13)
ax1.set_ylabel('I [$mA/m^2$]', fontsize = 13, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()
color = 'black'

ax2.set_ylabel('I [$mA/m^2$]', fontsize = 13, color=color)  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor=color)

ax1.errorbar(a_org_11mV["voltage / V"], a_org_11mV["current / A"]*1e2/Ao1, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "blue", capsize=3)
ax2.errorbar(a_org2_11["voltage / V"], a_org2_11["current / A"]*1e2/Ao2, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "green", capsize=3)

ax1.legend(['Kennlinie organisch 1'], loc = 2)
ax2.legend(["Kennlinie organisch 2"], loc = 4)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.subplots_adjust(top=0.85)

plt.title('U-I Kennlinien hell organisch', fontsize = 15)

plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [$mA/m^2$]', fontsize = 13)
plt.grid(True)
#plt.legend(['U-I Kennlinie', "Linearer Fit", "Plot  $U' = U - R_s I$", "logarithmischer Fit"], fontsize = 13)
plt.savefig('Plots/Kennlinie_Hell.png', dpi=300)
plt.clf()


#Si Kennlinie hell

I2 = a_an_11mV["current / A"][780]
I1 = a_an_11mV["current / A"][770]
U1 = a_an_11mV["voltage / V"][770]
U2 = a_an_11mV["voltage / V"][780]

U_OC = U1 + (U2-U1) * 0.33             #lineare Interpolation U_OC
j_sc = a_an_11mV["current / A"][570]/Aa
P = a_an_11mV["voltage / V"]*(-a_an_11mV["current / A"])

#Finde I_MPP
max = 0
j = 0
for i in range(570, 770):
    if P[i] > max:
        max = P[i]
        j = i
I_MPP = a_an_11mV["current / A"][j]
U_MPP = a_an_11mV["voltage / V"][j]
FF = -P[j]/U_OC*j_sc*Aa
p_ein = 100*11/33.2
Wirkungsgrad = (P[j]*1e3/Aa)/p_ein            #sollte 0.22 sein S.10

#Plot Si hell

plt.errorbar(a_an_11mV["voltage / V"][530:880], a_an_11mV["current / A"][530:880], xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "blue", capsize=3)  #Normaler Plot
plt.errorbar(U_MPP, I_MPP, fmt = "s", color = "black")

plt.title('U-I Kennlinie Si hell', fontsize = 15)
plt.text( 0.2, 0.35, "$U_{Oc} = $" + str(round(U_OC, 3)) + " V", fontsize = 14)   #U_OC
plt.text( 0.2, 0.25, "$j_{Sc} = $" + str(round(j_sc*1e3, 3)) + r"$\frac{mA}{cm^2}$", fontsize = 14)   #j_sc
plt.text( 0.2, 0.15, "$I_{MPP} = $" + str(round(I_MPP, 3)) + " A", fontsize = 14)   #I_s
plt.text( 0.2, 0.05, "$U_{MPP} = $" + str(round(U_MPP, 3)) + " V", fontsize = 14)
plt.text( 0.2, -0.05, "$FF = $" + str(round(FF, 3)), fontsize = 14)
plt.text( 0.2, -0.15, "$\eta = $" + str(round(Wirkungsgrad, 3)), fontsize = 14)
plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
#plt.grid(True)
plt.legend(['U-I Kennlinie', "MPP"], fontsize = 13)
plt.savefig('Plots/Hell_Si.png', dpi=300)
plt.clf()


#_______B____________


#Si Kennlinie hell
I2_5 = b_5["current / A"][760]
I1_5 = b_5["current / A"][759]
U1_5 = b_5["voltage / V"][759]
U2_5 = b_5["voltage / V"][760]
I2_20 = b_20["current / A"][770]
I1_20 = b_20["current / A"][769]
U1_20 = b_20["voltage / V"][769]
U2_20 = b_20["voltage / V"][770]
I2_32 = b_sonne["current / A"][760]
I1_32 = b_sonne["current / A"][759]
U1_32 = b_sonne["voltage / V"][759]
U2_32 = b_sonne["voltage / V"][760]


U_OC_5 = U1_5 + (U2_5-U1_5) * 0.33
j_sc_5 = b_5["current / A"][570]/Aa
U_OC_20 = U1_20 + (U2_20-U1_20) * 0.33
j_sc_20 = b_20["current / A"][570]/Aa
U_OC_32 = U1_32 + (U2_32-U1_32) * 0.33
j_sc_32 = b_sonne["current / A"][570]/Aa


#Plot Si Intensitäten
plt.errorbar(b_5["voltage / V"][530:900], b_5["current / A"][530:900], xerr = U_err, yerr = I_err,\
            linewidth = 2, capsize=3) #5 mV
plt.errorbar(a_an_11mV["voltage / V"][530:900], a_an_11mV["current / A"][530:900], xerr = U_err, yerr = I_err,\
            linewidth = 2, capsize=3)  #0.3 Sonnen
plt.errorbar(b_20["voltage / V"][530:900], b_20["current / A"][530:900], xerr = U_err, yerr = I_err,\
            linewidth = 2, capsize=3) #0.6 Sonnen
plt.errorbar(b_sonne["voltage / V"][530:900], b_sonne["current / A"][530:900], xerr = U_err, yerr = I_err,\
            linewidth = 2, capsize=3) #Sonne

plt.title('U-I Kennlinie Si hell', fontsize = 15)
'''plt.text( 0.5, -0.2, "$U_{Oc,5mV} = $" + str(round(U_OC_5, 3)) + " V", fontsize = 12)   #U_OC
plt.text( 0.5, -0.3, "$j_{Sc,5mV} = $" + str(round(j_sc_5*1e3, 3)) + r"$\frac{mA}{cm^2}$", fontsize = 12)   #j_sc
plt.text( 0.5, -0.4, "$U_{Oc, 11mV} = $" + str(round(U_OC, 3)) + " V", fontsize = 12)   #U_OC
plt.text( 0.5, -0.5, "$j_{Sc, 11mV} = $" + str(round(j_sc*1e3, 3)) + r"$\frac{mA}{cm^2}$", fontsize = 12)   #j_sc
plt.text( 0.5, -0.6, "$U_{Oc, 20mV} = $" + str(round(U_OC_20, 3)) + " V", fontsize = 12)   #U_OC
plt.text( 0.5, -0.7, "$j_{Sc, 20mV} = $" + str(round(j_sc_20*1e3, 3)) + r"$\frac{mA}{cm^2}$", fontsize = 12)   #j_sc
plt.text( 0.5, -0.8, "$U_{Oc, 32mV} = $" + str(round(U_OC_32, 3)) + " V", fontsize = 12)   #U_OC
plt.text( 0.5, -0.9, "$j_{Sc, 32mV} = $" + str(round(j_sc_32*1e3, 3)) + r"$\frac{mA}{cm^2}$", fontsize = 12) '''  #j_sc

plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
plt.grid(True)
plt.legend(['$U_{Oc,5mV} = $' + str(round(U_OC_5, 2)) + ' V,  ' + "$j_{Sc,5mV} = $" + str(round(j_sc_5*1e3, 2)) + r"$\frac{mA}{cm^2}$",\
 '$U_{Oc,11mV} = $' + str(round(U_OC, 2)) + ' V,  ' + "$j_{Sc,11mV} = $" + str(round(j_sc*1e3, 2)) + r"$\frac{mA}{cm^2}$",\
 '$U_{Oc,20mV} = $' + str(round(U_OC_20, 2)) + ' V,  ' + "$j_{Sc,20mV} = $" + str(round(j_sc_20*1e3, 2)) + r"$\frac{mA}{cm^2}$",\
  '$U_{Oc,32mV} = $' + str(round(U_OC_32, 2)) + ' V,  ' + "$j_{Sc,32mV} = $" + str(round(j_sc_32*1e3, 2)) + r"$\frac{mA}{cm^2}$"], fontsize = 12)
plt.savefig('Plots/Intensitäten.png', dpi=300)
plt.clf()

# Kurzschlussstrom über Intensität
#linearer Fit
def lin(x, m, a, b):
    return m*(x-a) + b

j_SC = pd.DataFrame()
j_SC["j_SC"] = [j_sc_5, j_sc, j_sc_20, j_sc_32]
j_SC["i"] = [500/32, 1100/32, 2000/32, 100]

params, params_cov = optimize.curve_fit(lin, j_SC["i"], j_SC["j_SC"])

#Plot
plt.title("Intensity Dependence of $j_{sc}$", fontsize = 18)

plt.errorbar(np.linspace(15, 100, 1000), lin(np.linspace(15, 100, 1000), params[0], params[1], params[2]), color = "black")
plt.errorbar(j_SC["i"], j_SC["j_SC"], fmt = "d", linewidth = 2, capsize=3, markersize='14', color = "blue")

'''plt.errorbar(500/32, j_sc_5, fmt = "s", markersize='16')
plt.errorbar(1100/32, j_sc, fmt = "d", markersize='16')
plt.errorbar(2000/32, j_sc_20, fmt = "x", markersize='16')
plt.errorbar(100, j_sc_32, fmt = ".", markersize='16')'''

plt.xlabel(r'$I_0 [\frac{mW}{cm^2}]$' , fontsize = 14)
plt.ylabel(r'j [$\frac{A}{cm^2}]$', fontsize = 14)
plt.legend(["linear fit", "data"])
plt.grid(True)
plt.savefig('Plots/j_Intensität.png', dpi=300)
plt.clf()

# U_Oc über Intensität
plt.title("Intensity Dependence of $U_{OC}$", fontsize = 18)

plt.errorbar(500/32, U_OC_5, fmt = "s", markersize='16')
plt.errorbar(1100/32, U_OC, fmt = "d", markersize='16')
plt.errorbar(2000/32, U_OC_20, fmt = "x", markersize='16')
plt.errorbar(3200/32, U_OC_32, fmt = ".", markersize='16')

plt.xlabel(r'$I_0 [\frac{mW}{cm^2}]$' , fontsize = 14)
plt.ylabel(r'$U_{OC} [V]$', fontsize = 14)

plt.savefig('Plots/U_Intensität.png', dpi=300)
plt.clf()

#_______________C_______________
