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
Kippung = pd.read_csv("Daten/Kippung.dat", sep = "\t\t", header = 0, \
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
I_err = 0    #Temperaturabhngig!!!!

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

plt.title('U-I characteristics Si dark', fontsize = 15)
plt.text( -0.25, 0.2, "$R_S = $" + str(round(params_lin[0], 3)) + "$\Omega$", fontsize = 14)   #R_s
plt.text( -0.25, 0.3, "$n = $" + str(round(params[0], 3)), fontsize = 14)   #n
plt.text( -0.25, 0.4, "$I_S = $" + str(round(params[1], 3)) + "A", fontsize = 14)   #I_s
plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
plt.grid(True)
plt.legend(['U-I characteristics', "Linear fit", "Plot  $U' = U - R_s I$", "Logarithmic fit"], fontsize = 13)
plt.savefig('Plots/Dunkel_Si.png', dpi=300)
plt.clf()


#Plot Characteristicsn dunkel
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
            linewidth = 2, color = "blue", capsize=3, label = "Characteristics organic 1")

ax1.errorbar(a_org2_du["voltage / V"], a_org2_du["current / A"]*1e5/Ao2, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "green", capsize=3, label = "Characteristics organic 2 x 100")

ax2.errorbar(a_si_du["voltage / V"], a_si_du["current / A"]/Aa, xerr = U_err, yerr = I_err,\
            linewidth = 2, color = "red", capsize=3)

plt.title('U-I characteristics dark', fontsize = 15)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.subplots_adjust(top=0.85)

plt.grid(True)
ax1.legend(loc = 2)
ax2.legend(["Characteristics anorganic"], loc = 4)
plt.savefig('Plots/Kennlinie_Dunkel.png', dpi=300)
plt.clf()


#Characteristics hell
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

ax1.legend(['Characteristics organic 1'], loc = 2)
ax2.legend(["Characteristics organic 2"], loc = 4)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.subplots_adjust(top=0.85)

plt.title('U-I Characteristics organic illuminated', fontsize = 15)

plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [$mA/m^2$]', fontsize = 13)
plt.grid(True)
#plt.legend(['U-I Characteristics', "Linearer Fit", "Plot  $U' = U - R_s I$", "logarithmischer Fit"], fontsize = 13)
plt.savefig('Plots/Kennlinie_Hell.png', dpi=300)
plt.clf()


#Si Characteristics hell

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

plt.errorbar(np.linspace(0, 0.8, 1000), np.linspace(0, 0, 1000), color = "black", linewidth = "1")

plt.title('U-I characteristics Si illuminated', fontsize = 15)
plt.text( 0.2, 0.35, "$U_{Oc} = $" + str(round(U_OC, 3)) + " V", fontsize = 14)   #U_OC
plt.text( 0.2, 0.25, "$j_{Sc} = $" + str(round(j_sc*1e3, 3)) + r"$\frac{mA}{cm^2}$", fontsize = 14)   #j_sc
plt.text( 0.2, 0.15, "$I_{MPP} = $" + str(round(I_MPP, 3)) + " A", fontsize = 14)   #I_s
plt.text( 0.2, 0.05, "$U_{MPP} = $" + str(round(U_MPP, 3)) + " V", fontsize = 14)
plt.text( 0.2, -0.05, "$FF = $" + str(round(FF, 3)), fontsize = 14)
plt.text( 0.2, -0.15, "$\eta = $" + str(round(Wirkungsgrad, 3)), fontsize = 14)
plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
plt.xlim(xmin=0, xmax = 0.8)
plt.grid(b = True, which = "major", axis = "y")
plt.legend(['U-I characteristics', "MPP"], fontsize = 13)
plt.savefig('Plots/Hell_Si.png', dpi=300)
plt.clf()


#_______B____________


#Si Characteristics hell
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

plt.title('U-I characteristics Si illuminated', fontsize = 15)
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
plt.title("Intensity dependence of $j_{sc}$", fontsize = 18)

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
plt.title("Intensity dependence of $U_{OC}$", fontsize = 18)

plt.errorbar(500/32, U_OC_5, fmt = "x", markersize='13', c = "blue")
plt.errorbar(1100/32, U_OC, fmt = "x", markersize='13', c = "blue")
plt.errorbar(2000/32, U_OC_20, fmt = "x", markersize='13', c = "blue")
plt.errorbar(3200/32, U_OC_32, fmt = "x", markersize='13', c = "blue")

plt.xlabel(r'$I_0 [\frac{mW}{cm^2}]$' , fontsize = 14)
plt.ylabel(r'$U_{OC} [V]$', fontsize = 14)
plt.grid(True)
plt.savefig('Plots/U_Intensität.png', dpi=300)
plt.clf()

#_______________C_______________
#Si Characteristics hell
#1
I2_1 = c1_11mV["current / A"][350]
I1_1 = c1_11mV["current / A"][349]
U1_1 = c1_11mV["voltage / V"][349]
U2_1 = c1_11mV["voltage / V"][350]



U_OC_1 = U1_1 + (U2_1-U1_1) * 0.5           #lineare Interpolation U_OC
j_sc_1 = c1_11mV["current / A"][150]/Aa
P_1 = c1_11mV["voltage / V"]*(-c1_11mV["current / A"])

#Finde I_MPP
max_1 = 0
j_1 = 0
for i in range(150, 349):
    if P_1[i] > max_1:
        max_1 = P_1[i]
        j_1 = i
I_MPP_1 = c1_11mV["current / A"][j_1]
U_MPP_1 = c1_11mV["voltage / V"][j_1]
FF_1 = -P_1[j_1]/U_OC_1*j_sc_1*Aa*6
p_ein = 100*11/33.2
Wirkungsgrad_1 = (P_1[j_1]*1e3/(Aa*6))/p_ein

#2
I2_2 = r68_s2["current / A"][350]
I1_2 = r68_s2["current / A"][349]
U1_2 = r68_s2["voltage / V"][349]
U2_2 = r68_s2["voltage / V"][350]

U_OC_2 = U1_2 + (U2_2-U1_2) * 0.5           #lineare Interpolation U_OC
j_sc_2 = r68_s2["current / A"][50]/Aa*6
P_2 = r68_s2["voltage / V"]*(-r68_s2["current / A"])

#Finde I_MPP
max_2 = 0
j_2 = 0
for i in range(150, 349):
    if P_2[i] > max_2:
        max_2 = P_2[i]
        j_2 = i
I_MPP_2 = r68_s2["current / A"][j_2]
U_MPP_2 = r68_s2["voltage / V"][j_2]
FF_2 = -P_2[j_2]/U_OC_2*j_sc_2*Aa*6
p_ein = 100*11/33.2
Wirkungsgrad_2 = (P_2[j_2]*1e3/(Aa*6))/p_ein

#3
I2_3 = r68_s3["current / A"][345]
I1_3 = r68_s3["current / A"][344]
U1_3 = r68_s3["voltage / V"][344]
U2_3 = r68_s3["voltage / V"][345]

U_OC_3 = U1_3 + (U2_3-U1_3) * 0.5           #lineare Interpolation U_OC
j_sc_3 = r68_s3["current / A"][50]/Aa*6
P_3 = r68_s3["voltage / V"]*(-r68_s3["current / A"])

#Finde I_MPP
max_3 = 0
j_3 = 0
for i in range(150, 344):
    if P_3[i] > max_3:
        max_3 = P_3[i]
        j_3 = i
I_MPP_3 = r68_s3["current / A"][j_3]
U_MPP_3 = r68_s3["voltage / V"][j_3]
FF_3 = -P_3[j_3]/U_OC_3*j_sc_3*Aa*6
p_ein = 100*11/33.2
Wirkungsgrad_3 = (P_3[j_3]*1e3/(Aa*6))/p_ein

#4
I2_4 = r73_s2["current / A"][340]
I1_4 = r73_s2["current / A"][339]
U1_4 = r73_s2["voltage / V"][339]
U2_4 = r73_s2["voltage / V"][340]

U_OC_4 = U1_4 + (U2_4-U1_4) * 0.5           #lineare Interpolation U_OC
j_sc_4 = r73_s2["current / A"][150]/Aa*6
P_4 = r73_s2["voltage / V"]*(-r73_s2["current / A"])

#Finde I_MPP
max_4 = 0
j_4 = 0
for i in range(150, 339):
    if P_4[i] > max_4:
        max_4 = P_4[i]
        j_4 = i
I_MPP_4 = r73_s2["current / A"][j_4]
U_MPP_4 = r73_s2["voltage / V"][j_4]
FF_4 = -P_4[j_4]/U_OC_4*j_sc_4*Aa*6
p_ein = 100*11/33.2
Wirkungsgrad_4 = (P_4[j_4]*1e3/(Aa*6))/p_ein

#Verbraucher
I2_Verb = c4_11verb["current / A"][219]
I1_Verb = c4_11verb["current / A"][218]
U1_Verb = c4_11verb["voltage / V"][218]
U2_Verb = c4_11verb["voltage / V"][219]

U_OC_Verb = U1_Verb + (U2_Verb-U1_Verb) * 0.5           #lineare Interpolation U_OC
j_sc_Verb = c4_11verb["current / A"][50]/Aa*12
P_Verb = c4_11verb["voltage / V"]*(-c4_11verb["current / A"])

#Finde I_MPP
max_Verb = 0
j_Verb = 0
for i in range(100, 218):
    if P_Verb[i] > max_Verb:
        max_Verb = P_Verb[i]
        j_Verb = i
I_MPP_Verb = c4_11verb["current / A"][j_Verb]
U_MPP_Verb = c4_11verb["voltage / V"][j_Verb]
FF_Verb = -P_Verb[j_Verb]/U_OC_Verb*j_sc_Verb*Aa*12
p_ein = 100*11/33.2
Wirkungsgrad_Verb = (P_Verb[j_Verb]*1e3/(Aa*12))/p_ein

print("\n\n\nFolgende Werte erhalten wir für die Kurzschlussströme:\n", "j_1: ", j_sc_1, "\n j_2: ", j_sc_2, "\n j_3: ", j_sc_3, \
"\n j_4: ", j_sc_4, "\nj_Verb: ", j_sc_Verb)
print("\n\nFolgende Werte erhalten wir für die Leerlaufspannungen:\n", "U_1: ", U_OC_1, "\n U_2: ", U_OC_2, "\n U_3: ", U_OC_3, \
"\n U_4: ", U_OC_4, "\nU_Verb: ", U_OC_Verb)
print("\n\nFolgende Werte erhalten wir für die Füllfaktoren:\n", "FF_1: ", FF_1, "\n FF_2: ", j_sc_2, "\n FF_3: ", j_sc_3, \
"\n FF_4: ", j_sc_4, "\nFF_Verb: ", j_sc_Verb)
print("\n\nFolgende Werte erhalten wir für die Wirkungsgrade:\n", "Wirkungsgrad_1: ", Wirkungsgrad_1, "\n Wirkungsgrad_2: ", Wirkungsgrad_2, "\n Wirkungsgrad_3: ", Wirkungsgrad_3, \
"\n Wirkungsgrad_4: ", Wirkungsgrad_4, "\nWirkungsgrad_Verb: ", Wirkungsgrad_Verb, "\n\n\n")

#Plot
plt.errorbar(c1_11mV["voltage / V"], c1_11mV["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")
plt.errorbar(r68_s2["voltage / V"], r68_s2["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")
plt.errorbar(r68_s3["voltage / V"], r68_s3["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")
plt.errorbar(r73_s2["voltage / V"], r73_s2["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")

plt.errorbar(c4_11verb["voltage / V"][50:219], c4_11verb["current / A"][50:219], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")

plt.title('U-I characteristics for different wirings', fontsize = 15)
plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
#plt.xlim(xmin=0, xmax = 0.8)
plt.grid(b = True, which = "major", axis = "y")
plt.legend(['Module', "$R_{ser} = 68 \Omega , R_{par} = 2.4 \Omega$", "$R_{ser} = 68 \Omega , R_{par} = 3.6 \Omega$",\
"$R_{ser} = 73 \Omega , R_{par} = 2.6 \Omega$", "Ventilator"], fontsize = 13)
plt.savefig('Plots/Verschaltung.png', dpi=300)
plt.clf()


#Verschattung
plt.errorbar(c1_11mV["voltage / V"], c1_11mV["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")
plt.errorbar(c3_v1_11["voltage / V"], c3_v1_11["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")
plt.errorbar(c3_v2_11["voltage / V"], c3_v2_11["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")
plt.errorbar(c3_v3_11["voltage / V"], c3_v3_11["current / A"], xerr = U_err, yerr = I_err,\
            linewidth = 1, fmt = ".")


plt.title('U-I characteristics for different shadings', fontsize = 15)
plt.xlabel('U [V]' , fontsize = 13)
plt.ylabel('I [A]', fontsize = 13)
#plt.xlim(xmin=0, xmax = 0.8)
plt.grid(b = True, which = "major", axis = "y")
plt.legend(['shading 0', "shading 1/6","shading 1/3", "shading 2/3"], fontsize = 13)
plt.savefig('Plots/Verschattung.png', dpi=300)
plt.clf()


#_______________D________________
T_err = 2
U_L_err = 0

params, params_cov = optimize.curve_fit(lin, Leerlaufspannung["T[C]"][0:5], Leerlaufspannung["U_L[mV]"][0:5])

plt.title("Temperature dependence of $U_{OC}$", fontsize = 18)

plt.errorbar(np.linspace(27, 55, 1000), lin(np.linspace(27, 55, 1000), params[0], params[1], params[2]), color = "black")
plt.errorbar(Leerlaufspannung["T[C]"], Leerlaufspannung["U_L[mV]"], xerr=T_err, yerr=U_L_err,\
fmt = "x", linewidth = 2, capsize=3, markersize='12', color = "blue")
plt.xlabel('T[K]' , fontsize = 14)
plt.ylabel('$U_{OC}$[mV]', fontsize = 14)
plt.legend(["linear fit", "data"])
plt.grid(True)
plt.savefig('Plots/Temperatur.png', dpi=300)
plt.clf()


#___________E____________
#fit
a = 2*np.pi/360
def cos(x, b, c):
    return b*np.cos(a*x) + c

params_ao, params_cov_ao = optimize.curve_fit(cos, Kippung["Winkel[rad]"], Kippung["I_K,ano[mA]"], p0 = [0.01, 65])
params_o, params_cov_o = optimize.curve_fit(cos, Kippung["Winkel[rad]"], Kippung["I_K,o[yA]"], p0 = [0.01, 65])


theta_err = 3

fig, ax1 = plt.subplots()
color = 'black'
ax1.set_xlabel('$\Theta [°]$' , fontsize = 13)
ax1.set_ylabel('I [$mA$]', fontsize = 13, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()
color = 'black'
ax2.set_ylabel(r'$I [\mu A$]', fontsize = 13, color=color)  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor=color)

ax1.errorbar(Kippung["Winkel[rad]"], Kippung["I_K,ano[mA]"], xerr = theta_err, yerr = I_err,\
            linewidth = 2, color = "blue", capsize=3, fmt = "x")
ax1.errorbar(np.linspace(0, 95, 1000), cos(np.linspace(0, 95, 1000), params_ao[0], params_ao[1]), color = "black", linestyle = "--")
ax2.errorbar(Kippung["Winkel[rad]"], Kippung["I_K,o[yA]"], xerr = theta_err, yerr = I_err,\
            linewidth = 2, color = "green", capsize=3, fmt = "x")
ax2.errorbar(np.linspace(0, 90, 1000), cos(np.linspace(0, 90, 1000), params_o[0], params_o[1]), color = "black", linestyle = "-.")

ax1.legend(['data non-organic', "fit cos"], loc = 3)
ax2.legend(["data organic", "fit cos"], loc = 1)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.subplots_adjust(top=0.85)

plt.title('Angle dependence of $I_{sc}$', fontsize = 15)


plt.grid(True)
#plt.legend(['U-I Characteristics', "Linearer Fit", "Plot  $U' = U - R_s I$", "logarithmischer Fit"], fontsize = 13)
plt.savefig('Plots/Kippung.png', dpi=300)
plt.clf()
