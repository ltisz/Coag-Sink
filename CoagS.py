#Coagulation Sink Calculator
#Lee Tiszenkel - University of Alabama Huntsville

import numpy as np
import re
import os
import sys

#Setting the path of the file to read to the same folder as the script
scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))
os.chdir(scriptdir)

#Constants
k = (1.381*(10**-23))           #Boltzmann constant
pi = 3.1415927                  #Pi
e = 2.718281828459045           #e
lamb = 0.0686                   #Mean free path in air (micrometers)
mu = (1.80*(10**-5))            #Viscosity of air (Pa * s)

#Variables
T = 298                         #Temperature
rho = 1                         #Density (g cm^-3)
h = 1                           #Correction factor in gamma calculation, taken to be 1 with particle measurements

def CoagS(dp1, dp2):
    dp1cm = dp1/10000000                                    
    dp2cm = dp2/10000000
    Kn1 = (2.*lamb)/(dp1/1000.)                              #Knudsen number (unitless)
    Cc1 = 1+(Kn1*(1.142+(0.58*(e**(-0.99/Kn1)))))          #Slip correction (unitless)
    #print "Cc " + str(dp1) + ": " + str(Cc1)
    D1 = ((k*T*Cc1)/(3.*pi*mu*(dp1/1000000000.)))*10000.      #Diffusivity in cm^2 s^-1
    Kn2 = (2*lamb)/(dp2/1000.)                              #Knudsen number (unitless)
    Cc2 = 1+(Kn2*(1.142+(0.58*(e**(-0.99/Kn2)))))          #Slip correction (unitless)
    D2 = ((k*T*Cc2)/(3.*pi*mu*(dp2/1000000000.)))*10000.      #Diffusivity in cm^2 s^-1
    
    Bi = D1/(k*T)                                          #Mechanical Mobilities (Millikan formula, see Tammet 1994)
    Bj = D2/(k*T)
    
    mi = rho * ((pi*((dp1/10000000.)**3.))/6.)/1000.            #Mass of particle with diameter dp1 in kg
    #print "Mass " + str(dp1) + ": " + str(mi)
    mj = rho * ((pi*((dp2/10000000.)**3.))/6.)/1000.            #Mass of particle with diameter dp2 in kg
    #print "Mass " + str(dp1) + ": " + str(mi)
    
    #gamma = 2*((Bi+Bj)/(dp1+dp2+(2*h)))*(((2*pi*k*T*mi*mj)/(mi+mj))**(1/2))     #Gamma term in coag. coefficient calculation
    
    #Kij = (2*pi*k*T*(Bi+Bj)*(dp1+dp2+2*h))/(1+gamma-((0.299*gamma)/((gamma**1.1)+0.64)))    #Coag. coefficient calc.
                                                                                             #Alpha (sticking coefficient) assumed to be 1 (Lee 2016)
    
    cBar1 = (((8.*k*T)/(pi*mi))**(1./2.))*100              #velocity in cm/s
    cBar2 = (((8.*k*T)/(pi*mj))**(1./2.))*100 
    
    zeta1 = ((8.*D1)/(pi*cBar1))
    zeta2 = ((8.*D2)/(pi*cBar2))
    
    g1 = (((2.**(1/2))/(3.*dp1cm*zeta1))*(((dp1cm+zeta1)**3.)-((dp1cm**2. + zeta1**2)**(3./2.))))-dp1cm
    g2 = (((2.**(1/2))/(3.*dp2cm*zeta2))*(((dp2cm+zeta2)**3.)-((dp2cm**2. + zeta2**2)**(3./2.))))-dp2cm
    
    beta = ((dp1+dp2)/(dp1+dp2+(2.*(((g1**2.)+(g2**2.))**(1./2.)))))+((8.*(D1+D2))/((((cBar1**2.)+(cBar2**2.))**(1./2.))*(dp1+dp2)))
    Kij = (2.*pi*(dp1cm+dp2cm)*(D1+D2))*(beta)
    #print "Kij for " + str(dp1) + " and " + str(dp2) + ": " + str(Kij)
    return Kij
    
filename='0921-SMPS-1.txt'
datastartline = 0

##Skipping the SMPS output header:
lookup = 'Start Time' #Column heading lines (first line of file we need)
with open(filename, 'rb') as myFile:
    for num,line in enumerate(myFile, 1):
        if lookup.encode() in line:
            print('found at line:', num)
            dataStartLine = num-1
    
##Getting list of particle sizes from SMPS output, made in to a list "y"
y = np.genfromtxt(filename, skip_header=dataStartLine, comments='Comment', delimiter=',', max_rows=1)       ##Grab column header line
yy = [s for s in list(y) if str(s) != 'nan']            ##Every line in column header line aside from particle sizes is string and comes through as 'nan'
tupleStart = list(y).index(yy[0])                                                                           ##Making our range of columns to grab "Z"
tupleEnd = tupleStart + len(yy)
Z=tuple(range(tupleStart,tupleEnd))
fullDiamList = np.array(yy)                                        ##The last step made a list, so this converts it back to a numpy array for plotting
concs = np.genfromtxt(filename, skip_header=dataStartLine+1, delimiter=',', usecols=Z)
#print concs[0]
print fullDiamList

#print fullDiamList[np.where(fullDiamList==11.8)[0][0]:len(fullDiamList)-1]
#print y

CoagSinkPerLine = []
CoagCoeffList = []

#Getting coagulation coefficients for all diameters in SMPS export
for sample in concs:
    fullCoagSink = []
    for diam in fullDiamList:
        greaterDiams = fullDiamList[np.where(fullDiamList==diam)[0][0]:len(fullDiamList)-1]
        for diam2 in greaterDiams:
            Kij = CoagS(diam, diam2)
            CoagCoeffList.append(Kij)
            #print "K for size " + str(diam) + " and " + str(diam2) + ": " + str(Kij)
            #print "Concentration in size bin " + str(diam2) + ": " + str(sample[np.where(fullDiamList==diam2)])
            if sample[np.where(fullDiamList==diam2)] == 'nan':
                CoagSink = 0
            else:
                CoagSink = Kij * sample[np.where(fullDiamList==diam2)]
                fullCoagSink.append(CoagSink)
    CoagSinkPerLine.append(sum(fullCoagSink))

for item in CoagSinkPerLine:
    print str(item[0]) + ' s^-1'
