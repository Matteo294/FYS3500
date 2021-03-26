from matplotlib import pyplot as plt 
import numpy as np 

# LSQ Coefficients (source Wikipedia)
a_v = 15.8 
a_s = 18.3 
a_c = 0.714 
a_symm = 23.2
a_p = 12

# Other constants
mp = 938.272 # Proton mass in MeV/c^2
me = 0.511 # Electron mass in MeV/c^2
mn = 939.565 # Neutron mass in MeV/c^2
MeV_to_u = 931.5 # conv factor from MeV/c^2 to atomic units (divide)

def B(A, Z):
    delta = 0
    if A%2 == 0 and Z%2 == 0:
        delta = a_p/np.sqrt(A)
    elif A%2 == 0 and Z%2 == 1:
        delta = -a_p/np.sqrt(A)
    return + a_v*A-a_s*A**(2/3)-a_c*Z*(Z-1)/A**(1/3)-a_symm*(A-2*Z)**2/A+delta

def Coulomb(r, a=1):
    return a/r

def protonPotential(r, a=1):
    if r<0.5:
        return 100
    elif r<5:
        return Coulomb(r, a) - 10
    else:
        return Coulomb(r, a)

# Repulsion vertical line
repulsion = (1000, -10)
xrepulsion = (0.5, 0.5)
# Well horizontal line
well = (-10, -10)
xwell = (0.5, 5)
# Attraction vertical line
attraction = (-10, 0)
xattraction = (5, 5)
# Horizonal line
horiz = (0, 0)
xhoriz = (5, 100)
# Coulomb curve
xline = np.linspace(0.01, 100, 1000)
yline = Coulomb(xline, a=10)
# Proton potential
zline = [protonPotential(xx, a=10) for xx in xline]

plt.plot(xrepulsion, repulsion, color='dodgerblue', 
        label="Nuclear potential", linewidth=2.0)
plt.plot(xwell, well, color='dodgerblue', linewidth=2.0)
plt.plot(xattraction, attraction, color='dodgerblue', linewidth=2.0)
plt.plot(xhoriz, horiz, color='dodgerblue', linewidth=2.0)
plt.plot(xline, yline, color='sandybrown', 
            label='Coulomb potential', linewidth=2.0)
plt.plot(xline, zline, '--', color='gray', 
            label='Proton interaction', linewidth=2.0)
plt.legend(fontsize=14)
plt.xlim(0, 15)
plt.ylim(-15, 30)
plt.xlabel('Distance', fontsize=16)
plt.ylabel('Potential', fontsize=16)
plt.xticks([], [])
plt.yticks([], [])
plt.show()

# 1e) neutron energy
print("118 Sn:", B(118, 50) - B(117, 50))
print("16 O:", B(16, 8) - B(15, 8))
print("119 Sn:", B(119, 50) - B(118, 50))
