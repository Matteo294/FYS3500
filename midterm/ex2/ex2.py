import numpy as np 
from matplotlib import pyplot as plt 

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
MeV_to_u = 931.5 # conversion factor from MeV/c^2 to atomic mass units (divide)


# Mass semi-empirical formula
def M(A, Z, force_even_Z=False, force_odd_Z=False):
    return Z*(mp+me) + (A-Z)*mn - B(A, Z, force_even_Z=force_even_Z, force_odd_Z=force_odd_Z) 

# Binding energy
def B(A, Z, force_even_Z=False, force_odd_Z=False):
    delta = 0
    if (A%2 == 0 and Z%2 == 0 and not force_odd_Z) or force_even_Z:
        delta = a_p/np.sqrt(A)
    elif (A%2 == 0 and Z%2 == 1) or force_odd_Z:
        delta = -a_p/np.sqrt(A)
    return a_v*A - a_s*A**(2/3) - a_c*Z*(Z-1)/A**(1/3) - a_symm*(A-2*Z)**2/A + delta

# Find most stable atom for fixed A
def most_stable(A):
    Zbest_rounded = int((mn - mp -me + 4*A*a_symm + a_c*A**(2/3)) / (2*a_c*A**(2/3) + 8*a_symm))
    if M(A, Zbest_rounded) > M(A, Zbest_rounded  + 1):
        return Zbest_rounded + 1
    else:
        return Zbest_rounded

print("Mass of 48Ca:", M(48, 20)/MeV_to_u)
print("Energy to extract a neutron from 44Ca:", (B(44, 20) - B(43, 20)))
print("Most stable Z for fixed A=136:", most_stable(136))

# Z values
Zline = np.linspace(52, 60, 1000)
Zvals = np.arange(52, 61, 1)

# Calculate masses line
Mvals_oddA = M(135, Zline)
Mvals_evenA_evenZ = [M(136, zz, force_even_Z=True) for zz in Zline]
Mvals_evenA_oddZ = [M(136, zz, force_odd_Z=True) for zz in Zline]
# Discrete masses value
M_oddA = [M(135, zz) for zz in Zvals]
M_evenA = [M(136, zz) for zz in Zvals]

# Odd A
plt.plot(Zline, Mvals_oddA)
plt.plot(Zvals, M_oddA, '.')
plt.xlabel('Z', fontsize=16)
plt.ylabel('Mass', fontsize=16)
plt.show()

# Even A
plt.plot(Zline, Mvals_evenA_evenZ, label="Even Z,N", markersize=14)
plt.plot(Zline, Mvals_evenA_oddZ, label="Odd Z,N", markersize=14)
plt.plot(Zvals, M_evenA, '.')
plt.xlabel('Z', fontsize=16)
plt.ylabel('Mass', fontsize=16)
plt.legend(fontsize=12)
plt.show()

