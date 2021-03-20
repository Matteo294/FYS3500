from matplotlib import pyplot as plt 
import numpy as np
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

T_halve_Cs = 9.27/60 # Half life time of 139Cs in hours
T_halve_Ba = 82.93/60 # Half life time of 139Ba int hours
lam_Cs = np.log(2)/T_halve_Cs # Decay width of 139Cs
lam_Ba = np.log(2)/T_halve_Ba # Decay width of 139Ba

N0 = 37e6/(lam_Cs/60/60)

print("N0:", int(N0))

# First decay
def N1(N0, t, lam1):
    return N0*np.exp(-lam1*t)

# Second decay
def N2(N0, t, lam1, lam2):
    return N0 * lam1 * (np.exp(-lam1*t) - np.exp(-lam2*t)) / (lam2 - lam1)

# Third decay
def N3(N0, t, lam1, lam2):
    return N0 * (lam1*(1-np.exp(-lam2*t)) - lam2*(1-np.exp(-lam1*t))) / (lam1 - lam2)

# Activity of the first reaction
def A1(N0, t, lam1):
    return lam1*N0*np.exp(-lam1*t)

def A2(N0, t, lam1, lam2):
    return -N0 * lam1/(lam2-lam1) * (lam2*np.exp(-lam2*t) - lam1*np.exp(-lam1*t))

def A3(N0, t, lam1, lam2):
    return -N0 * lam1*lam2/(lam1-lam2) * (np.exp(-lam2*t) - np.exp(-lam1*t))

# Time values to evaluate the numbers of species
time_ticks = np.linspace(0, 12, 10000)

# N(t)
N_Cs = N1(N0, time_ticks, lam_Cs)
N_Ba = N2(N0, time_ticks, lam_Cs, lam_Ba)
N_La = N3(N0, time_ticks, lam_Cs, lam_Ba)

# Activities
A_Cs = A1(N0, time_ticks, lam_Cs)
A_Ba = A2(N0, time_ticks, lam_Cs, lam_Ba)
A_La = A3(N0, time_ticks, lam_Cs, lam_Ba)

# Finding requested quantities
max_A_Ba = print("Maximum activity: %.3f at time %.3f hours" % (max(A_Ba), time_ticks[np.where(A_Ba == max(A_Ba))[0][0]]))
diff = np.abs(A_Ba - A_Cs)
print("Activities of Cs and Ba are equal at time %.3f hours" % (time_ticks[np.where(diff == min(diff))]))

# Plotting N(t)
plt.plot(time_ticks, N_Cs, label='$^{139}$Cs')
plt.plot(time_ticks, N_Ba, label='$^{139}$Ba')
plt.plot(time_ticks, N_La, label='$^{139}$La')
plt.xlabel('Time [h]', fontsize=14)
plt.ylabel('N', fontsize=14)
plt.legend(fontsize=14)
plt.show()

# Plotting activities
plt.plot(time_ticks, A_Cs, label='$^{139}$Cs')
plt.plot(time_ticks, A_Ba, label='$^{139}$Ba')
plt.plot(time_ticks, A_La, label='$^{139}$La')
plt.xlabel('Time [h]', fontsize=14)
plt.ylabel('Activity', fontsize=14)
plt.legend(fontsize=14)
plt.show()