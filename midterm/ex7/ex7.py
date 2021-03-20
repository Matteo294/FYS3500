from matplotlib import pyplot as plt 
import numpy as np
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

T_halve_Cs = 9.27/60 # Half life time of 139Cs in hours
T_halve_Ba = 82.93/60 # Half life time of 139Ba int hours
lambda_Cs = np.log(2)/T_halve_Cs # Decay width of 139Cs
lambda_Ba = np.log(2)/T_halve_Ba # Decay width of 139Ba

print(np.log(2), T_halve_Cs, lambda_Cs, lambda_Ba, lambda_Cs-lambda_Ba)
N0 = 1000

# First decay
def N1(N0, t, lambda1):
    return N0*np.exp(-lambda1*t)

# Second decay
def N2(N0, t, lambda1, lambda2):
    return N0 * lambda1 * (np.exp(-lambda1*t) - np.exp(-lambda2*t)) / (lambda2 - lambda1)

# Third decay
def N3(N0, t, lambda1, lambda2):
    print((lambda1 - lambda2))
    return N0 * (lambda1*(1-np.exp(-lambda2*t)) - lambda2*(1-np.exp(-lambda1*t))) / (lambda1 - lambda2)

# Time values to evaluate the numbers of species
time_ticks = np.linspace(0, 12, 10000)

N_Cs = N1(N0, time_ticks, lambda_Cs)
N_Ba = N2(N0, time_ticks, lambda_Cs, lambda_Ba)
N_La = N3(N0, time_ticks, lambda_Cs, lambda_Ba)

# Plotting
plt.plot(time_ticks, N_Cs, label='$^{139}$Cs')
plt.plot(time_ticks, N_Ba, label='$^{139}$Ba')
plt.plot(time_ticks, N_La, label='$^{139}$La')

plt.xlabel('Time [h]', fontsize=14)
plt.ylabel('N', fontsize=14)
plt.legend(fontsize=14)

plt.show()