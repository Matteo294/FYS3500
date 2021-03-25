import numpy as np 

p1 = [0.814, -13.810, -8.409]
p2 = [-27.649, -1.511, -20.665]
p3 = [38.632, 13.361, 44.775]
p4 = [-13.443, 2.514, -5.853]

m_Z = 91.188 #GeV
m_e = 5.1e-4 #Gev

def norm3(x):
    return np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

def Energy(p, m):
    return np.sqrt(m**2+ norm3(p)**2)

en_Z = Energy(p2, m_e) + Energy(p4, m_e)
print("Momentum of Z:", en_Z**2 - m_Z**2)

Etot = np.sqrt( (Energy(p1, m_e) + Energy(p2, m_e) + Energy(p3, m_e) + Energy(p4, m_e))**2 - norm3(p1+p2+p3+p4)**2 )

print("Higgs boson mass:", Etot)