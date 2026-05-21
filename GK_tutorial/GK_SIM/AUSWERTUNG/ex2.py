import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

data=np.genfromtxt(Path(__file__).parent / 'histo.dat')

particle_numbers=data[:,1].astype(int)
# plt.hist(particle_numbers,bins=np.arange(0,np.max(particle_numbers)))
# plt.show()

# generate unweighted histogram
weights = data[:,14]
w_uw=np.exp(weights)
weighted_c = np.bincount(
particle_numbers,
weights=w_uw
)

P = weighted_c/np.sum(weighted_c)
N=np.arange(0,np.max(particle_numbers)+1)
plt.bar(N,P)
#plt.show()

# calculate free energy
# F = k_B*T*ln(P)
F = - 0.9* np.log(P)
plt.plot(N,F)
plt.show()

# Determine the average energy per particle (fourth column in
# histo.dat) in the gas and the liquid phase. Explain the qualitative
# difference between the two.

energy=data[:,3]

average_energy_per_particle = energy/particle_numbers

n_divisor=np.mean(N)
avg_E_pp=[]
for i in N:
    if np.sum((particle_numbers==i))>0:
        avg_E_pp.append(np.mean(energy[particle_numbers==i]/i))
        
    else:
        avg_E_pp.append(0)

avg_E_pp=np.array(avg_E_pp)

gas_mask = N < n_divisor
liquid_mask = ~gas_mask
gas_particle_energy = np.sum(P[gas_mask]*avg_E_pp[gas_mask])/np.sum(P[gas_mask])

liquid_particle_energy = np.sum(P[liquid_mask]*avg_E_pp[liquid_mask])/np.sum(P[liquid_mask])

print(gas_particle_energy, liquid_particle_energy)

# gas_particle_energy is much bigger than liquid_particle_energy, because in the gas phase, particles are more energetic and have more freedom to move around, while in the liquid phase, particles are more closely packed and have less energy due to interactions with neighboring particles.

# Exercise 3
