import numpy as np
import matplotlib.pyplot as plt

# define physical constants
q0 = 1.6e-19

box_size = 50
steady_state_reached = 800

# directory = "/Users/amirhossein/Desktop/new_runs/test_1/"
directory = "/Users/amirhossein/research/test/"

################################################################################
# read the current data
################################################################################
filename = directory + "region_current.dat"
current = np.loadtxt(filename)
time = current[:,0]
current = current[:,1:]

# calculate and plot box average of the data
box = np.ones(box_size)
smooth_current = np.zeros((current.shape[0]-box.shape[0]+1,current.shape[1]))
for i in range(0,current.shape[1]):
	smooth_current[:,i] = np.convolve(current[:,i], box, mode='valid')/np.sum(box)

plt.figure()
plt.plot(smooth_current, linewidth=3)
plt.title("smoothed current density vs time step")
plt.draw()

avg_current = np.mean(current[steady_state_reached:,:],0)
avg_current = np.mean([abs(avg_current[0]), abs(avg_current[-1])])

print("current density[Coulomb*m^-2*second^-1] = %0.2e" % (avg_current))


################################################################################
# read the population profile data
################################################################################
filename = directory + "population_profile.dat"
population = np.loadtxt(filename)
distance = population[0,1:]
time = population[1:,0]
population = population[1:,1:]

box = np.ones(box_size)
smooth_population = np.zeros((population.shape[0]-box.shape[0]+1,population.shape[1]))
for i in range(0,population.shape[1]):
	smooth_population[:,i] = np.convolve(population[:,i], box, mode='valid')/np.sum(box)

# plt.plot(time, data, linewidth=3)
plt.figure()
plt.plot(smooth_population, linewidth=3)
plt.title("smoothed carrier density profile vs time step")
plt.draw()

avg_population = np.mean(population[steady_state_reached:,:],0)
plt.figure()
plt.plot(distance, avg_population, linewidth=3, marker='o')
plt.title("average population density vs distance")
plt.draw()

grad_population = abs((avg_population[0]-avg_population[-1])/(distance[-1]-distance[0]))
print("gradient of population density [m^-4] = %0.2e" %(grad_population))
print("diffusion coefficient [m^2/second] = %0.2e" %(avg_current/q0/grad_population))

plt.show()
