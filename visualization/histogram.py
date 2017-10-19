import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

filename = "/home/amirhossein/research/test/energy.dat"

data = np.loadtxt(filename)

hist, bin_edges = np.histogram(data,bins=20)

sns.set_style("white")
plt.ion()
fig = plt.figure()
ax = fig.gca()

ax.hist(data, bins=20, edgecolor="w")

input("Press Enter to exit...")