import numpy as np
import matplotlib.pyplot as plt
import os.path
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.signal import filtfilt
import argparse

import plotly.offline as plo
import plotly.graph_objs as go

# from bokeh.plotting import figure, output_file, show
import bokeh.plotting as bp

# customize matplotlib styles
mpl.rc('lines', linewidth=4) # default linewidth
mpl.rc('lines', dash_capstyle='round') # default dashed line capstyle
mpl.rc('lines', solid_capstyle='round') # default solid line capstyle
mpl.rc('xtick', labelsize=10) # default label size for xtick
mpl.rc('ytick', labelsize=10) # default label size for ytick
mpl.rc('axes', titlesize=13) # default title size

# # Physical constants
q0 = 1.6e-19  # electron charge in Coulomb
eV = 1.6e-19  # electron volt in Jouls

def plot_current(directory, box_size=200, color=None, ax=None):
  '''
  Plot current data
  '''
  if (ax is None):
      fig = plt.figure(figsize=(12,10))
      ax = fig.add_subplot(1,1,1)
  filename = os.path.join(directory, "region_current.dat")
  df = pd.read_csv(filename,header=None, delim_whitespace=True)
  time = df.iloc[:,0].values
  current = df.iloc[:,1:].values
  
  steady_current = np.average(current)
  print(f"steady state drain current: {steady_current:.6e}")

  assert (current.shape[0] > box_size), f"box_size must be smaller than the number of time steps: {current.shape[0]}"
  current = filtfilt(np.ones(box_size),box_size,current, axis=0)
  time = time[:current.shape[0]]
  
  ax.plot(time*1.e9,current, color=color)
  ax.set_title("Particle current vs time")
  ax.set_xlabel("Time [ns]")
  ax.set_ylabel("Current [#/sec.m$^2$]")
  return steady_current

def plot_population_vs_time(directory, box_size=10, color=None, ax=None):
  '''
  Plot population versus time at different cross sections
  '''
  filename = os.path.join(directory, "population_profile.dat")
  df = pd.read_csv(filename, header=None, sep=' ')
  time = df.iloc[:, 0].values
  population = df.iloc[:, 1:11].values

  assert (
      population.shape[0] > box_size), f"box_size must be smaller than the number of time steps: {population.shape[0]}"

  population = filtfilt(np.ones(box_size), box_size, population, axis=0)
  time = time[:population.shape[0]]

  if (ax is None):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
  ax.plot(time*1.e9, population, linewidth=4, color=color)
  _ = ax.set_title("Particle density vs time")
  _ = ax.set_xlabel("Time [ns]")
  _ = ax.set_ylabel("Population [#/m$^3^]")
  return ax

def plot_population_vs_position(directory, start=0, end=-1, color=None, ax=None):
  '''
  plot population density versus position in the simulation domain
  '''
  filename = os.path.join(directory, "population_profile.dat")
  df = pd.read_csv(filename, header=None, sep=' ')
  population = df.iloc[:, 1:11].values

  filename = os.path.join(directory, "scatterer_distribution.dat")
  df = pd.read_csv(filename, delim_whitespace=True)
  distance = df['position'].values
  density = df['density'].values

  avg_population = np.mean(population[start:end, :], 0)
  avg_population = np.divide(avg_population, density)

  if (ax is None):
      fig = plt.figure()
      ax = fig.add_subplot(1, 1, 1)
  ax.plot(distance*1.e9, avg_population, marker='o', markersize=10, color=color)
  ax.set_title("Particle profile vs position")
  ax.set_xlabel("Position [nm]")
  ax.set_ylabel("Population [A.U.]")

  dn = avg_population[0]-avg_population[-1]
  dx = distance[0]-distance[-1]
  density_gradient = dn/dx
  print(f"dn/dx = {density_gradient:0.2e}")
  return density_gradient

def plot_scatterer_stats(directory):
  '''
  plot distribution of scatterers versus position
  '''
  filename = os.path.join(directory, "scatterer_statistics.dat")
  df = pd.read_csv(filename)
  print(df)

  fig = plt.figure()
  ax1 = fig.add_subplot(3, 1, 1)
  ax1.plot(df['position'], df['distribution'])
  ax1.set_ylim(0)
  ax2 = fig.add_subplot(3, 1, 2)
  ax2.plot(df['position'], df['population'])
  ax2.set_ylim(0)
  ax3 = fig.add_subplot(3, 1, 3)
  ax3.plot(df['position'], df['density'])
  ax3.set_ylim(0)

  ax1.set_title("Scatterer distribution vs position")
  ax3.set_xlabel("position [m]")

  plt.show()

def read_cnt_coordinates(directory=None):
  '''
  Read the cnt coordinates and return the coordinates in a single np array
  r[dim,num_tube,num_point_in_tube]
  '''

  if directory is None:
    directory = os.path.expanduser("~/research/cnt_mesh_fiber")

  print('reading cnt coordinates...', end='', flush=True)

  filename = os.path.join(directory, 'single_cnt.pos.x.dat')
  x = np.loadtxt(filename, skiprows=2)
  filename = os.path.join(directory, 'single_cnt.pos.y.dat')
  y = np.loadtxt(filename, skiprows=2)
  filename = os.path.join(directory, 'single_cnt.pos.z.dat')
  z = np.loadtxt(filename, skiprows=2)

  r = np.stack((x, y, z))

  print('done!')
  return r

def plot_mesh():
  directory = os.path.expanduser("~/research/cnt_mesh_fiber")
  r = read_cnt_coordinates(directory)

  fig = plt.figure()
  ax = fig.add_subplot('111',projection='3d')
  ax.plot(r[0, 0, :], r[1, 0, :], r[2, 0, :])
  plt.show()

def histogram(r):
  r = r.reshape((3,-1))
  mincoor = np.min(r, axis=1)
  maxcoor = np.max(r, axis=1)
  print(f'''
         min coordinates: {mincoor}
         max coordinates: {maxcoor}
         ''')
  
  
  n = np.array([10, 10, 10])
  dr = (maxcoor - mincoor)/n

  r = r.T
  r = r - mincoor
  
  idx = r / dr
  idx = np.clip(idx, 0, n-1)

  idx = idx.astype(dtype='int')

  print(idx)

  return idx

  # def get_index(coor, mincoor, )

  # return r

  xmin, xmax = np.min(r[0, :]), np.max(r[0, :])
  ymin, ymax = np.min(r[1, :]), np.max(r[1, :])
  zmin, zmax = np.min(r[2, :]), np.max(r[2, :])

  nx, ny, nz = 10, 10, 10
  dx, dy, dz = (xmax-xmin)/nx, (ymax-ymin)/ny, (zmax-zmin)/nz

  mincoor = np.array([xmin, ymin, zmin])
  dr = np.array([dx, dy, dz])

  hist = np.zeros((nx, ny, nz), dtype='int')

  for i in range(r.shape[1]):
    ix = int((r[0, i]-xmin)/dx)
    iy = int((r[1, i]-ymin)/dy)
    iz = int((r[2, i]-zmin)/dz)
    ix = np.clip(ix, 0, nx-1)
    iy = np.clip(iy, 0, ny-1)
    iz = np.clip(iz, 0, nz-1)
    hist[ix, iy, iz] += 1
  
  return hist


if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--rates', help='plot scattering table rates', action='store_true')
  parser.add_argument('--histogram', help='plot histogram of populations', action='store_true')
  parser.add_argument('--all', help='plot all the results together', action='store_true')
  parser.add_argument('--coordinates', help='plot atom coordinates for donor and acceptor CNTs', action='store_true')
  parser.add_argument('--bokeh', help='test plotting with bokeh', action='store_true')

  args = parser.parse_args()

  if args.histogram:
    directory = os.path.expanduser("~/research/cnt_mesh_fiber")
    r = read_cnt_coordinates(directory)
    print(f'number of cnt coordinates: {r.shape}', flush=True)
    h = histogram(r)
    print(h.shape)

  
  if args.rates:
    d = os.path.expanduser("~/research/monte_carlo_fiber_test")
    filename = os.path.join(d,'scat_table.rates.dat')
    rates = np.loadtxt(filename, skiprows=4)
    file = open(filename)
    lines = file.readlines()
    dims = lines[2]
    dims = dims[:dims.find('\n')].split(',')
    dims = [int(d) for d in dims]
    rates = rates.reshape(dims, order='C')

    filename = os.path.join(d, 'scat_table.theta.dat')
    theta = np.loadtxt(filename)
    
    filename = os.path.join(d, 'scat_table.z_shift.dat')
    zshift = np.loadtxt(filename)

    filename = os.path.join(d, 'scat_table.axis_shift_1.dat')
    ax_shift_1 = np.loadtxt(filename)

    filename = os.path.join(d, 'scat_table.axis_shift_2.dat')
    ax_shift_2 = np.loadtxt(filename)

    print("plotting hopping rate versus angle:")

    i_ax_shift_1 = np.where(ax_shift_1 == 0)
    i_ax_shift_2 = np.where(ax_shift_2 == 0)

    print(f"axis shift 1: {ax_shift_1[i_ax_shift_1]}")
    print(f"axis shift 2: {ax_shift_2[i_ax_shift_2]}")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i_z_shift, zsh in enumerate(zshift):
      rates_vs_theta = rates[:, i_z_shift, i_ax_shift_1, i_ax_shift_2].reshape((-1))
      ax.plot(theta/np.pi, rates_vs_theta, linewidth=3, label=f"{zsh*1.e9:.2f} nm")
    ax.legend(title="z-shift")
    ax.set_title("Hopping rate versus angle")
    ax.set_ylabel("Rate [sec$^{-1}$]")
    ax.set_xlabel("Angle [$\\pi$]")
    # ax.set_yscale('log')

    print("plotting hopping rate versus z-shift:")

    i_ax_shift_1 = np.where(ax_shift_1 == 0)
    i_ax_shift_2 = np.where(ax_shift_2 == 0)

    print(f"axis shift 1: {ax_shift_1[i_ax_shift_1]}")
    print(f"axis shift 2: {ax_shift_2[i_ax_shift_2]}")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for i_theta, th in enumerate(theta):
      rates_vs_zshift = rates[i_theta, :, i_ax_shift_1, i_ax_shift_2].reshape((-1))
      ax.plot(zshift*1e9, rates_vs_zshift, linewidth=3, label=f"{th/np.pi:.2f}")
    ax.legend(title="Angle [$\\pi$]")
    ax.set_title("Hopping rate versus z-shift")
    ax.set_ylabel("Rate [sec$^{-1}$]")
    ax.set_xlabel("Distance [nm]")
    # ax.set_yscale('log')

    print("plotting hopping rate versus angle and z-shift:")

    i_ax_shift_1 = np.where(ax_shift_1 == 0)
    i_ax_shift_2 = np.where(ax_shift_2 == 0)

    print(f"axis shift 1: {ax_shift_1[i_ax_shift_1]}")
    print(f"axis shift 2: {ax_shift_2[i_ax_shift_2]}")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    rates_vs_theta_zshift = rates[:, :, i_ax_shift_1, i_ax_shift_2].reshape((len(theta), len(zshift)))
    ax.imshow(np.log(rates_vs_theta_zshift), extent=(np.min(zshift)*1e9, np.max(zshift)*1e9, np.min(theta)/np.pi, np.max(theta)/np.pi), aspect='auto', origin='lower', interpolation='bilinear')
    ax.set_title("Hopping rate versus angle and z-shift")
    ax.set_xlabel("z-shift [nm]")
    ax.set_ylabel("angle [$\\pi$]")

    plt.show() 

  if args.coordinates:
    directory = os.path.expanduser('~/research/monte_carlo_fiber_test')
    
    filename = os.path.join(directory, '3d_coord.180_angle.cnt1.dat')
    # filename = os.path.join(directory, '3d_coord.cnt1.dat')
    r1 = np.loadtxt(filename)

    filename = os.path.join(directory, '3d_coord.180_angle.cnt2.dat')
    # filename = os.path.join(directory, '3d_coord.cnt2.dat')
    r2 = np.loadtxt(filename)

    fig = plt.figure()
    ax = fig.add_subplot('111',projection='3d')
    
    ax.scatter(r1[:, 0], r1[:, 1], r1[:, 2], marker='o')
    ax.scatter(r2[:, 0], r2[:, 1], r2[:, 2], marker='o')

    # Create cubic bounding box to simulate equal aspect ratio
    ax.set_aspect('equal')
    all_r = np.concatenate((r1, r2), axis=0)
    max_range = (np.max(all_r,axis=0)-np.min(all_r,axis=0)).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5*(all_r[:, 0].max()+all_r[:, 0].min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5*(all_r[:, 1].max()+all_r[:, 1].min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5*(all_r[:, 2].max()+all_r[:, 2].min())

    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
      ax.plot([xb], [yb], [zb], 'w')

    plt.show()

    





  # # Give the input directories as a list
  directories = []
  directories.append(os.path.expanduser("~/research/monte_carlo_fiber_test"))
  

  if args.all:

    plot_scatterer_stats(directories[0])

    current = []

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    legends = []
    for d in directories:
        c = plot_current(d, ax=ax, box_size=1000)
        current.append(c)
        legends.append(os.path.basename(d))

    current = np.array(current)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    _ = plot_population_vs_time(directories[0],ax=ax, box_size=1000)
    plt.show()

    density_gradient = []

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for d in directories:
      grad = plot_population_vs_position(d, start=1,end=-1,ax=ax)
      density_gradient.append(grad)
    _ = ax.set_title("Particle profile vs position at various time intervals")

    density_gradient = np.array(density_gradient)
    plt.show()


if args.bokeh:
  # prepare some data
  x = [1, 2, 3, 4, 5]
  y = [6, 7, 2, 4, 5]

  filename = os.path.expanduser('~/Desktop/lines.html')

  # output to static HTML file

  bp.output_file(filename)

  # create a new plot with a title and axis labels
  p = bp.figure(title="simple line example", x_axis_label='x', y_axis_label='y')

  # add a line renderer with legend and line thickness
  # p.line(x, y, legend="Temp.", line_width=2)
  p.circle(x,y)

  # show the results
  bp.show(p)

# # # Calculate diffusion coefficient
# # $$ D = J/\big(\frac{dn}{dx}\big) $$
# # has units of [$\text{m}^2/\text{s}$]

# # In[49]:


# diff_coef = np.divide(current,density_gradient)
# x_axis = np.array(range(len(diff_coef)))
# # x_axis = np.array([5,10,15,20])

# fig = plt.figure(figsize=(10,9))
# ax = fig.add_subplot(1,1,1)
# ax.plot(x_axis, diff_coef[:len(x_axis)]*1.e4, linewidth=4, marker='o', markersize='14')

# ax.axhline(y=diff_coef[-1]*1.e4, linewidth=4, color='red', linestyle='dashed')

# ax.set_title("Diffusion coefficient", fontsize=30)
# ax.set_ylabel("Diffusion coeffitient [cm$^2$/s]",fontsize=24)
# ax.set_xlabel("Length of CNT section [nm]",fontsize=24)
# ax.tick_params(labelsize=20)
# ax.set_xlim([4.8,20.2])
# ax.legend(['fixed section length', 'random section length\nbetween 5 nm and 20 nm'], fontsize=24)


# # # Plot particle path

# # In[51]:


# # directories = [os.path.expanduser("~/research/monte_carlo_fiber")]

# fig = plt.figure(figsize=(20,10))
# ax = fig.add_subplot(111,projection='3d')
# data = []
# for directory in directories:
#     for i in range(11,20):
#         filename = os.path.join(directory,"particle_path."+str(i)+".dat")
#         if not os.path.exists(filename):
#             continue
#         x,y,z = np.loadtxt(filename, unpack=True)
#         ax.plot(z, x, y)
#         ax.set_aspect('equal')

#         trace = go.Scatter3d(
#             x=z, y=x, z=y,
#             marker=dict(
#                 size=0,
#                 colorscale='Viridis',
#             ),
#             line=dict(
#                 width=1
#             )
#         )

#         data.append(trace)


# layout = dict(
#     autosize=True,
#     title='Iris dataset',
#     scene=dict(
#         xaxis=dict(
#             gridcolor='rgb(255, 255, 255)',
#             zerolinecolor='rgb(255, 255, 255)',
#             showbackground=True,
#             backgroundcolor='rgb(230, 230,230)'
#         ),
#         yaxis=dict(
#             gridcolor='rgb(255, 255, 255)',
#             zerolinecolor='rgb(255, 255, 255)',
#             showbackground=True,
#             backgroundcolor='rgb(230, 230,230)'
#         ),
#         zaxis=dict(
#             gridcolor='rgb(255, 255, 255)',
#             zerolinecolor='rgb(255, 255, 255)',
#             showbackground=True,
#             backgroundcolor='rgb(230, 230,230)'
#         ),
#         camera=dict(
#             up=dict(
#                 x=0,
#                 y=0,
#                 z=1
#             ),
#             eye=dict(
#                 x=-1.7428,
#                 y=1.0707,
#                 z=0.7100,
#             )
#         ),
#         aspectratio = dict( x=1, y=1, z=0.03 ),
#         aspectmode = 'manual'
#     ),
# )

# fig = go.Figure(data=data, layout=layout)
# plo.plot(fig,filename=os.path.join(directory,"temp-pyplot.html"))


# # # Perform statistical analysis of particle path

# # In[22]:


# directories = [os.path.expanduser("~/research/monte_carlo_fiber_track_particle.1")]
# directories += [os.path.expanduser("~/research/monte_carlo_fiber_track_particle.2")]
# directories += [os.path.expanduser("~/research/monte_carlo_fiber_track_particle.3")]

# for directory in directories:
#     print(directory)
#     xmean, ymean, zmean = [], [], []
#     std = {'x':[],'y':[],'z':[]}
#     for file in os.listdir(directory):
#         if file.startswith("particle_path"):
#             fullpath = os.path.join(directory,file)
#             x,y,z = np.loadtxt(fullpath, unpack=True)
#             xmean += np.mean(x)
#             ymean += np.mean(y)
#             zmean += np.mean(z)
#             std['x'].append(np.std(x))
#             std['y'].append(np.std(y))
#             std['z'].append(np.std(z))
            
#     std['x'] = np.array(std['x'])
#     std['y'] = np.array(std['y'])
#     std['z'] = np.array(std['z'])
    
#     print("standard deviation is [{:1.2e},{:1.2e},{:1.2e}]".format(np.mean(std['x']), np.mean(std['y']), np.mean(std['z'])))

