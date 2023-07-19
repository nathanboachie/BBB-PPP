#File:		video_generator.py
#Date:		19.07.2023
#Author:	Armand Sieber
#Tag:		Short script to visualize simualtion output of bubble dynamics
#Copyright:	2023 EPFL. All rights reserved.


import numpy as np
import matplotlib.pyplot as plt


# ---- User-defined inputs ----
# -----------------------------

# Must be provided by user
file_name = 'unbounded_Rayleigh_Plesset_epsilon100.txt' # simulation output filename
gamma = 0.56 # stand-off distance
Nb = 40 # number of elements on bubble surface
Ns = 60 # number of elements on fluid-fluid interface
time_step = -1 # time step of interest


# ---- Loads  solution vectors ----
# ---------------------------------
coordinates = np.genfromtxt(file_name, skip_header = 10, skip_footer=1)

time = coordinates[:, 1] # time vector

r_b = coordinates[:, 2:Nb+3] # r-coordinates bubble surface
z_b = coordinates[:, Nb+3:2*Nb+4] #z-coordinates bubble surface

r_s = coordinates[:, 4*Nb+6:4*Nb+6+Ns+1] # r-coordinates fluid-fluid interface
z_s = coordinates[:, 4*Nb+6+Ns+1:4*Nb+6+2*Ns+2] # z-coordinates fluid-fluid interface

volume = coordinates[:, -1] # bubble volume

# ---- Display solutions ----
# ---------------------------

# bubble volume and equivalent radius temporal evolution
# ------------------------------------------------------

pi = 3.1415926535897932
radius = (3.0 / 4.0 * volume / pi)**(1.0 / 3.0)

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
plt.grid()
ax1.plot(time, volume, '-k', linewidth = 2, label = 'Volume')
ax1.set_ylabel('volume [-]',fontsize=12, color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xlabel('time [-]', fontsize=12)

ax2 = ax1.twinx()
ax2.plot(time, radius, '-r', linewidth = 2, label = 'Equ. radius')
ax2.set_ylabel('Equ. radius [-]',fontsize=12, color = 'red')
ax2.tick_params(axis='y', labelcolor='red')

# bubble and fluid-fluid interface shape
# --------------------------------------

# bubble
rb = r_b[time_step, :]
zb = z_b[time_step, :]
rb_flip = np.flipud(rb)
zb_flip = np.flipud(zb)
rb_tot = np.concatenate((rb, -rb_flip))
zb_tot = np.concatenate((zb, zb_flip))

# fluid-fluid interface
rs = r_s[-1, :]
zs = z_s[-1, :]
rs_flip = np.flipud(rs)
zs_flip = np.flipud(zs)
rs_tot = np.concatenate((rs, -rs_flip))
zs_tot = np.concatenate((zs, zs_flip))

# Region of interest
ymin = - gamma - 2.0
ymax = - gamma + 2.0
xmin = -2.0
xmax = 2.0

fig = plt.figure(2)
plt.axis('equal')
plt.ylim(ymin, ymax)
plt.xlim(xmin, xmax)
plt.grid()
plt.plot(rb_tot, zb_tot, '-r', linewidth = 2)
plt.plot(rs_tot, zs_tot, '-b', linewidth = 2)
plt.xlabel(r'$r$',fontsize=12)
plt.ylabel(r'$z$',fontsize=12)


plt.show()
