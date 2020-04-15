## python script for plotting earthquake building response in a region

import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import scipy.interpolate
from scipy import ndimage
rcParams['figure.figsize'] = 25, 20
rcParams['axes.linewidth'] = 3
rcParams['font.family'] = "sans-serif"
rcParams['font.sans-serif'] = "Arial"

N = 1000 #number of points for plotting/interpolation
x1, y1, z1 = np.genfromtxt(r'/Users/MMiah/Documents/mapping/gf_hypo_2019/fulldomain/analysis/pidratios/left_hypo_pid/40st_fn.txt', unpack=True)
x2, y2, z2 = np.genfromtxt(r'/Users/MMiah/Documents/mapping/gf_hypo_2019/fulldomain/analysis/pidratios/center_hypo_pid/40st_fn.txt', unpack=True)
#x3, y3, z3 = np.genfromtxt(r'/Users/MMiah/Documents/mapping/gf_hypo_2019/fulldomain/analysis/pidratios/right_hypo_pid/40st_fn.txt', unpack=True)

print max(z1/z2)
#Add to a dictionary the ratios of the PIDs and their locations in the domain
dic={k:v for k,v in zip(zip(x1,y1),z1/z2)}
print max(dic, key=dic.get) # Find maximum ratio in the domain
#items = dic.items()
# Declare two lists of elements containing the station indices for the subdomain
l1=[i for i in range(9,92)]  # contains station index along horizontal
l2=[i for i in range(11,21)]  # contains station index along normal
xs=[]
ys=[]
zs=[]
for i in range(len(x1)):
    if x1[i] in l1 and y1[i] in l2:
        xs.append(x1[i])
        ys.append(y1[i])
        zs.append(z1[i]/z2[i])
print max(zs)
#Add to a dictionary the ratios of the PIDs and their locations in the subdomain
subdic={k:v for k,v in zip(zip(xs,ys),zs)}
print max(subdic, key=subdic.get) # Find maximum ratio in the subdomain
##########################################################################################
xi = np.linspace(x1.min(), x1.max(), N)
yi = np.linspace(y1.min(), y1.max(), N)
zi = scipy.interpolate.griddata((x1, y1), z1/z2, (xi[None,:], yi[:,None]), method='cubic')

import matplotlib.colors as colors
#smooth=np.linspace(0.0,1,5, endpoint=True)
smooth=[0,1,2,3,6]
cmap = colors.ListedColormap(['blue','thistle','yellow','red'])
bounds=smooth
norm = colors.BoundaryNorm(bounds, cmap.N)
#plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
#fig = plt.figure()
im=plt.contourf(xi, yi, abs(zi), smooth, cmap=cmap, norm=norm, rotation=90)

# uncomment the following to add space between tick labels and the axes
#plt.tick_params(axis='both', which='major', pad=15)
plt.xlabel("Along the fault (km)", fontsize=50, labelpad=25)
plt.ylabel("Perpendicular to the fault (km)", fontsize=50, labelpad=45)
plt.title('PID ratios for 40st_fn (Left/Center hypo)', y=1.15, fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.axis('scaled')
ax=plt.gca()
ax.tick_params(axis='both', which='major', pad=18)
#plt.colorbar()

from matplotlib import ticker
import matplotlib.colors as colors

#cmap = colors.ListedColormap(['g','y','r'])
bounds=smooth
norm = colors.BoundaryNorm(bounds, cmap.N)
cb = plt.colorbar(im, fraction=0.018, norm=norm, boundaries=bounds, ticks=smooth,)
#tick_locator = ticker.MaxNLocator(nbins=4)
#cb.locator = tick_locator
#cb.update_ticks()
cb.set_label('PID(L)/PID(C)', size=50, labelpad=40)
cb.ax.tick_params(size=30,labelsize=40, pad = 15)
#plt.savefig('140st_fnormal.pdf')
plt.savefig('40st-fn-pid-ratios.pdf')
#plt.show()
