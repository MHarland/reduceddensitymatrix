import numpy as np
from mpltotex import PRLPlotConf

mplcfg = PRLPlotConf()
ax = mplcfg.ax
data =np.load('mutual_information.npy')
ax.plot(data[:,3], data[:,5], label = '$\mathrm{nn}$')
ax.plot(data[:,3], data[:,6], label = '$\mathrm{nnn}$')
ax.set_xlim(data[:,3].min(), data[:,3].max())
ax.set_xlabel('$U$')
ax.set_ylabel('$\mathcal{I}$')
ax.legend(**mplcfg.legendkwargs)
mplcfg.save('mutual_information')
