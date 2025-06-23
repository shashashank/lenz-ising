import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv, jv
import seaborn as sns
sns.set_theme()
plt.rcParams.update({
    "text.usetex": True,
    "font.family": 'computer-modern',
    "figure.dpi":300
})

reset = [0.01, 0.10, 0.14898, 10.00, 100.00]
Temp = [5.0]
for i in range(len(Temp)):
    T = Temp[i]
    r = reset[i]
    folder = '../bin/'
    resetData = np.loadtxt(folder+'resetData'+"{:.6f}".format(r)+'.txt')
    x = np.arange(0, 1, 0.01)
    plt.hist(resetData[:,0],bins=50, density=True, label='reset')
    plt.title(r'$P(E)$ vs $E$ for $r=$ '+str(r)+'; T='+str(T))
    plt.legend()
    plt.xlabel(r'$E$')
    plt.ylabel(r'$P(E)$')
    plt.show()
    plt.close()
