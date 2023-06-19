# Load, import.
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# Format plots.
import mpltex
plt.rcParams['figure.figsize'] = [7, 4]
plt.rcParams.update({'font.size': 18})

# Read all data.
shear = ['3', '15', '60', '300', '1200']
data_list = []
mom_list = []
for i in shear:
    with open("mom_"+str(i)+".txt", "r") as f:
        mom, err = f.read().split()
        mom_list.append(float(mom))
    data = pd.read_csv('vel_prof_'+str(i)+'_av.txt', delimiter = ' ', header = None, skiprows = 1)
    data.columns = ["r", "vx"]
    data_list.append(data)
#data.columns = ['r', 'shear_3', 'shear_15', 'shear_60', 'shear_300', 'shear_1200'] 

# Define function for fitting.
def func_lin(x, m, b):
    return m*x + b

# Define fitting routine (not necessary, but convenient).
def fit(x, y):
    popt, pcov = curve_fit(func_lin, x, y)
    perr = np.sqrt(np.diag(pcov))
    return popt[0], popt[1], perr[0], perr[1] # <- m, b, dm, db

vis = []
n = [3, 15, 60, 300, 1200]
freq = [1 / i for i in n]

for i in range(5):
    m, b, dm, db = fit(data_list[i]['r'][2:8], data_list[i]['vx'][2:8])
    print(r"dvx/dz: {:f}, error: {:f}".format(m, dm))
    print(r"Flux $j_x$: {:f}".format(mom_list[i]))
    print("Viscosity: {:f}".format( abs(mom_list[i] / m)))
    vis.append(abs(mom_list[i] / m))

with open("viscosities.txt", "w") as file:
    for i in range(5):
        file.write("{:f} {:f}\n".format(freq[i], vis[i]))

# Define function for plot.
@mpltex.acs_decorator
def plot_shear_profile():
    fig, ax = plt.subplots(1,1)
    #fig, ax = plt.subplots(1,1, sharex=True, figsize=(3.25,2))
    linestyle=mpltex.linestyles(hollow_styles=[],markers=[],lines=['solid'])

    # Plot.

    ax.plot(data_list[0]['r'][:11] - 1.5, data_list[0]['vx'][:11], label='0.33', **next(linestyle))
    ax.plot(data_list[1]['r'][:11] - 1.5, data_list[1]['vx'][:11], label='0.07', **next(linestyle))
    ax.plot(data_list[2]['r'][:11] - 1.5, data_list[2]['vx'][:11], label='$1.67 \cdot 10^-2$', **next(linestyle))
    ax.plot(data_list[3]['r'][:11] - 1.5, data_list[3]['vx'][:11], label='$3.33 \cdot 10^-3$', **next(linestyle))
    ax.plot(data_list[4]['r'][:11] - 1.5, data_list[4]['vx'][:11], label='$8.33 \cdot 10^-4$', **next(linestyle))

    ax.axhline(y=0, color='black', linestyle='--')

    # Settings for plot.
    #ax.grid()
    ax.set_xlim([0, 15])
    ax.set_xlabel('$z / \sigma$')
    ax.set_ylabel('$v_x (m / \epsilon)^{(1/2)}$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.8, 1), ncol=1)
    fig.tight_layout()
    
    # Save file.
    fig.savefig('shear_profile.png')   

# Above, the plot was only defined. Execute.
plot_shear_profile()
