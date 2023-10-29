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

with open("mom.txt", "r") as f:
    mom, err = f.read().split()
data = pd.read_csv('vel_prof.txt', delimiter = ' ', header = None, skiprows = 1)
data.columns = ["r", "vx"]


# Define function for fitting.
def func_lin(x, m, b):
    return m*x + b

# Define fitting routine (not necessary, but convenient).
def fit(x, y):
    popt, pcov = curve_fit(func_lin, x, y)
    perr = np.sqrt(np.diag(pcov))
    return popt[0], popt[1], perr[0], perr[1] # <- m, b, dm, db

vis = 0
m, b, dm, db = fit(data['r'][2:8], data['vx'][2:8])
print(r"dvx/dz: {:f}, error: {:f}".format(m, dm))
print(r"Flux $j_x$: {:f}".format(mom))
print("Viscosity: {:f}".format( abs(mom / m)))
vis = abs(mom / m)

# Define function for plot.
@mpltex.acs_decorator
def plot_shear_profile():
    fig, ax = plt.subplots(1,1)
    #fig, ax = plt.subplots(1,1, sharex=True, figsize=(3.25,2))
    linestyle=mpltex.linestyles(hollow_styles=[],markers=[],lines=['solid'])

    # Plot.

    ax.plot(data['r'][:11] - 1.5, data['vx'][:11], **next(linestyle))
    

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
