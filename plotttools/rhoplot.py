import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np
import h5py
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator,LogLocator)
import matplotlib as mpl
import matplotlib.ticker
from matplotlib import cm


# Get all file names matching the pattern "element*.txt"
file_names = glob.glob("output/step_0/element_*.txt")

# Remove the ".txt" extension from each file name and sort them
# sorted_files = sorted(file_names, key=lambda x: int(x.lstrip("output/step_0/element_").rstrip(".txt")))

x = []
y = []
rho = []

# Load the data into a DataFrame
for dire in file_names:
    df = pd.read_csv(dire, delim_whitespace=True)
    x.append(df['x'])
    y.append(df['y'])
    rho.append(df['u1'])

x = np.array(x).flatten()
y = np.array(y).flatten()
rho = np.array(rho).flatten()

# Create the color plot
plt.scatter(x, y,  s=100, marker='s', c=rho, cmap='viridis')
plt.colorbar(label=r'$\rho$')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
# plt.xlim(-5,5)
# plt.ylim(-5,5)
plt.savefig('plotttools/rho.pdf',bbox_inches='tight')
plt.clf()