import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np
import h5py
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator,LogLocator)
import matplotlib as mpl
import matplotlib.ticker
from matplotlib import cm
import matplotlib.tri as tri
import matplotlib.animation as animation

mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2





































# plot individual elements 

# dir_names = glob.glob("output/step_*")
# dir_names = sorted(dir_names, key=lambda x: int(x.lstrip("output/step_")))

# for i in np.arange(len(dir_names))[::200]:

#     file_names = glob.glob(dir_names[i] + "/element_*.txt")

#     x = []
#     y = []
#     u1 = []
#     u2 = []
#     u3 = []
#     u4 = []

#     for dire in file_names:

#         df = pd.read_csv(dire, sep='\s+')
#         x.append(df['x'])
#         y.append(df['y'])
#         u1.append(df['u0'])
#         u2.append(df['u1'])
#         u3.append(df['u2'])
#         u4.append(df['u3'])

#     gamma = 1.4

#     x = np.array(x).flatten()
#     y = np.array(y).flatten()
#     rho = np.array(u1).flatten()
#     u = np.array(u2).flatten() / rho
#     v = np.array(u3).flatten() / rho
#     E = np.array(u4).flatten() / rho
#     p   = rho * ( gamma - 1 ) * ( E - ( u**2 + v**2 ) / 2 )
#     H   = E + p / rho

#     triangles = tri.Triangulation(x, y)
#     plt.tripcolor(triangles, rho, cmap='viridis')
#     plt.colorbar(label=r'$\rho$')  
#     plt.xlabel(r'$x$')
#     plt.ylabel(r'y')
#     plt.savefig('plotttools/step_'+ str(i) +'.pdf',bbox_inches='tight')
#     plt.clf()




































import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import glob

# Get directory names
dir_names = glob.glob("output/step_*")
dir_names = sorted(dir_names, key=lambda x: int(x.lstrip("output/step_")))

# Loop through each frame
for i in range(len(dir_names))[::20]:

    print(dir_names[i])

    # Load data for the current frame
    file_names = glob.glob(dir_names[i] + "/element_*.txt")
    x, y, u1, u2, u3, u4, t = [], [], [], [], [], [], []
    for dire in file_names:
        df = pd.read_csv(dire, sep='\s+')
        x.append(df['x'])
        y.append(df['y'])
        u1.append(df['u0'])
        u2.append(df['u1'])
        u3.append(df['u2'])
        u4.append(df['u3'])
        t.append(list(df['time'])[0])

    # Flatten the data
    x = np.array(x).flatten()
    y = np.array(y).flatten()
    rho = np.array(u1).flatten()
    u = np.array(u2).flatten() / rho
    v = np.array(u3).flatten() / rho
    E = np.array(u4).flatten() / rho
    gamma = 1.4
    p = rho * (gamma - 1) * (E - (u**2 + v**2) / 2)
    H = E + p / rho

    # Create the triangular mesh
    triangles = tri.Triangulation(x, y)

    # Create subplots with shared axes
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(10, 8))

    # Plot each variable in a separate subplot
    im1 = axs[0, 0].tripcolor(triangles, rho, cmap='viridis')
    axs[0, 0].set_title(r'$\rho$')
    fig.colorbar(im1, ax=axs[0, 0])

    im2 = axs[0, 1].tripcolor(triangles, u, cmap='viridis')
    axs[0, 1].set_title(r'$u$')
    fig.colorbar(im2, ax=axs[0, 1])

    im3 = axs[1, 0].tripcolor(triangles, p, cmap='viridis')
    axs[1, 0].set_title(r'$p$')
    fig.colorbar(im3, ax=axs[1, 0])

    im4 = axs[1, 1].tripcolor(triangles, v, cmap='viridis')
    axs[1, 1].set_title(r'$v$')
    fig.colorbar(im4, ax=axs[1, 1])

    plt.suptitle(f'Time: {t[0]}')

    # Set common labels
    for ax in axs.flat:
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    plt.savefig('plotttools/step_' + dir_names[i].lstrip("output/step_") + '.pdf', bbox_inches='tight')

    # Close the plot to avoid memory leaks
    plt.close(fig)



























# Function to generate data for each frame of the animation
def update(frame):
    print(frame)
    # Load data for the current frame
    file_names = glob.glob(dir_names[frame] + "/element_*.txt")
    x, y, rho = [], [], []
    for dire in file_names:
        df = pd.read_csv(dire, sep='\s+')
        x.append(df['x'])
        y.append(df['y'])
        rho.append(df['u0'])

    # Flatten the data
    x = np.array(x).flatten()
    y = np.array(y).flatten()
    rho = np.array(rho).flatten()

    # Create the triangular mesh
    triangles = tri.Triangulation(x, y)

    # Clear the plot
    plt.clf()

    # Plot the triangular mesh with color based on sample data
    plt.tripcolor(triangles, rho, cmap='viridis')
    plt.colorbar(label=r'$\rho$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')

# Get directory names
dir_names = glob.glob("output/step_*")
dir_names = sorted(dir_names, key=lambda x: int(x.lstrip("output/step_")))

# Create a figure
fig = plt.figure()

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(dir_names), interval=100)

# Save the animation as a video file
ani.save('plotttools/animation.mp4', writer='ffmpeg')



