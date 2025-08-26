

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import pandas as pd 
import seaborn as sns




def sigmoid_function(x, a, b, s, ec_mean_size, ec_size_optimum, ec_size_range):
	return s*(1.0-(ec_mean_size-ec_size_optimum)/ec_size_range)/(1.0+np.exp(-a*(x-b)))




sigmoid_a = 0.5
sigmoid_b = 10
sigmoid_s = 1

ec_size_range = 10 - 1



x_toPlot = np.linspace(0 , 50 , 1000)
y_toPlot = np.linspace(0 , ec_size_range , 1000)
z_toPlot = np.zeros((len(x_toPlot) , len(y_toPlot)))
for x in range(len(x_toPlot)):
	for y in range(len(y_toPlot)):
		z_toPlot[x][y] = sigmoid_function(x_toPlot[x], sigmoid_a, sigmoid_b, sigmoid_s, y_toPlot[y], 1, ec_size_range)
		#print(x_toPlot[x],y_toPlot[y],z_toPlot[x][y])



# Plot 
plasma = mpl.colormaps['rocket_r']
newcolors = plasma(np.linspace(0 , 0.6, 256))
newcmp = ListedColormap(newcolors)

plt.imshow(z_toPlot.T , cmap = newcmp , vmin = 0 , vmax = sigmoid_s , interpolation = None , aspect = 'auto')

cbar = plt.colorbar(ticks = [0, 0.5, 1])
cbar.ax.set_yticklabels(['0' , "s/2" , 's'] , fontsize = 14)
cbar.set_label("Selective advantage", fontsize = 15 , labelpad = 17 , rotation = 270)


# Formatting
plt.gca().invert_yaxis()

custom_xticks = [0 , 200 , 400 , 600 , 800 , 1000]
custom_xticklabels = ["0" , "10" , "20" , "30" , "40" , "50"]
plt.gca().set_xticks(custom_xticks)
plt.gca().set_xticklabels(custom_xticklabels)


custom_yticks = [0 , 500 , 1000]
custom_yticklabels = ["0" , "5" , "10"]
plt.gca().set_yticks(custom_yticks)
plt.gca().set_yticklabels(custom_yticklabels)

plt.tick_params(labelsize = 12)


plt.xlabel("ecDNA copy number" , fontsize = 15 , labelpad = 10)
plt.ylabel("Distance from optimal ecDNA size" , fontsize = 15 , labelpad = 10)
#.set_xticklabels(["0" , "s"])


# axes["C"].set_xlabel("ecDNA copy number" , fontsize = 10)
# axes["D"].set_xlabel("ecDNA copy number" , fontsize = 10)
# axes["A"].set_ylabel("Cell replication rate" , fontsize = 10)
# axes["C"].set_ylabel("Cell replication rate" , fontsize = 10)

# axes["A"].set_yticks([0, sigmoid_s])
# axes["A"].set_yticklabels(
# 	["0" , "s"])

# axes["A"].legend()
plt.tight_layout()

plt.savefig("./sigmoid_fitnessFunction_withecDNAsizeDependence.png" , dpi = 600 , transparent = False , format = 'png')










