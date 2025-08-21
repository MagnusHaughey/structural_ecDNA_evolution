

import numpy as np 
import matplotlib.pyplot as plt 



def sigmoid_function(x, a, b, s):
	return s/(1.0+np.exp(-a*(x-b)))




sigmoid_a = [0.2, 0.5, 1, 2]
sigmoid_b = [5, 10, 15, 20]
sigmoid_s = 0.5


axes_map = {
	5: "A",
	10: "B",
	15: "C",
	20: "D"
}


mosaic = '''AB
	CD'''
fig, axes = plt.subplot_mosaic(mosaic = mosaic , figsize = (7,5) , sharex = True , sharey = True)

for b in sigmoid_b:
	axes[axes_map[b]].set_title("b = {}".format(b))

	for a in sigmoid_a:
		x_toPlot = np.linspace(0 , 50 , 1000)
		y_toPlot = [sigmoid_function(val, a, b, sigmoid_s) for val in x_toPlot]

		axes[axes_map[b]].plot(x_toPlot , y_toPlot , lw = 2 , label = "a = {}".format(a))



axes["C"].set_xlabel("ecDNA copy number" , fontsize = 10)
axes["D"].set_xlabel("ecDNA copy number" , fontsize = 10)
axes["A"].set_ylabel("Cell replication rate" , fontsize = 10)
axes["C"].set_ylabel("Cell replication rate" , fontsize = 10)

axes["A"].set_yticks([0, sigmoid_s])
axes["A"].set_yticklabels(
	["0" , "s"])

axes["A"].legend()
plt.tight_layout()

plt.savefig("./sigmoid_fitnessFunction.png" , dpi = 600 , transparent = False , format = 'png')