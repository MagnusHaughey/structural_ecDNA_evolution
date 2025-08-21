


import numpy as np 
import matplotlib.pyplot as plt 
import glob
import seaborn as sns 
import pandas as pd 





def sigmoid_function(x, a, b, s):
	return s/(1.0+np.exp(-a*(x-b)))



s_color_map = {
	0: "#173F3F",
	0.1: "#20639B",
	0.5: "#3CAEA3",
	1: "#F6D55C",
	2: "#ED553B"
}



all_s = [0 , 0.1 , 0.5 , 1 , 2]
all_a = [0.5 , 1 , 2]
all_b = [5 , 10 , 20]
maxSize = 100000


dataToPlot = []
for a in all_a:
	for  b in all_b:

		plt.clf()
		fig , ax = plt.subplot_mosaic("AB" , figsize = (8 , 5) , width_ratios = [2,3])

		for s in all_s:


			# Loop over all files for this particular selection strength
			for file in glob.glob("/Users/ahw621/Documents/Projects/ecDNA/structural_ecDNA_evolution/2025_08_04/RESULTS/Nmax={}_initialCopyNumber=1_s={}_sigmoidA={}_sigmoidB={}/*/tissue.csv".format(maxSize , s , a , b)):

				# Read file and add data to cumulative vector
				ecDNA_cn , multiplicity = np.loadtxt(file , unpack = True , ndmin = 2 , delimiter = "," , dtype = int)


				# Skip ecDNA-cold tumours 
				if (ecDNA_cn[0] == 0) and (multiplicity[0] == maxSize): continue


				# Unpack data into individual cell data to compute median, mean etc
				singleCellData = []
				for pair in zip(ecDNA_cn , multiplicity):
					singleCellData += [pair[0] for _ in range(pair[1])]


				dataToPlot.append([s , 1/(1+s) , np.mean(singleCellData)])


				# Now plot selection function for given s value
				x_toPlot = np.linspace(0 , 50 , 1000)
				y_toPlot = [sigmoid_function(val, a, b, s) for val in x_toPlot]

				ax["B"].plot(x_toPlot , y_toPlot , lw = 3 , color = s_color_map[s])




		# Turn into dataframe and plot 
		dataToPlot_df = pd.DataFrame(dataToPlot , columns = ["s" , "1/(1+s)" , "Median ecDNA copy number"])


		boxPlot = sns.boxplot(data = dataToPlot_df , 
			x = "1/(1+s)" ,
			hue = "1/(1+s)" ,	
			y = "Median ecDNA copy number" , 
			width = 0.8 ,
			gap = 0.4 ,
			linewidth = 1.5 , 
			legend = False ,
			#color = "#eb4034" ,
			palette = list(s_color_map.values())[::-1],
			fill = True ,
			showfliers = False ,
			ax = ax["A"])


		sns.stripplot(data = dataToPlot_df , 
			x = "1/(1+s)" ,
			hue = "1/(1+s)" ,	
			y = "Median ecDNA copy number" ,
			legend = False ,
			#dodge = True ,
			#color = "#eb4034" ,
			palette = list(s_color_map.values())[::-1],
		 	edgecolor = "black" ,
			linewidth = 1 , 
			size = 6 ,
			zorder = 10 ,
			ax = ax["A"] )







		# Formatting
		ax["A"].tick_params(labelsize = 15 , bottom = False , width = 1.5)
		ax["A"].set_xlabel("1/(1+s)" , fontsize = 15 , labelpad = 10)
		ax["A"].set_ylabel("Median ecDNA copy number" , fontsize = 15 , labelpad = 10)

		ax["A"].set_xticks([i for i in range(len(boxPlot.get_xticklabels()))] , [round(float(t.get_text()) , 2) for t in boxPlot.get_xticklabels()] , rotation = 45 , ha = 'right')

		#boxPlot.set_xticklabels([float(t.get_text()) for t in boxPlot.get_xticklabels()])




		# Now plot selection function
		x_toPlot = np.linspace(0 , 50 , 1000)
		y_toPlot = [sigmoid_function(val, a, b, s) for val in x_toPlot]

		ax["B"].plot(x_toPlot , y_toPlot , lw = 3 , color = "#ba0000" , label = "a = {}".format(a))



		ax["B"].tick_params(labelsize = 15 , bottom = False , width = 1.5)
		ax["B"].set_xlabel("ecDNA copy number" , fontsize = 15 , labelpad = 30)
		ax["B"].set_ylabel("Selective advantage" , fontsize = 15)

		ax["B"].set_yticks(all_s)



		# More formatting
		ax["A"].spines[['right', 'top']].set_visible(False)
		ax["B"].spines[['right', 'top']].set_visible(False)

		for axis in ['bottom','left']:
			ax["A"].spines[axis].set_linewidth(1.5)
			ax["B"].spines[axis].set_linewidth(1.5)


		plt.tight_layout()

		plt.subplots_adjust(wspace = 0.5)
		plt.savefig("/Users/ahw621/Documents/Projects/ecDNA/structural_ecDNA_evolution/2025_08_04/PLOTS/Nmax={}_initialCopyNumber=1_sigmoidA={}_sigmoidB={}_meanCopyNumber_vs_selection.png".format(maxSize , a , b) , dpi = 600 , transparent = False , format = 'png')


		print("Plotted: /Users/ahw621/Documents/Projects/ecDNA/structural_ecDNA_evolution/2025_08_04/PLOTS/Nmax={}_initialCopyNumber=1_sigmoidA={}_sigmoidB={}_meanCopyNumber_vs_selection.png".format(maxSize , a , b))




