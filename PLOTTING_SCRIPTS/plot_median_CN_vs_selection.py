


import numpy as np 
import matplotlib.pyplot as plt 
import glob
import seaborn as sns 
import pandas as pd 


all_s = [0 , 0.1 , 0.2 , 0.5 , 1 , 2]
maxSize = 1000000


dataToPlot = []
for s in all_s:

	# Loop over all files for this particular selection strength
	for file in glob.glob("./RESULTS/Nmax={}*_s={}/*/tissue.csv".format(maxSize , s)):

		# Read file and add data to cumulative vector
		ecDNA_cn , multiplicity = np.loadtxt(file , unpack = True , ndmin = 2 , delimiter = "," , dtype = int)


		# Skip ecDNA-cold tumours 
		if (ecDNA_cn[0] == 0) and (multiplicity[0] == maxSize): continue


		# Unpack data into individual cell data to compute median, mean etc
		singleCellData = []
		for pair in zip(ecDNA_cn , multiplicity):
			singleCellData += [pair[0] for _ in range(pair[1])]


		dataToPlot.append([s , 1/(1+s) , np.mean(singleCellData)])




# Turn into dataframe and plot 
dataToPlot_df = pd.DataFrame(dataToPlot , columns = ["s" , "1/(1+s)" , "Median ecDNA copy number"])


fig , ax = plt.subplots(figsize = (4 , 6))
boxPlot = sns.boxplot(data = dataToPlot_df , 
	x = "1/(1+s)" ,	
	y = "Median ecDNA copy number" , 
	width = 0.8 ,
	gap = 0.4 ,
	linewidth = 1.5 , 
	legend = True ,
	color = "#eb4034" ,
	fill = True ,
	showfliers = False )


sns.stripplot(data = dataToPlot_df , 
	x = "1/(1+s)" ,	
	y = "Median ecDNA copy number" ,
	legend = False ,
	dodge = True ,
	color = "#eb4034" ,
 	edgecolor = "black" ,
	linewidth = 1 , 
	size = 6 ,
	zorder = 10 )







# Formatting
plt.tick_params(labelsize = 15 , bottom = False , width = 1.5)
plt.xlabel("1/(1+s)" , fontsize = 15 , labelpad = 10)
plt.ylabel("Median ecDNA copy number" , fontsize = 15 , labelpad = 10)

plt.xticks([i for i in range(len(boxPlot.get_xticklabels()))] , [round(float(t.get_text()) , 2) for t in boxPlot.get_xticklabels()] , rotation = 45 , ha = 'right')

#boxPlot.set_xticklabels([float(t.get_text()) for t in boxPlot.get_xticklabels()])


plt.gca().spines[['right', 'top']].set_visible(False)

for axis in ['bottom','left']:
	plt.gca().spines[axis].set_linewidth(1.5)



plt.tight_layout()
plt.savefig("./PLOTS/maxSize={}_meanCopyNumber_vs_selection.png".format(maxSize) , dpi = 600 , transparent = False , format = 'png')








