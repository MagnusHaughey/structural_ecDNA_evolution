

import numpy as np
import matplotlib.pyplot as plt 
import sys



# Take either single simulation file, or combine a list of files to look at averaged copy number distribution
singleCellData = []
for file in sys.argv[1:]:

	# Read file and add data to cumulative vector
	ecDNA_cn , multiplicity = np.loadtxt(file , unpack = True , ndmin = 2 , delimiter = "," , dtype = int)


	# Skip ecDNA-cold tumours 
	fileNameSplit = file.split("/")
	parameters = [word for word in fileNameSplit if ("Nmax" in word)][0]
	maxSize = parameters.split("_")[0]
	maxSize = int(maxSize.split("=")[-1])
	if (ecDNA_cn[0] == 0) and (multiplicity[0] == maxSize): continue


	# Unpack data into individual cell data to compute median, mean etc
	for pair in zip(ecDNA_cn , multiplicity):
		singleCellData += [pair[0] for _ in range(pair[1])]



# Plot copy number distribution
fig, ax = plt.subplots(figsize = (6 , 5))

maxCopyNumber = np.max(singleCellData)
bins = np.linspace(0 , maxCopyNumber , maxCopyNumber+1)
plt.hist(singleCellData , bins = bins , color = "#4528b8" , density = True)

plt.xlim([0 , maxCopyNumber + 10])
#plt.xticks([0,25,50])
#plt.yticks([0 , 0.05 , 0.1 , 0.15] , ["0.00" , "0.05" , "0.10" , "0.15"])

plt.ylabel("Counts (normalized)" , fontsize = 25)
plt.xlabel("ecDNA copies" , fontsize = 25)

plt.yscale("log")


### Formatting
for axis in ['bottom','left']:
	plt.gca().spines[axis].set_linewidth(1.5)

for axis in ['top','right']:
	plt.gca().spines[axis].set_linewidth(0)


plt.tick_params(labelsize = 20 , width = 1.5)

plt.tight_layout()
plt.subplots_adjust(right = 0.95)

if (len(sys.argv[1:]) == 1):
	outFile = sys.argv[1] + ".png"
else:
	outFile = "multiple_files_histogram.png"

plt.savefig(outFile , dpi = 500 , transparent = False  , format = 'png')








