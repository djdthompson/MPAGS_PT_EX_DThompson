# Particle Techniques Analysis Code
Folder containing python code and corresponding root files for analysis.
To run analysis user will require the following python packages:
* NumPy
* MatPlotLib
* itertools
* scipy.stats norm
* uproot

Most python installations through conda will have all of the above bar uproot as default.
To install this version of uproot (not uproot4 i believe) do:
```
conda install -c conda-forge uproot
```

# Documentation
The discussion of the results are on overleaf at https://www.overleaf.com/read/mrxnbgzfnhhc if you would prefer to clone the project:
```
git clone https://git.overleaf.com/600eed8e000e33094946b07a
```
I will occasionaly upload the latest pdf to this folder with a date to indicate version.

# Exercise 1
Open the analysis file, use lines 10-17 to choose the inital momentum and mag field setting, then run the analysis. It will print some events that have been skipped and produce a histogram of momentum with a gaussian fit. If you would like more details about certain events you can use the
```
plotFits(n,momentum,p=None)
```
function to plot event n with momentum an fitting values given. This will give you a plot of z-x and z-y to investigate the event.

# Exercise 2
Open the analysis file, use line 20 to choose the particle type, change the 0/1 in line 129 to choose whether to print a xy diagram of depositions in the calorimeter and change the settings according to the comments in lines 192 and onwards to get the correct fits. Then run the analysis.
