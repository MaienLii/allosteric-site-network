# Allosteric Site Network
## Setup Environment
 We set up the environment using Anaconda.
 
 **Requirements**  
python == 3.7.12  
enspara == 0.1.1 （developed by bowman-lab）

## Input files need to be prepared
**trajectory_file**  *traj.nc*
**reference.pdb**  *protein.nc*
**pocket_residue.txt** *pocket.txt* (predicted by fpocket)

## 1 sidechain_deviation
Calculate the absolute deviation of each residue's side chain in each frame based on the trajectory file ***traj.nc***, and output the results to a file named ***sidechain_deviation.txt**.*


## 2 Cluster-MSMs-State
2.1 Perform **k-centers** and **k-medoids** clustering on the ***sidechain_deviation.txt*** file and select the appropriate **number of clusters** based on the implied timescales.

 2.2 Select the appropriate **lag time** based on the implied timescales to build the **Markov model**, and extract the **equilibrium distribution** of each microstate along with the side chain deviation values of its ***representative frames***.


## 3 Site-Network
Set **four thresholds** to divide the entire trajectory's side chain deviation values into five **equal-frequency intervals**. Input the files ***weight2.txt***, ***repre2.txt***, and the corresponding residue sequences of the predicted sites from ***pocket.txt***. Calculate the **allosteric network** between the sites, and rank the sites based on **eigenvector centrality**.

