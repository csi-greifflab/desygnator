## desygnator
__desYgnator__ is a computational method for comparing AIRR-seq samples. deYgnator stands for DEtection of SYstematic differences in GeneratioN of Adaptive immune recepTOr Repertoires


This repository corresponds to the following [research paper](https://www.biorxiv.org/content/10.1101/2021.04.19.440409v1):  
__Individualized VDJ recombination predisposes the available Ig sequence space.__ _bioRxiv_ (2021). Andrei Slabodkin, Maria Chernigovskaya, Ivana Mikocziova,  Rahmad Akbar, Lonneke Scheffer, Milena Pavlović, Igor Snapkov, Brij Bhushan Mehta, Habib Bashour, Jose Gutierrez-Marcos, Ludvig M. Sollid, Ingrid Hobaek Haff, Geir Kjetil Sandve, Philippe A. Robert, Victor Greiff

The "notebooks" directory contains the code (in the form of Jupyter notebooks) for the analysis described in the paper. 
To run it, you need to download the following [zip archive](https://doi.org/10.11582/2021.00089), unpack it and move the "models" subdirectory it into the root directory of the cloned repository. 
The unpacked folders will contain the [IGoR](https://github.com/qmarcou/IGoR) model parameter files corersponding to the fitted models analyzed in the paper, along with the fasta files these model parameters were inferred from. See, for example, `notebooks/greiff_cell_reports_pre.ipynb` for an example of parsing such files in Python.


