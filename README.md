# GammarusPulex_DEB-GUTS-T
This repository holds the code used to derive and analyze the scientific data presented in the ECOMOD publication: "How relevant are temperature corrections of toxicity parameters in population models for environmental risk assessment of chemicals?" (currently under review after resubmission)

The simulation raw data are stored in a Mendeley repository: Mangold-Döring, Annika; Buddendorf, Willem Bas; van den Brink, Paul; Baveco, Johannes M (2024), “GammarusPulex_DEB_GUTS_T_simulations”, Mendeley Data, V1, doi: 10.17632/8rymxsx46j.1

The folder SQGam.zip holds the DEB-IBM model written in the programming language Squeak 5.3 (www.squeak.org). To use the model, one needs a basic understanding of Squeak. First, unzip the folder and open the image file named “run_Sq6Gam.bat”. Further instructions exist to run the model in the image (see screenshot). 
 
![image](https://github.com/user-attachments/assets/253fb7b2-ad2c-4570-b21f-2ee318e458fb)

There are two R scripts to analyze the simulated models. DEB_GUTS_T_functions.R holds all the functions needed to read the data from the folder structure created by the Squeak program, and DEB_GUTS_T_analysis.R holds the script to execute the functions and create the figures presented in the publication. 
