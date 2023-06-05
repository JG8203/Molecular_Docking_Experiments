# Research Plan

## Objective
The objective of this research is to investigate the antimicrobial and inhibitory potential of E. tirucalli phytochemicals by examining their antimicrobial activity on Xanthomonas oryzae Ddl enzyme through molecular docking.

## Specific Objectives
In the conduct of this study, the researchers aim to achieve the following specific objectives:

1. Determine the binding energies/affinities of the following phytochemicals (ligands) present in E. tirucalli and the control group through molecular docking simulations:
   - Arjunolic Acid
   - Eriodictyol
   - Afzelin
   - Scopoletin
   - 3,3′,4-trimethyl ellagic acid
   - Gallic acid
   - Piperic Acid
   - Terramycin (control)

2. Investigate the protein-ligand interactions of the seven phytochemicals towards Ddl enzyme in Xanthomonas oryzae.

3. Evaluate the drug-likeness and toxicity risk profiles of the seven phytochemicals with the best hits using SwissADME and Osiris Property Explorer.

## Hypothesis
Ho: There is no significant antimicrobial potential and inhibition of the phytochemicals of E. tirucalli against Xanthomonas oryzae Ddl enzyme as observed through molecular docking.

## Methodology

### Preparation of Ddl Enzyme
1. Download the 3D model of Ddl enzyme from the RCSB Protein Data Bank (https://www.rcsb.org/structure/1HXE) in the .pdb file format.
2. Prepare the compound by removing unnecessary characteristics such as water molecules and adding hydrogen atoms.
3. Convert the .pdb file of the target protein into the .pdbqt file format and save it in the desired file directory path.
4. Gather the coordinates of the active site where the docking will take place.

### Preparation of E. tirucalli Phytochemicals
1. Download the 3D models of the E. tirucalli phytochemicals from the National Library of Medicine (https://pubchem.ncbi.nlm.nih.gov).
2. Convert the compounds from their original file type (.sdf) to the .pdb file format using Avogadro software.
3. Load the .pdb files of the E. tirucalli phytochemicals into the Autodock Vina software.
4. Convert the compounds into the .pdbqt file format through the "Torsion Tree" settings.
5. Perform the same preparation procedures for the control compound (terramycin).

### Molecular Docking
1. Use the Autodock Vina software for simulating the E. tirucalli phytochemicals and control compound against Ddl enzyme.
2. Place the 3D models of the E. tirucalli phytochemicals, control compound, and Ddl enzyme (in .pdbqt format) in the Autodock Vina software.
3. Set up the grid parameters to define the exact location for the docking process.
4. Save the coordinates (x, y, z values) of the grid box in a separate text file.
5. Run the docking simulation by executing the command prompt for Windows and providing the directory and necessary dimensions.
6. After the simulation is completed, save the results in a text file for further analysis.

### Protein-Ligand Interaction Analysis
1. Identify the top-ranked configuration with the best binding affinity among the eight E. tirucalli phytochemicals from the simulation.
2. Analyze the docking results through visualization using PyMOL for 3D model visualization and Ligplot for 2D model visualization.
3

. Conduct Protein-Ligand Interaction Analysis to evaluate the interactions between the proteins and ligands.

### Drug Likeness and Toxicity Preparation
1. Illustrate the Drug Likeness and Toxicity Preparation levels of the niclosamide derivatives and the control compound.
2. Use the SwissADME website (http://www.swissadme.ch/) and Osiris Property Explorer for this analysis.

## Conclusion
The research plan outlines the objectives, methodology, and analysis approach for investigating the antimicrobial and inhibitory potential of E. tirucalli phytochemicals against Xanthomonas oryzae Ddl enzyme through molecular docking. The plan includes the preparation of Ddl enzyme and E. tirucalli phytochemicals, molecular docking simulations, protein-ligand interaction analysis, and evaluation of drug likeness and toxicity risk profiles. The hypothesis to be tested is also stated.
