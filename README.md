# Gly-pipe
Software for the creation of new protein surface pockets.

Gly-pipe acts on protein structure files (.pdb) by substituting a heavy amino-acid with a glycine, a process known as glycinization.

Gly-pipe version 1.0

*******************************************************************************************************************

1) INSTALLATION

Gly-pipe is distributed as a python script and does not need particular installation procedures itself.
Simply extract the archive in the destination directory.

*******************************************************************************************************************

2) DEPENDENCIES

Gly-pipe runs with any python3 or python2 distribution. A python2 distribution is mandatory for SADIC software.

SADIC (Simple Atom Depth Index Calculator) is needed in order to calculate atom depth indexes. 
- It can be retrieved at: 
	http://sadic.sourceforge.net/
- To install SADIC, a python2 distribution is needed.
- Numarray package is needed in order to build and install SADIC. Numarray can be downloaded at:
	https://sourceforge.net/projects/numpy/files/Old%20Numarray/
- Numarray should be extracted and installed using the "setup.py" script included in the Numarray package.
- SADIC can then be installed using the "setup.py" script included in the SADIC package.

POPS (solvent accessible surface areas of proteins and nucleic acids) is needed in order to calculate the SASA of each residue.
- It can be downloaded at:
	https://github.com/Fraternalilab/POPS
- Documentation is available at:
	http://fraternalilab.github.io/POPS/

Pymol is needed, in order to mutate and minimize the structures.
- It can be found at:
	https://pymol.org/2/

Fpocket software is needed in order to estimate the pockets.
- It can be retrieved at:
	https://github.com/Discngine/fpocket

The python scikit-learn package is needed for the DS prediction.
- Described at:
	https://scikit-learn.org/
- Installation instructions at:
	https://scikit-learn.org/stable/install.html

*******************************************************************************************************************

3) USAGE

Gly-pipe can be run as any python script by calling: 
	python GlyPipe.py [args]

To analyze a protein structure, its ".pdb" file should be provided. 
- The file must be saved in the "structures" subdirectory.
- Then, to analyze it with GlyPipe just call:
	python GlyPipe.py [name] 
- [name] is the filename of the structure to be analyzed, left out the ".pdb" extension

For instance, if you want to analyze the HEWL (Henn Egg White Lysozime) structure:
- Search it on the PDB (https://www.rcsb.org/) through its pdb identifier "9lyz".
- Download the corresponding pdb structure file ("9lyz.pdb").
- Put the file in the "structures" subdirectory inside your Gly-pipe directory.
- Call:
	python GlyPipe.py 9lyz
- Using the default parameters in the GlyPipe.py script, you should obtain the following results:

Partial results are available for consultation in the other subdirectories.
- SADIC results, in which atom depth indexes can be found, are stored in "results_sadic".
- POPS results, where the SASA of each residue is stored, are saved in "results_pops".
- Mutant structures are saved in "mutated".
- Minimized vesrions of the mutant structures are saved in "minimized".
- Fpocket results are stored in "results_fpocket".

Energy minimization is needed to check the stability of each structure produced by GlyPipe.
- This is a resource-demanding step, and can be deactivated by setting the "do_minimize" parameter to False inside the script.
- Sometimes, results are meaningful also without minimizing the structure.
- GlyPipe can be used to build mutant structures which will be minimized later on a more powerful calculator.
- GlyPipe can bypass this step, by setting the parameter "already_minimized" to True.
	This works only if all the minimized structure versions exist in the "minimized" directory.
- GlyPipe can be stopped before minimization
	This is done by setting "mutate_only" to True.
	The mutant versions can be collected in "mutated" and minimized later and/or on a more powerful calculator.
	The minimized version of each mutant structure should then be put into "minimized".
	By setting "already_minimized" to True, GlyPipe will then estimate the pockets inside these structures and calculate their DS.

Fpocket estimation parameters can be changed by modifying the "fpocket_params" string in "GlyPipe.py".
- Some predefined versions are provided as commented lines inside the script.
- To choose your own pocket estimation parameters, check the fpocket manual.
- Changing the pocket estimation parameters significantly will decrease the accuracy of the DS predictor.

*******************************************************************************************************************
