# pyright: reportUnusedImport=false, reportMissingImports=false, reportGeneralTypeIssues=false, reportOptionalIterable=false, reportOptionalMemberAccess=false, reportOptionalOperand=false, reportInvalidStringEscapeSequence=false
# pylint: disable-all
#coding=utf-8
import sys
import os
import shutil
import re
import math
import pymol
import numpy as np
import sadic
#call: python GlyPipe.py "PDB identifier"

#import classes
sys.path.insert(0, "classes/")
from atom_class import Atom
from residue_class import Residue
from pocket_class import Pocket

#constants
pdb_code = sys.argv[1]
sadic_version = sys.argv[2] # can be "v1" or "v2"
target_path = "structures/"+pdb_code+".pdb"
sasa_min = 5 # %
sasa_max = 20 # %
atom_to_residue = {}  #dictionary which tells the residue index (of res_list) corresponding to each atom id
dist_threshold = 5.00 #maximum distance (in angstroms) of the closest pocket from the glycinization target
pseudoatom_radius = 3.00 #radius for pseudoatoms created with pymol
max_targets = 1      #set to 0 to look for all the possible glycinization targets
best_depth_index = 0.7   #the glycinization targets closer to this value will have better ranking in the target selection
out_path = "results_glypipe/"
network_path = "network/network.pkl"
feature_scaler_path = "network/feature_scaler.pkl"
already_minimized = True       #if True, Gly-pipe will not minimize the mutant structures, searching for pre-minimized structures in the "minimized" subdirectory
do_minimize = False        #parameter to toggle structure minimization on/off (checked only if <already_minimized> is False) 

#list of aminoacids which can be substituted by Gly-pipe
glycinizable_aminoacids = ["ASN", "GLN", "ILE", "LEU", "MET", "PHE", "THR", "TRP", "TYR", "VAL"]

#fpocket parameters: uncomment the version to be used and comment the other ones
fpocket_params = "-i 15 -m 3.0 -M 4.0 -p 0.3"                                   #base params for Gly-pipe
#fpocket_params = "-i 25 -m 2.8 -M 4.0 -p 0.3 -r 5.0 -s 3.5 -n 3"               #mid-depth pockets
#fpocket_params = "-i 20 -m 2.5 -M 3.2 -p 0.3"                                  #deep pockets
#fpocket_params = "-i 25 -m 3.2 -M 5.5 -p 0.2"                                  #external pockets
#fpocket_params = "-i 30 -m 3.0 -M 5.0 -p 0.2 -r 5.2 -s 4.0 -n 4"               #large pockets

#lists of atoms whose d.i. must be compared to that of the residue CA, in order to determine the orientation of the lateral chain
sidechain_dict = {	"ASN" : ["CG", "OD1", "ND2"] ,
				"GLN" : ["CD", "OE1", "NE2"] ,
				"HIS" : ["ND1", "CD2", "CE1", "NE2"] ,
				"ILE" : ["CG1", "CG2", "CD1"] ,
				"LEU" : ["CG", "CD1", "CD2"] ,
				"MET" : ["CG", "SD", "CE"] ,
				"PHE" : ["CD1", "CD2", "CE1", "CE2", "CZ"] ,
				"SER" : ["CB", "OG"] ,
				"THR" : ["OG1", "CG2"] ,
				"TRP" : ["CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"] ,
				"TYR" : ["CD1", "CD2", "CE1", "CE2", "CZ", "OH"] ,
				"VAL" : ["CG1", "CG2"] }

#dictionary which associates an atom name to the corresponding value of Van Der Walls Radius
#taken from EPOS-BP software (PARSE.siz)
vdwr_dict =  {  "C": 1.7 , "CA": 1.7, "CB": 1.7, "CG": 1.7, "CG1": 1.7, "CG2": 1.7, "CD": 1.7, "CD1": 1.7, "CD2": 1.7, "CE": 1.7,
		"CE1": 1.7, "CE2": 1.7, "CE3": 1.7 , "CZ": 1.7, "CZ2": 1.7, "CZ3": 1.7, "CH2": 1.7, "N": 1.5, "NE": 1.5, "NE1": 1.5, 
		"NE2": 1.5, "NH1": 1.5, "NH2": 1.5, "ND1": 1.5, "ND2": 1.5, "NZ": 1.5, "O": 1.4, "OH": 1.4, "OG": 1.4, "OG1":1.4, 
		"OD1": 1.4, "OD2": 1.4, "OE1": 1.4, "OE2": 1.4, "OXT": 1.4, "SD": 1.85, "SG": 1.85}

#list of fpocket descriptor labels
descriptor_labels = [ 	"Score", "Number of Alpha Spheres", "Total SASA", "Polar SASA", "Apolar SASA", "Volume", 
			"Mean local hydrophobic density", "Mean alpha sphere radius", "Mean alp. sph. solvent access", 
			"Apolar alpha sphere proportion", "Hydrophobicity score", "Volume score", "Polarity score", "Charge score",
			"Proportion of polar atoms", "Alpha sphere density", "Cent. of mass - Alpha Sphere max dist", "Flexibility",
			"Druggability Score" ]

if sadic_version not in {"v1", "v2"}:
	raise Exception("Invalid SADIC version")

#execution

#check existence of the structure file for the protein in input
if not os.path.exists(target_path):
	print("Please provide the structure file for protein\""+pdb_code+"\"")

#run Gly-pipe
print("Running Gly-pipe on "+target_path)
#call sadic
print("Calling SADIC algorithm in order to calculate the Depth Index of each atom")
try:
	if sadic_version == "v1":
		os.system("python scripts/sadic.py "+target_path)
	elif sadic_version == "v2":
		sadic_result = sadic.sadic(target_path)
		sadic_result.save_txt(pdb_code+"_di.txt")

except:
	sys.exit("Exception caught during the execution of SADIC algorithm")
#move result file to SADIC result directory
os.rename(pdb_code+"_di.txt" , "results_sadic/"+pdb_code+"_di.txt")
#call pops for SASA calculation
print("Extracting the Solvent Accessible Surface Area of each residue")
try:
	os.system("pops --pdb structures/"+pdb_code+".pdb --atomOut")
except:
	sys.exit("Exception caught during surface analysis")
os.rename("pops.out","results_pops/"+pdb_code+"_sasa.txt")
try:
	os.remove("popsb.out")
	os.remove("sigma.out")
except:
	sys.exit("Exception while removing extra output files produced by pops")	
	
#analysis of results produced by sadic and pops
pops_file = open("results_pops/"+pdb_code+"_sasa.txt","r")
pops_text = pops_file.read()
pops_file.close()
pops_lines = pops_text.splitlines(pops_text.count("\n"))
del pops_text
res_list = list()
#define identifier
res_id = -1000
#first element (useful only for list initialization)
res_index = -1
#skip first three lines (header)
for i in range(3,len(pops_lines)):
	cells = re.split("[\s]+",pops_lines[i])
	#check if the line corrsponds to an entry
	if len(cells) < 14:
		continue
	#ignore heteroatoms and DNA
	if cells[3] in ["HET", "DA", "DG", "DC", "DT"]:
		continue
	#check whether the line corresponds to an atom belonging to the current residue, or a new residue must be built
	if int(res_id)!=int(cells[5]):
		atom_id = int(cells[1])
		res_id = cells[5]
		res_index = res_index+1
		res_list.append(Residue(cells[3], res_id, cells[4]))
		res_list[res_index].insert(Atom(cells[7],cells[2],cells[1]))	
		atom_to_residue[atom_id]=res_index
	else:
		atom_id = int(cells[1])
		res_list[res_index].insert(Atom(cells[7],cells[2],cells[1]))
		atom_to_residue[atom_id]=res_index
#sadic algorithm results
sadic_file = open("results_sadic/"+pdb_code+"_di.txt","r")
sadic_text = sadic_file.read()
sadic_file.close()
sadic_lines = sadic_text.splitlines(sadic_text.count("\n"))
for i in range(1,len(sadic_lines)):
	#each line corresponds to an atom
	cells = re.split("[\s]+",sadic_lines[i])
	cell_index = len(cells) - 2 if sadic_version == "v1" else 1
	depth_index = cells[cell_index]
	atom_id = int(cells[0])
	#skip secondary conformers (not reported by pops, thus not inserted in <atom_to_residue>)
	if atom_id not in atom_to_residue.keys():
		continue
	#set DI
	res_index = atom_to_residue[atom_id]
	res_list[res_index].setDIForAtom(atom_id, depth_index)
#scan PDB file in order to retrieve residue positions
pdb_file = open("structures/"+pdb_code+".pdb","r")
pdb_text = pdb_file.read()
pdb_file.close()
pdb_lines = pdb_text.splitlines(pdb_text.count("\n"))
for l in pdb_lines:
	#only "ATOM" lines, avoiding secondary conformers and nucleic acids
	if l[0:4]=="ATOM" and l[16] in [" ", "A"] and l[17:20] not in [" DA", " DG", " DC", " DT"]:
		#skip hydrogens
		if l[12]=="H" or l[13]=="H":
			continue
		#collect atom data
		atom_id = int(l[6:11])
		px = float(l[31:38])
		py = float(l[39:46])
		pz = float(l[47:54])
		#find the residue this atom belongs to
		res_index = atom_to_residue[atom_id]
		#update residue instance with position data for this atom
		res_list[res_index].setPositionFor(atom_id, px, py, pz)

#use acquired atom positions to calculate residue positions
for r in res_list:
	r.calculatePosition()

#select from <res_list>: Phe, Tyr and Trp residues, having SASA greater or equal than the minimum value, and the side-chain oriented upwards
print("Selecting potential targets for glycinization")
gt_list = []
for r in res_list:
	#check exposition to the solvent
	if r.getName() in glycinizable_aminoacids and r.checkExposition(sasa_min, sasa_max):
		#check side-chain orientation
		upwards = True
		alpha_carbon_di = r.getAtom("CA").getDepthIndex()
		for atom_name in sidechain_dict.get(r.getName()):
			if r.getAtom(atom_name) is None:
				upwards = False #uncertain configuration of the sidechain is assumed to be not good
				continue
			if r.getAtom(atom_name).getDepthIndex() < alpha_carbon_di:
				upwards = False
		#select the residue only if its side-chain is oriented upwards
		if upwards:
			gt_list.append(r)

#rank glycinization targets by depth index
gt_scores = []
for r in gt_list:
	gt_scores.append(abs(r.getDepthIndex() - best_depth_index))
#get sorted indices list (best scores are the smallest, so they will be at the top of the list)
sorted_indices = list(np.argsort(gt_scores))
#if <max_targets> is set to 0 (or smaller), set it to the total number of targets
if max_targets<=0:
	max_targets =len(sorted_indices)
#use the sorted indices to sort the glycinization target list and the scores list
fgt_list = []
fgt_scores = []
i = 0
#stop if all the targets were added to the list or the maximum number of targets has been reached
while (i < len(sorted_indices)) and (i < max_targets):
	fgt_list.append( gt_list[ sorted_indices[i] ] )
	fgt_scores.append( gt_scores[ sorted_indices[i] ] )
	i += 1
		
#print glycinization targets
print("Glycinization targets, by distance to ideal DI:")
for i in range(len(fgt_list)):
	print(fgt_list[i].toString() +"   "+ str(fgt_scores[i]))

#print the results to file
out_file = open(out_path+pdb_code+".txt", 'w')
out_file.write("Glycinization prediction of "+pdb_code+"\n")
out_file.write("Minimum sasa = "+str(sasa_min)+"\n")
#for each glycinization target
for fgt in fgt_list:
	#write target entry
	out_file.write("\nTarget = "+fgt.toString()+"\n")
out_file.flush()
out_file.close()

#create directory to store mutated versions of the structure (if it doesn't exist)
if not os.path.exists("mutated/"+pdb_code+"/"):
	os.makedirs("mutated/"+pdb_code+"/")

#open pymol
pymol.pymol_argv = ['pymol','-qc']
pymol.finish_launching()

#create mutated structures
for fgt in fgt_list:
	print("Substituting "+fgt.toString()+" with GLY")
	#reinitialize pymol environment
	pymol.cmd.reinitialize()
	#load natural structure
	pymol.cmd.load(target_path, pdb_code)
	#remove all the heteroatoms
	pymol.cmd.remove('hetatm')
	#launch mutagenesis wizard
	pymol.cmd.wizard("mutagenesis")
	pymol.cmd.do("refresh_wizard")
	#perform substitution
	pymol.cmd.get_wizard().set_mode("GLY")
	pymol.cmd.get_wizard().do_select(str(fgt.getNumber())+"/")
	pymol.cmd.get_wizard().apply()
	pymol.cmd.set_wizard()
	#add a pseudoatom representing the pocket, centered on the position of <fgt> and with radius equal to <pseudoatom_radius>
	pymol.creating.pseudoatom(pdb_code,'','AP1','AP','1','P','','X',pseudoatom_radius,1,0.0,0.0,'tv_yellow','',fgt.getPosition(),1)
	#save mutated structure
	pymol.cmd.save("mutated/"+pdb_code+"/"+pdb_code+"_"+fgt.getFilenameString()+".pdb", "all", -1, "pdb")

#quit pymol
pymol.cmd.quit()

#create protein subfolder for fpocket results
if not os.path.exists("results_fpocket/"+pdb_code+"/"):
	os.makedirs("results_fpocket/"+pdb_code+"/")
#copy wild-type (natural) structure
print("Creating a copy of the wild-type structure for fpocket")
shutil.copyfile(target_path,"results_fpocket/"+pdb_code+"/"+pdb_code+".pdb")

#manage energy minimization procedure
if already_minimized:
	#copy minimized structures into fpocket results folder
	print("Creating copies of the minimized structure versions for fpocket")
	for fgt in fgt_list:
		#build filename "structure id" string for <fgt>
		fgt_filename = pdb_code+"_"+fgt.getFilenameString()+".pdb"
		#check whether the pre-minimized structure is actually available in the "minimized" subdirectory
		if not os.path.exists("minimized/"+pdb_code+"/"+fgt_filename):
			sys.exit("Please provide minimized structure for mutant version "+fgt_filename+" in the \"minimized\" subdirectory")
		#copy minimized mutant structure
		shutil.copyfile("minimized/"+pdb_code+"/"+fgt_filename,"results_fpocket/"+pdb_code+"/"+fgt_filename)
else:
	if do_minimize:
		#create minimized structure versions
		print("Calculating minimized versions of the protein structure")
		for fgt in fgt_list:
			#build filename "structure id" string for <fgt>
			fgt_filename = pdb_code+"_"+fgt.getFilenameString()+".pdb"
			try:
				os.system("python scripts/minimize.py mutated/"+pdb_code+"/"+fgt_filename+" minimized/"+pdb_code+"/"+fgt_filename)
			except:
				sys.exit("Exception caught during energy minimization")
		#copy minimized structures into fpocket results folder
		print("Creating copies of the minimized structure versions for fpocket")
		for fgt in fgt_list:
			#build filename "structure id" string for <fgt>
			fgt_filename = pdb_code+"_"+fgt.getFilenameString()+".pdb"
			#copy minimized structure
			shutil.copyfile("minimized/"+pdb_code+"/"+fgt_filename,"results_fpocket/"+pdb_code+"/"+fgt_filename)		
	else:
		#copy mutant structures into fpocket results folder, without minimizing them
		print("Creating structure copies for fpocket")
		shutil.copyfile(target_path,"results_fpocket/"+pdb_code+"/"+pdb_code+".pdb")
		for fgt in fgt_list:
			#build filename "structure id" string for <fgt>
			fgt_filename = pdb_code+"_"+fgt.getFilenameString()+".pdb"
			#copy mutant structure
			shutil.copyfile("mutated/"+pdb_code+"/"+fgt_filename,"results_fpocket/"+pdb_code+"/"+fgt_filename)

#run fpocket on each structure
os.system("fpocket "+fpocket_params+" -f results_fpocket/"+pdb_code+"/"+pdb_code+".pdb")
for fgt in fgt_list:
	#run on <fgt> mutated structure
	os.system("fpocket "+fpocket_params+" -f results_fpocket/"+pdb_code+"/"+pdb_code+"_"+fgt.getFilenameString()+".pdb")

#collect fpocket results
for fgt in fgt_list:
	#acquire pockets for each mutant structure
	b_name = pdb_code+"_"+fgt.getFilenameString()
	b_file = open("results_fpocket/"+pdb_code+"/"+b_name+"_out/"+b_name+"_info.txt","r")
	b_file_text = b_file.read()
	b_file.close()
	b_file_lines = b_file_text.splitlines(b_file_text.count("\n"))
	del b_file_text
	pocket_list = None
	current_pocket = None
	for l in b_file_lines:
		#empty line -> finalize pocket data
		if len(l)==1:
			#in case <current_pocket> is not initialized, just skip blank line
			if current_pocket is not None:
				#if this is the first pocket, <pocket_list> must be initialized
				if pocket_list is None:
					pocket_list = [current_pocket]
					current_pocket = None
				else:
					pocket_list.append(current_pocket)
					current_pocket = None					
			continue
		#pocket declaration line
		if re.match("Pocket", l):
			c = re.split("[ ]", l)
			current_pocket = Pocket(c[1])
			continue
		#otherwise, this is a descriptor line
		c = re.split("[\s]+", l)
		key = c[1]
		i = 2
		#add line parts to <key> until the symbol ":" is reached, in order to build the descriptor string
		while c[i] not in [ ":", "score:", "atoms:", "dist:" ]:
			key = key + " " +c[i]
			i = i + 1
		#if ":" symbol is attached to the final word of the dscriptor, it must be removed, and the word appended to <key>
		if c[i] in [ "score:", "atoms:", "dist:" ]:
			key = key + " " + c[i][:-1]
		#value always follows the ":" symbol
		value = float(c[i+1])
		#add <key>, <value> pair to the pocket's descriptor dictionary
		current_pocket.insert(key, value)
	del b_file_lines
	#for each pocket in the list, scan the corresponding file in order to acquire the position of its vertices
	for p in pocket_list:
		#open the file where pocket vertices are declared
		vert_file = open("results_fpocket/"+pdb_code+"/"+b_name+"_out/pockets/pocket"+str(p.getID()-1)+"_vert.pqr","r")
		vert_text = vert_file.read()
		vert_file.close()
		vert_lines = vert_text.splitlines(vert_text.count("\n"))
		del vert_text
		#scan lines and acquire the coordinates of each pocket vertex
		for vl in vert_lines:
			#only "ATOM" lines must be considered
			if vl[0:4] == "ATOM":
				p.addVertex( float(vl[31:38]), float(vl[39:46]), float(vl[47:54]))
	#search for the pocket which is closest to <fgt>
	closest_pocket = None
	min_dist = None
	for p in pocket_list:
		if closest_pocket is None:
			closest_pocket = p
			min_dist = p.getMinDistance(fgt.getPosition())
			continue
		dist = p.getMinDistance(fgt.getPosition())
		if dist < min_dist :
			closest_pocket = p 
			min_dist = dist
	#if the pocket is located too far away, the mutation was not effective
	if min_dist > dist_threshold :
		print("No pocket was detected for "+fgt.toString())
		print("DEBUG: Distance = "+str(min_dist))
		#append glycinization results to the output file
		out_file = open(out_path+pdb_code+".txt", 'a')
		out_file.write("\n\n\nThe glycinization of "+fgt.toString()+" did not create a new pocket\n")
		out_file.flush()
		out_file.close()
	#else, give Fpocket's DS as an approximate evaluation of the pocket and suggest the structure for PockDrug analysis
	else:
		fpocket_ds = closest_pocket.getValue("Druggability Score")
		#print results
		print("Pocket |"+str(closest_pocket.getID())+"| detected for "+fgt.toString()+" with Fpocket's DS = "+str(fpocket_ds))
		print("DEBUG: Distance = "+str(min_dist))
		#append glycinization results to the output file
		out_file = open(out_path+pdb_code+".txt", 'a')
		out_file.write("\n\n\nPocket created by the glycinization of "+fgt.toString()+" :\n")
		out_file.write("Pocket ID = "+str(closest_pocket.getID())+"\n")
		out_file.write("Fpocket's Druggability Score = "+str(fpocket_ds)+"\n")
		out_file.write("Volume = "+str(closest_pocket.getValue("Volume"))+"\n")
		if already_minimized or do_minimize:
			out_file.write("Path to structure for PockDrug analysis : \"minimized/"+pdb_code+"/"+pdb_code+"_"+fgt.getFilenameString()+".pdb\"")
		else:
			out_file.write("Path to structure for PockDrug analysis : \"mutated/"+pdb_code+"/"+pdb_code+"_"+fgt.getFilenameString()+".pdb\"")
		out_file.flush()
		out_file.close()
		
