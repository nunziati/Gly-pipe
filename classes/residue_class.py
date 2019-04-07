#coding=utf-8
import sys
import os
import math

#import atom class
from atom_class import Atom

#dictionary of nominal ASA values for each aminoacid
asa_dict = { 'GLY': 1.0 , 'ALA': 67.0 , 'ARG':196.0 , 'ASN':113.0 , 'ASP':106.0 , 'CYS':104.0 , 'GLN':144.0, 'GLU':138.0 , 'HIS':151.0 , 'ILE':140.0 , 'LEU':137.0 , 'LYS':167.0 , 'MET':160.0 , 'PHE':175.0 , 'PRO':105.0 , 'SER':80.0 , 'THR':102.0 , 'TRP':217.0 , 'TYR':187.0 , 'VAL':117.0  }

#residue object class
class Residue:
	
	#constructor
	def __init__(self, res_name, res_number, res_chain):
		self.atoms = []            #list of atoms
		self.name = res_name       #residue code (f.i.: glycine -> gly)
		self.number = res_number   #position in chain sequence
		self.chain = res_chain     #chain identifier
		self.position = None       #position in 3D space, defined in the PDB file	

	#inserts an atom belonging to this residue
	def insert(self, atom):
		self.atoms.append(atom)
	
	#sets the position of the atom corresponding to <number> (which should belong to this residue)
	def setPositionFor(self, atom_id, px, py, pz):
		for a in self.atoms:
			if a.getNumber()==int(atom_id):
				a.setPosition(px, py, pz)
				return
		sys.exit("Atom "+atom_id+" does not correspond to any atom of "+self.toString())

	#calculates each coordinate of the residue position as the average of the corresponding coordinate values of the atoms
	def calculatePosition(self):
		coordinate_sum = [ 0, 0, 0 ]
		for a in self.atoms:
			ap = a.getPosition()
			coordinate_sum[0] = coordinate_sum[0] + ap[0]
			coordinate_sum[1] = coordinate_sum[1] + ap[1]
			coordinate_sum[2] = coordinate_sum[2] + ap[2]
		self.position = (coordinate_sum[0]/len(self.atoms), coordinate_sum[1]/len(self.atoms), coordinate_sum[2]/len(self.atoms))
	
	#returns the string descriptor of this residue
	def toString(self):
		return(self.name+" "+self.chain+" "+self.number)

	#returns the string used to refer to this residue in file names
	def getFilenameString(self):
		return(self.name.lower()+self.chain+str(self.number))
	
	#sets the depth index of the atom corresponding to <number> (which should belong to this residue)
	def setDIForAtom(self, number, di):
		for a in self.atoms:
			if a.getNumber()==int(number):
				a.setDepthIndex(di)
				return
		#otherwise, the atom does not belong to this residue. Debug section below allows to detect such cases.
		#DEBUG START
		#sys.exit("Atom "+str(number)+" does not correspond to any atom of "+self.toString())
		#DEBUG STOP

	#returns the residue name (f.i. "TRP" for Triptophan)
	def getName(self):
		return self.name

	#returns the residue number (its position in the protein chain sequence)
	def getNumber(self):
		return self.number

	#returns the depth index of the residue as the average depth index of its atoms
	def getDepthIndex(self):
		DI_sum = 0
		for a in self.atoms:
			DI_sum = DI_sum + a.getDepthIndex()
		return float(DI_sum)/len(self.atoms)
	
	#returns the coordinates in 3D space of this residue
	def getPosition(self):
		if(self.position):
			return self.position
		sys.exit("No position set for Residue "+self.toString())

	#returns the distance between the positions of the two residues
	def getDistanceFrom(self, other_residue):
		return math.sqrt( (self.position[0]-other_residue.position[0])*(self.position[0]-other_residue.position[0]) + (self.position[1]-other_residue.position[1])*(self.position[1]-other_residue.position[1]) + (self.position[2]-other_residue.position[2])*(self.position[2]-other_residue.position[2]) )

	#returns the distance between the alpha-carbons of the two residues
	def getDistanceFromCA(self, other_residue):
		sp = self.getAtom("CA").getPosition()
		op = other_residue.getAtom("CA").getPosition()
		return math.sqrt( (sp[0]-op[0])*(sp[0]-op[0]) + (sp[1]-op[1])*(sp[1]-op[1]) + (sp[2]-op[2])*(sp[2]-op[2]) )
	
	#returns the distance between the position of this residue and the closest point of the other residue's 3D shape (a collection of spheres, one for each residue atom, centered on the atom coordinates and with radius equal to Van Der Walls radius of the atom)
	def getMinDistanceVDW(self, res):
		sp = self.getPosition()
		min_dist = None
		#scan the list of atoms of the other residue
		for a in res.atoms:
			#get atom position
			ap = a.getPosition()
			#calculate distance between this atom and <self>, taking the van der Walls radius of each atom into account
			raw_dist = math.sqrt( (sp[0]-ap[0])*(sp[0]-ap[0]) + (sp[1]-ap[1])*(sp[1]-ap[1]) + (sp[2]-ap[2])*(sp[2]-ap[2]) )	
			net_dist = raw_dist - vdwr_dict[a.getName()] - vdwr_dict["C"]
			#check if this is the minimum distance so far
			if min_dist is None:
				min_dist = net_dist
				continue
			if net_dist < min_dist:
				min_dist = net_dist
		#check that the distance was actually calculated
		if min_dist is None:
			sys.exit("Unable to calculate the distance between "+ self.toString()+" and "+ res.toString())
		#return the final distance value
		return min_dist


	#tells whether at least one atom of this residue has a SASA greater than the specified <minimum_sasa_value>
	def hasExposedAtoms(self, minimum_sasa_value):
		for a in self.atoms:
			if a.getSasa()>=minimum_sasa_value:
				return True
		return False

	#check percentage SASA exposition of this residue
	def checkExposition(self, min_perc_sasa, max_perc_sasa):
		sasa_sum = 0		
		for a in self.atoms:
			sasa_sum += a.getSasa()
		perc_sasa = 100 * sasa_sum / asa_dict[self.name]
		return ( perc_sasa >= min_perc_sasa and perc_sasa <= max_perc_sasa )
	
	#returns this residue atom corresponding to <name> (f.i. "CA" for alpha carbon)
	def getAtom(self, name):
		for a in self.atoms:
			if a.getName()==name:
				return a
		return None

	#check whether the sidechain is pointing to the position <pos>
	def isPointingTo(self, pos):
		#get distance of alpha carbon from pos
		ca_dist = self.getAtom("CA").getDistanceFromCoordinateVector(pos)
		for atom_name in sidechain_dict.get(self.getName()):
			#if at least one sidechain atom is more distant from pos than alpha carbon is, return "False"
			if self.getAtom(atom_name).getDistanceFromCoordinateVector(pos) > ca_dist:
				return False
		#otherwise return "True"		
		return True
