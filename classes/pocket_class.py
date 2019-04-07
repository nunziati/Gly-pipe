#coding=utf-8
import sys
import os
import math

#import atom class
from atom_class import Atom

#pocket object class
class Pocket:
	
	#constructor
	def __init__(self, pocket_id):
		self.id = int(pocket_id)  #pocket identification number
		self.dictionary = {}      #dictionary of pocket descriptors
		self.vertices = None      #list of vertices (pseudoatoms)

	#insert a key, value pair in the descriptor dictionary	
	def insert(self, key, value):
		self.dictionary[key] = value

	#add vertex
	def addVertex(self, xx, yy, zz):
		if self.vertices is None:
			self.vertices = [Atom(0.0, "STP", str(1))]
		else:
			self.vertices.append(Atom(0.0, "STP", str(len(self.vertices)+1)))
		self.vertices[len(self.vertices)-1].setPosition(xx, yy, zz)

	#returns the center of this pocket in 3D space
	def getCenter(self):	
		x_sum = 0.0
		y_sum = 0.0
		z_sum = 0.0
		for v in self.vertices:
			v_pos = v.getPosition()
			x_sum = x_sum + v_pos[0]
			y_sum = y_sum + v_pos[1]
			z_sum = z_sum + v_pos[2]
		return [x_sum/len(self.vertices), y_sum/len(self.vertices), z_sum/len(vertices)]

	#returns the minimum euclidean distance in 3D space between <point> and this pocket
	def getMinDistance(self, point):
		if self.vertices is None:
			sys.exit("Error: cannot calculate distances from pocket "+self.getID()+" as its vertices are not set")
		min_dist = None
		for v in self.vertices:
			vp = v.getPosition()
			dist = math.sqrt((vp[0]-point[0])*(vp[0]-point[0])+(vp[1]-point[1])*(vp[1]-point[1])+(vp[2]-point[2])*(vp[2]-point[2]))
			if min_dist is None:
				min_dist = dist
				continue
			if dist < min_dist:
				min_dist = dist
		return min_dist

	#returns the pocket id
	def getID(self):
		return self.id
	
	#returns the descriptor dictionary
	def getDictionary(self):
		return self.dictionary

	#returns the value corresponding to key, generates a "key error" if key is not found
	def getValue(self, key):
		return self.dictionary.get(key)

	#returns the position of this pocket
	def getPosition(self):
		return self.position

	#returns the pocket's feature vector, where the feature labels are those specified in <labels>
	def getFeatureVector(self, labels):
		fv = []
		for l in labels:
			fv.append(self.dictionary.get(l))
		return fv

	#returns a string representation of the pocket descriptors
	def toString(self):
		des_str = ""
		for k in self.dictionary.keys():
			des_str = des_str + k + " : " + str(self.dictionary.get(k))+"\n"
		return des_str

	#just like "toString", but only for the keys in input
	def keysToString(self, input_keys):
		des_str = ""
		for k in input_keys:
			des_str = des_str + k + " : " + str(self.dictionary.get(k))+"\n"
		return des_str
