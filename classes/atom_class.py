#coding=utf-8
import sys
import os
import math

#atom object class
class Atom:
	
	#constructor
	def __init__(self, sasa, name, number):
		self.sasa = float(sasa)   #Solvent Accessible Surface Area
		self.di = -1              #Depth Index
		self.name = name          #atom name (f.i. CA -> alpha carbon)
		self.number = int(number) #atom sequential number
		self.position = None      #atom position in 3D space
	
	#sets the depth index of this atom
	def setDepthIndex(self, depth_index):
		self.di = float(depth_index)

	#sets the coordinates in 3D space of this atom
	def setPosition(self, px, py, pz):
		self.position = [px, py, pz]	

	#returns the solvent accessible surface area
	def getSasa(self):
		return self.sasa

	#returns the depth index
	def getDepthIndex(self):
		return self.di

	#returns the atom name (f.i. CA -> alpha carbon)
	def getName(self):
		return self.name
	
	#returns the atom number (global identifier)
	def getNumber(self):
		return self.number
	
	#returns the coordinates in 3D space of this atom
	def getPosition(self):
		return self.position

	#returns the distance between <self.position> and the coordinates stored in <cv>
	def getDistanceFromCoordinateVector(self, cv):
		return math.sqrt( (self.position[0]-cv[0])*(self.position[0]-cv[0]) + (self.position[1]-cv[1])*(self.position[1]-cv[1]) + (self.position[2]-cv[2])*(self.position[2]-cv[2]) )
