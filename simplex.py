import numpy as np
import math

class Simplex(object):
	# vertices should have shape nxd
	# basis should have shape dxd
	def __init__(self, vertices, basis=np.array([[1,0,0],[0,1,0],[0,0,1]])):
		self.v = np.array(vertices)
		self.b = np.array(basis)
		self.c = np.zeros(len(vertices[0]))
		for vertex in vertices:
			for i in range(len(vertex)):
				self.c[i] += vertex[i]/float(len(vertices))

	def translate(self, tvec):
		assert len(tvec) == len(self.v[0]), "Translation dimensionality {0} does not match simplex dimensionilty {1}".format(len(tvec),len(self.v[0]))
		for i in range(len(self.v)):
			for j in range(len(self.v[i])):
				self.v[i][j] += tvec[j]
		for i in range(len(self.c)):
			self.c[i] += tvec[i]

	# Rotate counter-clockwise about the centroid
	# Takes in an angle in degrees or radians, assumes radians if no unit argument is included,
	# and a 2d rotation plane specified by tuple or list containing 2 integers
	# which represent the two axes defining the plane
	def cRotate(self, angle, plane, unit='rad'):
		if unit == 'deg':
			angle = math.radians(angle)
		cos = math.cos(angle)
		sin = math.sin(angle)
		rotation_matrix = np.zeros((len(self.v[0]),(len(self.v[0]))))
		for i in range(len(rotation_matrix)):
			for j in range(len(rotation_matrix[i])):
				if i == j:
					if i == plane[0] or i == plane[1]:
						rotation_matrix[i][j] = cos
					else:
						rotation_matrix[i][j] = 1.0
				else:
					if i == plane[0] and j == plane[1]:
						rotation_matrix[i][j] = -sin
					elif i == plane[1] and j == plane[0]:
						rotation_matrix[i][j] = sin
		print(rotation_matrix)
		offset = self.c - np.dot(rotation_matrix,self.c)
		for i in range(len(self.v)):
			self.v[i] = np.dot(rotation_matrix,self.v[i]) #+ offset

	# def checkRegular(self):
	# 	total = 0.0
	# 	edges = 0.0
	# 	isRegular = True
	# 	for i in range(len(self.v)):
	# 		for j in range(i+1,len(self.v)):
	# 			total += np.linalg.norm(self.v[i]-self.v[j])
	# 			# print(np.linalg.norm(self.v[i]-self.v[j]))
	# 			edges += 1
	# 	for i in range(len(self.v)):
	# 		for j in range(i+1,len(self.v)):
	# 			if math.fabs(np.linalg.norm(self.v[i]-self.v[j]) - total/edges) > 0.000000001:
	# 				# print(math.fabs(np.linalg.norm(self.v[i]-self.v[j]) - total/edges))
	# 				isRegular = False
	# 	if isRegular:
	# 		print("All Good!")
	# 	else:
	# 		print("Bad Transformation!")









