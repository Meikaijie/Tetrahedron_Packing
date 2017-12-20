import numpy as np
import math

class Simplex(object):
	# vertices should have shape nxd
	# basis should have shape dxd
	def __init__(self, vertices, basis=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])):
		v = map(lambda x: map(lambda y: float(y),x),vertices)
		self.v = np.array(v)
		b = map(lambda x: map(lambda y: float(y),x),basis)
		self.b = np.array(b)
		self.c = np.zeros(len(vertices[0]))
		for vertex in vertices:
			for i in range(len(vertex)):
				self.c[i] += vertex[i]/float(len(vertices))

	# Translate the tetrahedron by performing element-wise addition of tvec to
	# each vertex in self.v
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
		if math.fabs(cos) < 1e-14:
			cos = 0.0
		sin = math.sin(angle)
		if math.fabs(sin) < 1e-14:
			sin = 0.0
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
		# print(rotation_matrix)
		offset = self.c - np.dot(rotation_matrix,self.c)
		for i in range(len(self.v)):
			self.v[i] = np.dot(rotation_matrix,self.v[i]) + offset

	# Print vertices for debugging
	def to_string(self):
		s = ""
		for i in range(len(self.v)):
			for j in range(len(self.v[i])):
				s += str(self.v[i][j]) + ' '
		return s[:-1]







