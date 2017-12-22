import numpy as np
import plotly.offline as offline
import plotly.graph_objs as go
import subprocess
from convexPolygon import Simplex

def plot_tetrahedra(tetrahedra_list, filename):
	mesh_list = []
	for tetrahedron in tetrahedra_list:
		mesh = go.Mesh3d(x=tetrahedron.v[:,0],
						 y=tetrahedron.v[:,1],
						 z=tetrahedron.v[:,2],
						 i=[0, 0, 0, 1],
						 j=[1, 2, 3, 2],
						 k=[2, 3, 1, 3])
		mesh_list.append(mesh)
	layout = go.Layout(xaxis=go.XAxis(title='x'),
					   yaxis=go.YAxis(title='y'))
	data = go.Data(mesh_list)
	fig = go.Figure(data=data, layout=layout)
	offline.plot(fig, filename=filename)

# Determines whether or not two tetrahedra intersect each other
# Returns True if they do, False otherwise
def collision_detection(tetra1, tetra2):
	l = np.linalg.norm(tetra1.v[0]-tetra1.c)
	# If the centroids of 2 tetrahedra are not at least this close
	# then there is no need to perform the expensive check
	if np.linalg.norm(tetra1.c-tetra2.c) > 2*l:
		return False
	# Slow probably because I create a new clojure instance
	# Using a python implementation will likely speed this up
	else:
		command = ["clojure", "tetrahedron-intersect.clj"]
		for i in range(len(tetra1.v)):
			for j in range(len(tetra1.v[i])):
				command.append(str(tetra1.v[i][j]))
		for i in range(len(tetra2.v)):
			for j in range(len(tetra2.v[i])):
				command.append(str(tetra2.v[i][j]))
		collided = subprocess.check_output(command)
		return collided == 'true\n'

def cubeContainerVolume(tetlist):
	corners = [[None,None] for _ in range(len(tetlist[0].c))]
	for tet in tetlist:
		for vertex in tet.v:
			for i in range(len(vertex)):
				if vertex[i] < corners[i][0] or corners[i][0] == None:
					corners[i][0] = vertex[i]
				if vertex[i] > corners[i][1] or corners[i][1] == None:
					corners[i][1] = vertex[i]
	volume = 1
	for edge in corners:
		volume *= (edge[1]-edge[0])
	return volume

# Generates a list of numtetras tetrahedra with random positions and orientations
def generate_tetrahedra(numtetras, l=1.0):
	# Initialize a regular tetrahedron with a vertex at [0,0,0]
	vertices = [[0.0,0.0,0.0],[l/np.sqrt(5),0.0,2*l/np.sqrt(5)],[2*l/np.sqrt(5),l/np.sqrt(5),0.0],[0.0,2*l/np.sqrt(5),l/np.sqrt(5)]]
	tetralist = []
	tetralist.append(Simplex(vertices))
	multiplier = 2

	# Perform random transformations on this tetrahedron until we have numtetras
	for i in range(numtetras-1):
		t = Simplex(vertices)
		choice = np.random.choice(3)
		trans = []
		for j in range(3):
			if j == choice:
				trans.append(np.random.uniform(multiplier*l,(multiplier+1)*l))
			else:
				trans.append(np.random.uniform((multiplier+1)*l))
		t.translate(trans)
		t.Rotate(np.random.uniform(360),(0,1),'deg')
		t.Rotate(np.random.uniform(360),(1,2),'deg')
		t.Rotate(np.random.uniform(360),(0,2),'deg')
		tetralist.append(t)
		multiplier += 2
	return tetralist

def getCOM(tetralist):
	COM = np.zeros(len(tetralist[0].c))
	for tet in tetralist:
		COM += tet.c
	COM /= len(tetralist)
	return COM

def randomizedPacking(numtetras, initstepscale=0.1, stepscalereduction=0.5, initrotrange=360, rotreduction=0.5, stepthreshold=1e-6, rotationthreshold=0.1, l=1.0, max_iters=1000):
	tetralist = generate_tetrahedra(numtetras)
	scales = [initstepscale]*len(tetralist)
	rotranges = [initrotrange]*len(tetralist)
	converged = [False]*len(tetralist)
	centerOfMass = getCOM(tetralist)
	for iteration in range(max_iters):
		print(iteration)
		updatedlist = []
		for i, tet in enumerate(tetralist):
			# preliminary checks
			useTranslation = True
			useRotation = True
			if scales[i] < stepthreshold and rotranges[i] < rotationthreshold:
				updatedlist.append(tet)
				converged[i] = True
				continue

			# get translation
			dist = np.linalg.norm(centerOfMass-tet.c)
			step_size = scales[i]*dist
			tvec = np.random.uniform(-step_size,step_size,len(tet.c))
			if np.linalg.norm(centerOfMass-tet.c-tvec)>dist:
				tvec *= -1
			translatedTet = Simplex(tet.v)
			translatedTet.translate(tvec)
			for j in range(len(tetralist)):
				if j == i:
					continue
				else:
					if collision_detection(translatedTet,tetralist[j]):
						useTranslation = False
						scales[i] *= stepscalereduction
						break

			# get rotation
			degrees = np.random.uniform(rotranges[i])
			plane = (0,np.random.choice(range(1,len(tet.c))))
			rotatedTet = Simplex(tet.v)
			rotatedTet.Rotate(degrees,plane,'deg')
			maxdist = 0
			maxrdist = 0
			for j in range(len(tet.v)):
				d = np.linalg.norm(centerOfMass-tet.v[j])
				rd = np.linalg.norm(centerOfMass-rotatedTet.v[j])
				if d > maxdist:
					maxdist = d
				if rd > maxrdist:
					maxrdist = rd
			if maxrdist > maxdist:
				degrees *= -1
			rotatedTet = Simplex(tet.v)
			rotatedTet.Rotate(degrees,plane,'deg')
			for j in range(len(tetralist)):
				if j == i:
					continue
				else:
					if collision_detection(rotatedTet,tetralist[j]):
						useRotation = False
						rotranges[i] *= rotreduction
						break

			# apply correct update
			updatedTet = Simplex(tet.v)
			if useTranslation:
				updatedTet.translate(tvec)
			if useRotation:
				updatedTet.Rotate(degrees,plane,'deg')
			updatedlist.append(updatedTet)

		tetralist = updatedlist[:]
		if iteration % 10 == 9:
			cubeVolume = cubeContainerVolume(tetralist)
			print "container volume: ", cubeVolume
			tetVolume = l**3.0/(6*2**0.5)*numtetras
			print "density: ", tetVolume/cubeVolume
		if not False in converged:
			break
	return tetralist







#TODO: MC Search, (Maybe Gradient Descent)