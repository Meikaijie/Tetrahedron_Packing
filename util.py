import numpy as np
import plotly.offline as offline
import plotly.graph_objs as go
import subprocess
import time
from convexPolygon import Simplex

# Plots a group of tetrahedra in a 3-D cartesian coordinate space
def plot_tetrahedra(tetrahedra_list, filename):
	mesh_list = []
	# Create plotly tetrahedra meshes
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
# Takes in 2 3-D simplex objects and a subprocess running a clojure REPL
# Returns True if they do, False otherwise
def collision_detection(tetra1, tetra2, subproc):
	l = np.linalg.norm(tetra1.v[0]-tetra1.c)
	# If the centroids of 2 tetrahedra are not at least this close
	# then there is no need to perform the expensive check
	if np.linalg.norm(tetra1.c-tetra2.c) > 2*l:
		return False
	else:
		# Uses a persistent subprocess for many repeated collision detections
		stringv1 = tetra1.to_string()
		stringv2 = tetra2.to_string()
		command = '(prn (intersect-tetrahedra? '+stringv1+' '+stringv2+'))'
		subproc.stdin.write(command)
		subproc.stdin.write('(flush)')
		subproc.stdout.readline() # flush out empty line
		collided = subproc.stdout.readline()
		subproc.stdout.readline() # flush out empty line

		# Preliminary integration of provided clojure implementation for
		# tetrahedron overlap/collision detection
		# Run collision detection with new clojure REPL
		# closes after algorithm is done
		# very slow due to REPL boot time
		# //////////////////////////////////////////////////
		# command = ["clojure", "tetrahedron-intersect.clj"]
		# for i in range(len(tetra1.v)):
		# 	for j in range(len(tetra1.v[i])):
		# 		command.append(str(tetra1.v[i][j]))
		# for i in range(len(tetra2.v)):
		# 	for j in range(len(tetra2.v[i])):
		# 		command.append(str(tetra2.v[i][j]))
		# collided = subprocess.check_output(command)
		# //////////////////////////////////////////////////

		return collided == 'true\n'

# Get volume of smallest axis aligned rectangular prism containing all tetrahedra
# Used to observe convergence; not a good measure of packing quality because
# the box is axis aligned meaning the volume is greatly affected by tetrahedra orientations
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
# initialmultiplier determines the minimum tetrahedron lengths away that two tetrahedra can be from each other
# multiplierincrement determines how much to increase the possible placement space per additional tetrahedron
# l determines the desired edge length of all generated tetrahedra
def generate_tetrahedra(numtetras, initialmultiplier=2, multiplierincrement=2, l=1.0):
	# Initialize a regular tetrahedron with a vertex at [0,0,0]
	vertices = [[0.0,0.0,0.0],[l/np.sqrt(5),0.0,2*l/np.sqrt(5)],[2*l/np.sqrt(5),l/np.sqrt(5),0.0],[0.0,2*l/np.sqrt(5),l/np.sqrt(5)]]
	tetralist = []
	tetralist.append(Simplex(vertices))
	multiplier = initialmultiplier

	# Perform random transformations on copies of this tetrahedron until we have numtetras
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
		multiplier += multiplierincrement
	return tetralist

# Get "center of mass" of all tetrahedra,
# calculated as the average of all tetrahedra centroid positions
def getCOM(tetralist):
	COM = np.zeros(len(tetralist[0].c))
	for tet in tetralist:
		COM += tet.c
	COM /= len(tetralist)
	return COM

# Packing algorithm that randomly initializes numtetras tetrahedra in space
# then compresses them by performing random translations and rotations that 
# move the tetrahedra towards the collective "center of mass"
def randomizedGuidedPacking(numtetras, genmult=2, geninc=2, initstepscale=0.9, stepscalereduction=0.95, initrotrange=360, rotreduction=0.95, stepthreshold=1e-6, rotationthreshold=0.1, l=1.0, max_iters=1000, verbose=False):
	assert(initstepscale>0 and initstepscale<=1)
	assert(stepscalereduction>0 and stepscalereduction<1)
	assert(stepthreshold>0)
	assert(rotreduction>0 and rotreduction<1)

	# Can be removed if python implementation of collision detection is used
	subproc = subprocess.Popen('clojure', stdout=subprocess.PIPE, stdin=subprocess.PIPE)
	subproc.stdin.write('(clojure.main/load-script "tetrahedron-intersect.clj")')
	subproc.stdin.write('(flush)')
	subproc.stdout.readline() # flush out useless line
	subproc.stdout.readline() # flush out useless line

	tetralist = generate_tetrahedra(numtetras,genmult,geninc,l)
	scales = [initstepscale]*len(tetralist)
	rotranges = [initrotrange]*len(tetralist)
	converged = [False]*len(tetralist)
	centerOfMass = getCOM(tetralist)
	for iteration in range(max_iters):
		updatedlist = []
		for i, tet in enumerate(tetralist):
			# preliminary convergence check
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
			# flip diretion of translation if it moves the tetrahedron away from COM
			if np.linalg.norm(centerOfMass-tet.c-tvec)>dist:
				tvec *= -1
			translatedTet = Simplex(tet.v)
			translatedTet.translate(tvec)
			# reject translation if it results in a collision and scale back step size
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					if collision_detection(translatedTet,updatedlist[j],subproc):
						useTranslation = False
						scales[i] *= stepscalereduction
						break
				else:
					if collision_detection(translatedTet,tetralist[j],subproc):
						useTranslation = False
						scales[i] *= stepscalereduction
						break

			tetUsed = tet
			if useTranslation:
				tetUsed = translatedTet

			# get rotation
			degrees = np.random.uniform(rotranges[i])
			plane = (0,np.random.choice(range(1,len(tetUsed.c))))
			rotatedTet = Simplex(tetUsed.v)
			rotatedTet.Rotate(degrees,plane,'deg')
			netdist = 0
			netrdist = 0
			for j in range(len(tetUsed.v)):
				d = np.linalg.norm(centerOfMass-tetUsed.v[j])
				rd = np.linalg.norm(centerOfMass-rotatedTet.v[j])
				netdist += d
				netrdist += rd
			# flip direction of rotation if it increases net vertex distance from COM
			if netrdist > netdist:
				degrees *= -1
			rotatedTet = Simplex(tetUsed.v)
			rotatedTet.Rotate(degrees,plane,'deg')
			# reject rotation if it results in a collision and scale back rotation size
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					if collision_detection(rotatedTet,updatedlist[j],subproc):
						useRotation = False
						rotranges[i] *= rotreduction
						break
				else:
					if collision_detection(rotatedTet,tetralist[j],subproc):
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

			# update center of mass
			if i != len(tetralist)-1:
				centerOfMass = getCOM(updatedlist+tetralist[i+1:])
			else:
				centerOfMass = getCOM(updatedlist)

		# update tetralist, print density calculation, check convergence condition
		tetralist = updatedlist[:]
		if verbose and iteration % 10 == 9:
			print("Iteration: " + str(iteration+1))
			cubeVolume = cubeContainerVolume(tetralist)
			print "box container volume: ", cubeVolume
			tetVolume = l**3.0/(6*2**0.5)*numtetras
			print "density: ", tetVolume/cubeVolume
		if not False in converged:
			break
	return tetralist

# Packing algorithm that randomly initializes numtetras tetrahedra in space
# then compresses them by performing random translations and rotations that 
# move the tetrahedra towards the collective "center of mass"
# V2 adds more dynamic rotations
def randomizedGuidedPackingV2(numtetras, genmult=2, geninc=2, initstepscale=0.9, stepscalereduction=0.95, initrotrange=360, rotreduction=0.95, stepthreshold=1e-6, rotationthreshold=0.1, l=1.0, max_iters=1000, verbose=False):
	assert(initstepscale>0 and initstepscale<=1)
	assert(stepscalereduction>0 and stepscalereduction<1)
	assert(stepthreshold>0)
	assert(rotreduction>0 and rotreduction<1)

	# Can be removed if python implementation of collision detection is used
	subproc = subprocess.Popen('clojure', stdout=subprocess.PIPE, stdin=subprocess.PIPE)
	subproc.stdin.write('(clojure.main/load-script "tetrahedron-intersect.clj")')
	subproc.stdin.write('(flush)')
	subproc.stdout.readline() # flush out useless line
	subproc.stdout.readline() # flush out useless line

	tetralist = generate_tetrahedra(numtetras,genmult,geninc,l)
	scales = [initstepscale]*len(tetralist)
	rotranges = [initrotrange]*len(tetralist)
	converged = [False]*len(tetralist)
	centerOfMass = getCOM(tetralist)
	for iteration in range(max_iters):
		updatedlist = []
		for i, tet in enumerate(tetralist):
			# preliminary convergence check
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
			# flip diretion of translation if it moves the tetrahedron away from COM
			if np.linalg.norm(centerOfMass-tet.c-tvec)>dist:
				tvec *= -1
			translatedTet = Simplex(tet.v)
			translatedTet.translate(tvec)
			# reject translation if it results in a collision and scale back step size
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					if collision_detection(translatedTet,updatedlist[j],subproc):
						useTranslation = False
						scales[i] *= stepscalereduction
						break
				else:
					if collision_detection(translatedTet,tetralist[j],subproc):
						useTranslation = False
						scales[i] *= stepscalereduction
						break

			tetUsed = tet
			if useTranslation:
				tetUsed = translatedTet

			# get rotation
			degreeslist = []
			planelist = []
			for j in range(len(tet.c)):
				for k in range(j+1,len(tet.c)):
					degrees = np.random.uniform(rotranges[i])
					plane = (j,k)
					rotatedTet = Simplex(tetUsed.v)
					rotatedTet.Rotate(degrees,plane,'deg')
					netdist = 0
					netrdist = 0
					for n in range(len(tetUsed.v)):
						d = np.linalg.norm(centerOfMass-tetUsed.v[n])
						rd = np.linalg.norm(centerOfMass-rotatedTet.v[n])
						netdist += d
						netrdist += rd
					# flip direction of rotation if it increases net vertex distance from COM
					if netrdist > netdist:
						degrees *= -1
					degreeslist.append(degrees)
					planelist.append(plane)
			rotatedTet = Simplex(tetUsed.v)
			for j in range(len(degreeslist)):
				rotatedTet.Rotate(degreeslist[j],planelist[j],'deg')
			# reject rotation if it results in a collision and scale back rotation size
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					if collision_detection(rotatedTet,updatedlist[j],subproc):
						useRotation = False
						rotranges[i] *= rotreduction
						break
				else:
					if collision_detection(rotatedTet,tetralist[j],subproc):
						useRotation = False
						rotranges[i] *= rotreduction
						break

			# apply correct update
			updatedTet = Simplex(tet.v)
			if useTranslation:
				updatedTet.translate(tvec)
			if useRotation:
				for j in range(len(degreeslist)):
					updatedTet.Rotate(degreeslist[j],planelist[j],'deg')
			updatedlist.append(updatedTet)

			# update center of mass
			if i != len(tetralist)-1:
				centerOfMass = getCOM(updatedlist+tetralist[i+1:])
			else:
				centerOfMass = getCOM(updatedlist)

		# update tetralist, print density calculation, check convergence condition
		tetralist = updatedlist[:]
		if verbose and iteration % 10 == 9:
			print("Iteration: " + str(iteration+1))
			cubeVolume = cubeContainerVolume(tetralist)
			print "box container volume: ", cubeVolume
			tetVolume = l**3.0/(6*2**0.5)*numtetras
			print "density: ", tetVolume/cubeVolume
		if not False in converged:
			break
	return tetralist

# Packing algorithm that randomly initializes numtetras tetrahedra in space
# then compresses them by performing random translations and rotations that 
# move the tetrahedra towards the collective "center of mass"
# V3 adds local optimization bias to potentially reduce "dog piling" in center
def randomizedGuidedPackingV3(numtetras, genmult=2, geninc=2, initstepscale=0.9, stepscalereduction=0.95, initrotrange=360.0, rotreduction=0.95, stepthreshold=1e-6, rotationthreshold=0.1, l=1.0, max_iters=1000, verbose=False):
	assert(initstepscale>0 and initstepscale<=1)
	assert(stepscalereduction>0 and stepscalereduction<1)
	assert(stepthreshold>0)
	assert(rotreduction>0 and rotreduction<1)

	# Can be removed if python implementation of collision detection is used
	subproc = subprocess.Popen('clojure', stdout=subprocess.PIPE, stdin=subprocess.PIPE)
	subproc.stdin.write('(clojure.main/load-script "tetrahedron-intersect.clj")')
	subproc.stdin.write('(flush)')
	subproc.stdout.readline() # flush out useless line
	subproc.stdout.readline() # flush out useless line

	tetralist = generate_tetrahedra(numtetras,genmult,geninc,l)
	scales = [initstepscale]*len(tetralist)
	rotranges = [initrotrange]*len(tetralist)
	converged = [False]*len(tetralist)
	centerOfMass = getCOM(tetralist)
	for iteration in range(max_iters):
		updatedlist = []
		for i, tet in enumerate(tetralist):
			# preliminary convergence check
			useTranslation = True
			useRotation = True
			if scales[i] < stepthreshold and abs(rotranges[i])< rotationthreshold:
				updatedlist.append(tet)
				converged[i] = True
				continue

			# get translation
			dist = np.linalg.norm(centerOfMass-tet.c)
			step_size = scales[i]*dist
			tvec = np.random.uniform(-step_size,step_size,len(tet.c))
			# flip diretion of translation if it moves the tetrahedron away from COM
			if np.linalg.norm(centerOfMass-tet.c-tvec)>dist:
				tvec *= -1
			# add local optimization bias
			nearest = None
			minDist = None
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					d = np.linalg.norm(updatedlist[j].c-tet.c)
					if minDist == None or d<minDist:
						minDist = d
						nearest = updatedlist[j]
				else:
					d = np.linalg.norm(tetralist[j].c-tet.c)
					if minDist == None or d<minDist:
						minDist = d
						nearest = tetralist[j]
			# if only one tetrahedron
			if nearest == None:
				nearest = tet
			bias = nearest.c - tet.c
			for j in range(len(bias)):
				tvec[j] += scales[i]*np.random.uniform(min(bias[j],0),max(bias[j],0))
			translatedTet = Simplex(tet.v)
			translatedTet.translate(tvec)
			# reject translation if it results in a collision and scale back step size
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					if collision_detection(translatedTet,updatedlist[j],subproc):
						useTranslation = False
						scales[i] *= stepscalereduction
						break
				else:
					if collision_detection(translatedTet,tetralist[j],subproc):
						useTranslation = False
						scales[i] *= stepscalereduction
						break

			tetUsed = tet
			if useTranslation:
				tetUsed = translatedTet
			# fine nearest after possible translation
			nearest = None
			minDist = None
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					d = np.linalg.norm(updatedlist[j].c-tetUsed.c)
					if minDist == None or d<minDist:
						minDist = d
						nearest = updatedlist[j]
				else:
					d = np.linalg.norm(tetralist[j].c-tetUsed.c)
					if minDist == None or d<minDist:
						minDist = d
						nearest = tetralist[j]
			# if only one tetrahedron
			if nearest == None:
				nearest = tetUsed

			# get rotation
			degreeslist = []
			planelist = []
			for j in range(len(tet.c)):
				for k in range(j+1,len(tet.c)):
					degrees = np.random.uniform(rotranges[i])
					plane = (j,k)
					rotatedTet = Simplex(tetUsed.v)
					rotatedTet.Rotate(degrees,plane,'deg')
					netdist1 = 0
					netrdist1 = 0
					netdist2 = 0
					netrdist2 = 0
					for n in range(len(tetUsed.v)):
						d1 = np.linalg.norm(centerOfMass-tetUsed.v[n])
						rd1 = np.linalg.norm(centerOfMass-rotatedTet.v[n])
						d2 = np.linalg.norm(nearest.c-tetUsed.v[n])
						rd2 = np.linalg.norm(nearest.c-rotatedTet.v[n])
						netdist1 += d1
						netrdist1 += rd1
						netdist2 += d2
						netrdist2 += rd2
					# flip direction of rotation if it increases net vertex distance from COM
					if netrdist1 > netdist1 and netrdist2 > netrdist2:
						degrees *= -1
					degreeslist.append(degrees)
					planelist.append(plane)

			rotatedTet = Simplex(tetUsed.v)
			for j in range(len(degreeslist)):
				rotatedTet.Rotate(degreeslist[j],planelist[j],'deg')
			# reject rotation if it results in a collision and scale back rotation size
			for j in range(len(tetralist)):
				if j == i:
					continue
				elif j < i:
					if collision_detection(rotatedTet,updatedlist[j],subproc):
						useRotation = False
						rotranges[i] *= rotreduction
						break
				else:
					if collision_detection(rotatedTet,tetralist[j],subproc):
						useRotation = False
						rotranges[i] *= rotreduction
						break

			# apply correct update
			updatedTet = Simplex(tet.v)
			if useTranslation:
				updatedTet.translate(tvec)
			if useRotation:
				for j in range(len(degreeslist)):
					updatedTet.Rotate(degreeslist[j],planelist[j],'deg')
			updatedlist.append(updatedTet)

			# update center of mass
			if i != len(tetralist)-1:
				centerOfMass = getCOM(updatedlist+tetralist[i+1:])
			else:
				centerOfMass = getCOM(updatedlist)

		# update tetralist, print density calculation, check convergence condition
		tetralist = updatedlist[:]
		if verbose and iteration % 10 == 9:
			print("Iteration: " + str(iteration+1))
			cubeVolume = cubeContainerVolume(tetralist)
			print "box container volume: ", cubeVolume
			tetVolume = l**3.0/(6*2**0.5)*numtetras
			print "density: ", tetVolume/cubeVolume
		if not False in converged:
			break
	return tetralist