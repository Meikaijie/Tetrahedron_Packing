import numpy as np
import plotly.offline as offline
import plotly.graph_objs as go
import subprocess

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

def collision_detection(tetra1, tetra2):
	l = np.linalg.norm(tetra1.v[0]-tetra1.c)
	# If the centroids of 2 tetrahedra are not at least this close
	# then there is no need to perform the expensive check
	if np.linalg.norm(tetra1.c-tetra2.c) > 2*l:
		return False
	# Slow probably because I create a new clojure process
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

#TODO: MC Search, (Maybe Gradient Descent)