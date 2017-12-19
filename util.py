import numpy as np
import plotly.offline as offline
import plotly.graph_objs as go

def plot_tetrahedra(tetrahedra_list, filename):
	mesh_list = []
	for tetrahedron in tetrahedra_list:
		mesh = go.Mesh3d(x=tetrahedron.v[:,0],
						 y=tetrahedron.v[:,1],
						 z=tetrahedron.v[:,2],
						 i=[0, 0, 0, 1],
						 j=[1, 2, 3, 2],
						 k=[2, 3, 1, 3],
						 color='blue')
		mesh_list.append(mesh)
	layout = go.Layout(xaxis=go.XAxis(title='x'),
					   yaxis=go.YAxis(title='y'))
	data = go.Data(mesh_list)
	fig = go.Figure(data=data, layout=layout)
	offline.plot(fig, filename=filename)
