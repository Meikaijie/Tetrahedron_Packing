import numpy as np
from simplex import Simplex
from util import plot_tetrahedra

def main():
	vertices = [[0,0,0],[1,0,2],[2,1,0],[0,2,1]]
	tetrahedron1 = Simplex(vertices)
	tetrahedron2 = Simplex(vertices)
	# tetrahedron2.translate([1,1,1])
	tetrahedron2.cRotate(25,(0,1),'deg')
	# print(tetrahedron2.v)
	tetrahedron3 = Simplex(vertices)
	tetrahedron3.translate([3,3,3])
	tetrahedron4 = Simplex(vertices)
	tetrahedron4.translate([5,5,5])
	tetrahedron4.cRotate(25,(0,1),'deg')
	tets = [tetrahedron1,tetrahedron2,tetrahedron4]
	plot_tetrahedra(tets, 'test1')

if __name__ == "__main__":
	main()