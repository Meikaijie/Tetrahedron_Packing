import numpy as np
from convexPolygon import Simplex
from util import plot_tetrahedra, collision_detection

def main():
	vertices = [[0,0,0],[1,0,2],[2,1,0],[0,2,1]]
	tetrahedron1 = Simplex(vertices)
	tetrahedron2 = Simplex(vertices)
	tetrahedron2.cRotate(25,(0,1),'deg')
	tetrahedron3 = Simplex(vertices)
	tetrahedron3.translate([1.1,1.1,1.1])
	tetrahedron4 = Simplex(vertices)
	tetrahedron4.translate([5,5,5])
	tetrahedron4.cRotate(25,(0,1),'deg')
	tets = [tetrahedron1,tetrahedron2,tetrahedron3,tetrahedron4]
	plot_tetrahedra(tets, 'test1')
	print collision_detection(tetrahedron1, tetrahedron2)
	print collision_detection(tetrahedron1, tetrahedron3)

if __name__ == "__main__":
	main()