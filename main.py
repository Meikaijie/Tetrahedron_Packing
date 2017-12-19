import numpy as np
from simplex import Simplex
from util import plot_tetrahedra

def main():
	vertices = [[0,0,0],[1,0,2],[2,1,0],[0,2,1]]
	tetrahedron1 = Simplex(vertices)
	tetrahedron2 = Simplex(vertices)
	tetrahedron2.cRotate(360,(0,1),'deg')
	tetrahedron3 = Simplex(vertices)
	tetrahedron3.translate([3,3,3])
	tetrahedron4 = Simplex(vertices)
	tetrahedron4.translate([3,3,3])
	tetrahedron4.cRotate(90,(0,1),'deg')
	tets = [tetrahedron1,tetrahedron2,tetrahedron3,tetrahedron4]
	plot_tetrahedra(tets, 'test1')

if __name__ == "__main__":
	main()