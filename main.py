import numpy as np
from convexPolygon import Simplex
import util

def main():
	# vertices = [[0,0,0],[1,0,2],[2,1,0],[0,2,1]]
	# tetrahedron1 = Simplex(vertices)
	# tetrahedron2 = Simplex(vertices)
	# tetrahedron2.Rotate(25,(0,1),'deg')
	# # tetrahedron2.translate([3,3,3])
	# tetrahedron2.Rotate(90,(1,2),'deg')
	# tetrahedron2.Rotate(180,(0,2),'deg')
	# tetrahedron3 = Simplex(vertices)
	# tetrahedron3.translate([1.1,1.1,1.1])
	# tetrahedron4 = Simplex(vertices)
	# tetrahedron4.translate([5,5,5])
	# tetrahedron4.Rotate(25,(0,1),'deg')
	# tets = [tetrahedron1,tetrahedron2,tetrahedron3,tetrahedron4]
	# util.plot_tetrahedra(tets, 'test1')
	# print util.collision_detection(tetrahedron1, tetrahedron2)
	# print util.collision_detection(tetrahedron1, tetrahedron3)
	# tets2 = util.generate_tetrahedra(100)
	# util.plot_tetrahedra(tets2, 'test2')
	packedTets = util.randomizedPacking(4, initstepscale=0.9, stepscalereduction=0.95, rotreduction=0.95, max_iters=1000)
	util.plot_tetrahedra(packedTets, 'test_is9_ir360_sr95_rr95_iter1000')

if __name__ == "__main__":
	main()