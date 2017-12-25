import numpy as np
from convexPolygon import Simplex
import util

def main():
	# Single Test
	packedTets = util.randomizedGuidedPacking(4, genmult=2, geninc=2, initstepscale=0.2, stepscalereduction=0.95, rotreduction=0.95, max_iters=1000)
	util.plot_tetrahedra(packedTets, 'test_tets4_is2_ir360_sr95_rr95_iter1000')

	# Aggregate statistic test
	# bestTets = None
	# bestDensity = 0
	# densities = []
	# for _ in range(100):
	# 	packedTets = util.randomizedGuidedPacking(2, genmult=2, geninc=2, initstepscale=0.9, stepscalereduction=0.95, rotreduction=0.95, max_iters=1000)
	# 	density = 1.0/(6*2**0.5)*2/util.cubeContainerVolume(packedTets)
	# 	densities.append(density)
	# 	if density > bestDensity:
	# 		bestTets = packedTets
	# util.plot_tetrahedra(bestTets, 'test_tets2_is2_ir360_sr95_rr95_iter1000')
	# print(sorted(densities))

if __name__ == "__main__":
	main()