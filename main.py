import argparse
import util

def parse_args():
	parser = argparse.ArgumentParser(description="Tetrahedron Packing Simulation")
	parser.add_argument("-n", "--numtetras", type=int, help="number of tetrahedra to pack, valid at or above 1", default=4)
	parser.add_argument("-is", "--initstepscale", type=float, help="initial step size coefficient, valid between 0 and 1", default=0.9)
	parser.add_argument("-sr", "--stepscalereduction", type=float, help="step scale decay rate, valid between 0 and 1", default=0.999)
	parser.add_argument("-st", "--stepthreshold", type=float, help="minimum allowed step scale, valid above 0", default=1e-6)
	parser.add_argument("-ir", "--initrotationrange", type=float, help="initial permitted range of a single rotation move", default=360.0)
	parser.add_argument("-rr", "--rotreduction", type=float, help="rotation range decay rate, valid between 0 and 1", default=0.999)
	parser.add_argument("-rt", "--rotationthreshold", type=float, help="minimum allowed rotation range", default=0.1)
	parser.add_argument("-l", "--length", type=float, help="edge length of the tetrahedra", default=1.0)
	parser.add_argument("-i", "--iterations", type=int, help="max number of iterations", default=5000)
	parser.add_argument("-v", "--verbose", action='store_true')
	parser.add_argument("-b", "--bias", help="induce nearest neighbor attraction bias between tetrahedra to avoid dogpile clustering", action='store_true')
	parser.add_argument("-np", "--noplot", help="disable plotting tetrahedra with plotly", action='store_true', default=False)
	return parser.parse_args()

def main():
	args = parse_args()
	packedTets = None
	filename = "graph_n"+str(args.numtetras)+"_is"+str(args.initstepscale)+"_sr"+str(args.stepscalereduction)+"_st"+str(args.stepthreshold)+"_ir"+str(args.initrotationrange)+"_rr"+str(args.rotreduction)+"_rt"+str(args.rotationthreshold)+"_l"+str(args.length)+"_i"+str(args.iterations)
	if args.bias:
		filename += "_withbias"
		packedTets = util.randomizedGuidedPackingV3(args.numtetras, initstepscale=args.initstepscale, stepscalereduction=args.stepscalereduction,
					initrotrange=args.initrotationrange, rotreduction=args.rotreduction, stepthreshold=args.stepthreshold,
					rotationthreshold=args.rotationthreshold, l=args.length, max_iters=args.iterations, verbose=args.verbose)
	else:
		packedTets = util.randomizedGuidedPackingV2(args.numtetras, initstepscale=args.initstepscale, stepscalereduction=args.stepscalereduction,
					initrotrange=args.initrotationrange, rotreduction=args.rotreduction, stepthreshold=args.stepthreshold,
					rotationthreshold=args.rotationthreshold, l=args.length, max_iters=args.iterations, verbose=args.verbose)
	output = ""
	with open("packing_output.txt","wb") as f:
		for basis in packedTets[0].b:
			line = ""
			for coord in basis:
				line += str(coord)+" "
			line += "\n"
			f.write(line)
			output += line
		for tet in packedTets:
			line = ""
			for vertex in tet.v:
				for coord in vertex:
					line += str(coord)+" "
			line += "\n"
			f.write(line)
			output += line
	print(output)
	if not args.noplot:
		util.plot_tetrahedra(packedTets, filename+".html")


if __name__ == "__main__":
	main()