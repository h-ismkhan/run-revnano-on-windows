""" Naive simple algorithm to recover staple routings

A staple is run along the scaffold and split at the longest
matching sites. The first longest match is used in each case

Performance to be contrasted with full REVNANO algorithm.

The naive algorithm is DETERMINISTIC and has NO FREE PARAMETERS
Hence, it just needs to be run once, per origami

Ben Shirt-Ediss, 2023
"""

import sys
sys.path.append("../")

import os

from REVNANO import dna
from REVNANO import revnano
from REVNANO import parameters
from REVNANO import stage3
from REVNANO import stage4

from contactmap import iohandler as cio






#
#
# NAIVE STAPLE ROUTING ALGORITHM
#
#

def find_longest_match(subseq53, origami, MU_MIN) :
	""" Returns (i35, runlen) of the longest match, or False if none

		5------------>3 staple
	3<---------------------5 scaffold
		^i35
	"""

	run = dna.get_run(origami, subseq53)
		# i53: number of consecutive matching sites

	mu_max = max(run.values())

	if mu_max < MU_MIN :
		return False 						# all longest matches too short

	for i35, runlen in run.items() :
		if runlen == mu_max :				# quit on finding first longest match
			break

	return i35, runlen




def route_staple(staple_id, origami, MU_MIN) :
	""" Returns staple route, from staple 5' to staple 3'
	"""

	# note: staples which self-cross are not accounted for

	subseq53 = origami["staples53"][staple_id]

	route = [] 	# list of (i53,runlen) tuples
	
	while True :

		lm = find_longest_match(subseq53, origami, MU_MIN)
		
		if lm != False :
			# longest match does exist for this staple subsection
			i35, runlen = lm
			i53 = origami["scaffold_nt"] - 1 - i35
			route.append( [i53,runlen] )

			# prepare next subsection of staple for matching
			subseq53 = subseq53[runlen:]

			if len(subseq53) == 0 :
				break 					# all of staple is matched
		else :
			return False  				# no matches for this finite staple subsection
										# therefore, staple cannot be routed
	return route




def route_all_staples(origami, MU_MIN) :
	# note: origami dict is updated by reference 
	# (origami["footprints"] and origami["staples_placed"])

	print("Placing staples by longest sequence alignment...")

	for staple_id in range(0,origami["num_staples"]) :
		
		route = route_staple(staple_id, origami, MU_MIN)

		if route == False :		# staple could not be routed
			print("--> staple %d could not be routed" % staple_id)
			continue

		# update St dict with staple route
		origami["staples_placed"][staple_id] = route

		# update footprints dict with staple route
		for section in route :
			i53s, runlen = section			
			for i53r in range(i53s, i53s-runlen, -1) :
				i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53r)
				origami["footprints"][i53].add(staple_id)

	print("[Done]\n", flush=True)























#
#
# EXECUTE
#
#

def reverse_engineer(origami, MU_MIN, verbose = False) :

	if verbose == False :
		sys.stdout = open(os.devnull, 'w')	# block print statements

	route_all_staples(origami, MU_MIN)

	# now run revnano stage 3 and 4
	# to fix overlaps, recompose staples, and make the final contact map

	try :	
		stage3.process(origami)
		contactmap = stage4.process(origami)

	except ValueError as errormsg :
		if verbose == False : sys.stdout = sys.__stdout__
		return None, str(errormsg) 	# without string cast, errormsg is a ValueError object

	if verbose == False : sys.stdout = sys.__stdout__	# re-enable print statements
	return contactmap, None




if __name__ == "__main__" :

	origami_name = str(sys.argv[1]).replace(".rev", "")

	print("\n")
	print("-----------------------------------------------------------")
	print("NAIVE reverse engineering: %s" % origami_name)
	print("-----------------------------------------------------------")
	print("")

	origami = revnano.read_rev_file(origami_name, "_assets/")

	contactmap, errormsg = reverse_engineer(origami, parameters.MU_MIN)

	if errormsg != None :
		print("\n%s\n" % errormsg)
		quit()

	# output stats
	c1, staples_placed = origami["stage4_result"]

	print("*********** FINAL SUMMARY ***********")
	print("Origami: %s" % (origami_name))
	percent_staples_placed = (staples_placed / origami["num_physical_staples"]) * 100
	print("%d of %d staples placed (%2.1f%%)" % (staples_placed, origami["num_physical_staples"], percent_staples_placed))
	print("%d of %d scaffold bases hybridised\n" % (c1, origami["scaffold_nt"]))

	# write contact map
	cio.write_origami_contactmap(contactmap, "_assets", "%s_naive.csv" % origami_name, madeby="NAIVE")
	print("Reverse engineered origami contact map saved to --> _assets/%s_naive.csv" % (origami_name))
	print("*************************************\n")



