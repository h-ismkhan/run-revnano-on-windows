"""REVNANO : Reverse Engineering Scaffolded Origami Designs from Sequence Information

Ben Shirt-Ediss, 2020-2023
"""

import sys
sys.path.append("../")

import os.path
import time
import pickle

from REVNANO import parameters
from REVNANO import stage0
from REVNANO import stage1
from REVNANO import stage2
from REVNANO import stage3
from REVNANO import stage4
from REVNANO import dna
from REVNANO import iohandler

from contactmap import iohandler as cio
from contactmap import contactutils




#
#
# REVNANO API
#
#

def read_rev_file(origami_name, assets_dir, verbose = False) :

	if verbose == False :
		# all print to screen output is suppressed: https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
		sys.stdout = open(os.devnull, 'w')

	if not assets_dir.endswith("/") :
		assets_dir += "/"		

	filein_name = "%s%s.rev" % (assets_dir, origami_name)

	if not os.path.exists(filein_name) :
		print("\nInitialisation Error: Origami sequences file '%s.rev' does not exist in the _assets/ directory. Stop.\n" % origami_name)
		if verbose == False : sys.stdout = sys.__stdout__
		return None

	try :
		filein = iohandler.load_input_csv(filein_name)	
	except ValueError as errormsg :
		print("\n%s\n" % errormsg)
		if verbose == False : sys.stdout = sys.__stdout__
		return None	

	scaffold_type, scaffold53, staples53raw = filein

	# The origami dictionary holds information about the state of reverse engineering, at all stages of REVNANO
	origami = {} 

	# 1. staples with interior loopouts and dangles (process first)
	# 		(* characters are removed from staple sequences)

	origami["num_physical_staples"] = len(staples53raw)
	
	# .. for recording staples with dangles
	origami["dangles5"] = {} 		# staple id: dangle sequence at 5'
	origami["dangles3"] = {} 		# staple id: dangle sequence at 3'
	
	# .. for recording staples with interior loopouts
	origami["extra_staples"] = {}	# staple id: list of id's of extra (virtual) staples the stub of this staple links to, via loopouts
	origami["loopouts"] = {}		# staple id: list of loopout subsequences

	staples53 = []
	staples53_extra = []

	for staple_id, seq53 in enumerate(staples53raw) :
		dangle5, dangle3, midlist = dna.dangles_loopouts(seq53)

		if dangle5 != None :
			origami["dangles5"][staple_id] = dangle5	# staple has marked dangle at 5' end
		if dangle3 != None :
			origami["dangles3"][staple_id] = dangle3	# staple has marked dangle at 3' end

		staples53.append(midlist[0][1])	

		if len(midlist) > 1 :	# staple has marked interior loopouts, and this staple needs breaking into virtual extra staples
								# SSSSLLLLSSSSS

			origami["extra_staples"][staple_id] = []
			origami["loopouts"][staple_id] = []
			
			for i, item in enumerate(midlist) :
				seqtype, seq = item

				if seqtype == "substaple" :
					if i == 0 :
						continue 	# first stub sequence of staple already added above with staples53.append(midlist[0][1]) 
					else :
						# add substaple as a new independent staple to route, at end of list
						staples53_extra.append(seq)
						origami["extra_staples"][staple_id].append( origami["num_physical_staples"] - 1 + len(staples53_extra) )

				elif seqtype == "loopout" :
					origami["loopouts"][staple_id].append(seq)

	staples53.extend(staples53_extra)

	# tell user which loopout staples where split into virtual staples
	if len(origami["extra_staples"]) > 0 :
		print("Staples with interior loopouts now separated into independent staples:")
		
		for staple_id in origami["extra_staples"] :
			print("%d --> %r" % (staple_id, [staple_id] + origami["extra_staples"][staple_id]))

		print("")
	
	# 2. basic origami attributes

	origami["scaffold_type"] = scaffold_type 				# "LINEAR" or "CIRCULAR"
	origami["scaffold53"] = scaffold53
	origami["scaffold_nt"] = len(scaffold53)
	
	origami["staples53"] = staples53 			# list of staple sequences, index is staple id
	origami["num_staples"] = len(staples53) 	# note: this number also includes virtual staples, created when there are staples with loopouts
	
	# 3. main data structures
	
	# -- routing trees
	origami["staple_routing_trees"] = None

	# -- staple footprints
	origami["footprints"] = {}		# Current staple footprints
									# scaffold base index i53: set() of staple ids claiming this base
									# re-made in stage 3 to be
									# scaffold base index i53: set() of tuples (staple id, section_id)

	for i53 in range(0, origami["scaffold_nt"]) :
		origami["footprints"][i53] = set()

	# -- staples placed
	origami["staples_placed"] = {} 	# Staples placed
									# staple id : list of staple section lists 
									# each staple section list has form [i53 of staple 5 prime, run length, whether this section is HARD SET]
									# note that i53 needs decrementing by run length bases to get to the staple 3'

	for staple_id in range(0, origami["num_staples"]) :
		origami["staples_placed"][staple_id] = []
	
	# 4. staples excluded

	origami["staple_ids_with_1definite_route"] = [] 					# at stage 0
	origami["staple_ids_with_no_tree"] = [] 							# at stage 0
	origami["staple_ids_with_zero_route"] = []							# at stage 1 or 2
	origami["staple_ids_with_persistent_ambiguous_route"] = []			# after stage 2
	origami["staple_ids_with_overlap_problems"] = []					# at stage 3

	origami["G1_staples"] 	 = None

	# 5. stage results

	origami["stage0_result"] = None
	origami["stage1_result"] = None
	origami["stage2_result"] = None
	origami["stage3_result"] = None
	origami["stage4_result"] = None

	if verbose == False : sys.stdout = sys.__stdout__		
		# re-enable print()

	return origami







def reverse_engineer(origami, MU_MIN, SIGMA, BETA, revnano_deterministic = True, verbose = False) :
	"""Run REVNANO. Returns reverse engineered origami contact map, and error message (if any)

	Possible errors are:
		Error #1: [Stage 0] A definite 1-route staple overlaps or fully nests a definite 1-route staple already placed
		Error #2: [Stage 2] Not enough staples placed to use shortest paths to place remaining >1 route staples
		Error #3: [Stage 3] No staples placed, no staple overlaps exist
		Error #4: [Stage 3] Three-staple overlap found
		Error #5: [Stage 3] Some staple-staple overlaps could not be resolved

	Note that the passed origami dictionary is also updated by reference, with results of each stage
	"""

	if verbose == False :
		# all print to screen output is suppressed: https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
		sys.stdout = open(os.devnull, 'w')


	if parameters.SIGMA_1ROUTE < SIGMA :
		errormsg = "Initialisation Error: It must be that parameters.SIGMA_1ROUTE >= SIGMA. Stop."
		if verbose == False : sys.stdout = sys.__stdout__
		return None, errormsg

	try :		
		# note: origami dict is updated by reference in all cases
		stage0.process(origami, MU_MIN, SIGMA, BETA, revnano_deterministic)
		stage1.process(origami, BETA, revnano_deterministic) 			
		stage2.process(origami, BETA, revnano_deterministic)
		stage3.process(origami)
		contactmap = stage4.process(origami)

	except ValueError as errormsg :
		if verbose == False : sys.stdout = sys.__stdout__
		return None, str(errormsg) 	# without string cast, errormsg is a ValueError object

	if verbose == False : sys.stdout = sys.__stdout__		
	return contactmap, None







def make_modified_contact_maps_for_guide_schematic_display(origami_name, assets_dir, verbose = False) :
	""" Assumes original origami contact map exists in assets_dir

		Returns
		0: if contactmap contains no 1nt domains, and thus can already be displayed as a guide schematic
		1: if D and R contact maps are made
		2: if D and R could not be made because a 1nt hybridised domain existed on the origami scaffold
	"""

	if verbose == False :
		# all print to screen output is suppressed: https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
		sys.stdout = open(os.devnull, 'w')

	if not assets_dir.endswith("/") :
		assets_dir += "/"

	# load original contact map
	contactmap = cio.read_origami_contactmap("%s%s.csv" % (assets_dir, origami_name))		

	I = contactutils.scaffold_1nt_domain_locations(contactmap)
	if len(I) > 0 :
		print("!! WARNING: GUIDE SCHEMATIC DISPLAY ISSUES !!")
		print("--> In the reverse engineered origami contact map, scaffold bases %r are 1nt scaffold domains (where the scaffold 5' has base index 0)" % I)
		print("--> With 1nt domains, the origami contact map cannot be converted to a domain-level graph and hence it cannot be displayed as an origami guide schematic")
		print("--> Two alternative APPROXIMATE contact maps will be made for this origami which will allow it to be displayed as an APPROXIMATE guide schematic:")

		try :
			print("1) A contact map with some staples removed, to enlarge the 1nt domains:")

			contactmapR, staple_ids_removed = contactutils.remove_staples_to_remove_1nt_scaffold_domains(contactmap)
			# throws ValueError if a 1nt scaffold domain hybridised to a staple is found: this case should not happen
			# contactmapR: junction ambiguity is accurate, but large parts of the nanostructure may not display well if whole staples missing

			print("--> Staples %r removed" % staple_ids_removed)
			cio.write_origami_contactmap(contactmapR, assets_dir, "%sR.csv" % origami_name, madeby="REVNANO (staples removed to remedy 1nt domains)")
			print("--> Contact map saved to --> %s%sR.csv" % (assets_dir, origami_name))

			print("2) A contact map with all 1nt domains deleted from the scaffold (shorter scaffold):")

			contactmapD = contactutils.delete_1nt_scaffold_domains(contactmap)
			# contactmapD: junction ambiguity may not be correct around 1nt domain deletion sites, but all staples in contact map will be present

			print("--> Scaffold bases %r deleted from origami" % I)
			cio.write_origami_contactmap(contactmapD, assets_dir, "%sD.csv" % origami_name, madeby="REVNANO (scaffold 1nt domains deleted)")
			print("--> Contact map saved to --> %s%sD.csv" % (assets_dir, origami_name))	

		except ValueError as exp :
			print("Failed. A 1nt scaffold domain hybridised to a staple was found.")
			if verbose == False : sys.stdout = sys.__stdout__
			return 2

		if verbose == False : sys.stdout = sys.__stdout__
		return 1
	else :
		if verbose == False : sys.stdout = sys.__stdout__
		return 0	






































#
#
# REVNANO CLI
#
#

if __name__ == "__main__" :
	"""Usage from command line:

	revnano.py origami_name

	Note: running from the command line uses the default values
	of SIGMA and BETA in parameters.py. These may not be optimal for all origamis.
	"""

	print("\n")
	print("-----------------------------------------------------------")
	print("REVNANO Constraint Programming Solver v%1.1f" % parameters.VERSION)
	print("Origami Staple/Scaffold Sequences ----> Origami Contact Map\n")
	print("Please cite our 2023 paper:")
	print("%s" % parameters.PAPER_DOI)
	print("-----------------------------------------------------------")

	# --------------------
	# Check for correct number of arguments
	# --------------------

	MU_MIN = parameters.MU_MIN
	SIGMA = parameters.SIGMA
	BETA = parameters.BETA

	if len(sys.argv) not in [2,5] :
		print("\nUsage: revnano.py origami_name")	# uses default parameters in parameters.py
		print("or: revnano.py origami_name MU_MIN SIGMA BETA\n")
		quit()

	if len(sys.argv) == 5 :
		try :
			MU_MIN = int(sys.argv[2])
			SIGMA = int(sys.argv[3])
			BETA = float(sys.argv[4])
		except ValueError :
			print("\nError: One of the supplied parameters MU_MIN SIGMA BETA is the wrong data type.\n")
			quit()

		if MU_MIN not in [5,6] :
			print("\nError: Parameter MU_MIN must be 5 or 6.\n")
			quit()
		if not ((SIGMA >= 0) and (SIGMA <= 10)) :
			print("\nError: Parameter SIGMA must be between 0 and 10.\n")
			quit()
		if not ((BETA > 0) and (BETA <= 1)) :
			print("\nError: Parameter BETA must be greater than 0 and less than 1.\n")
			quit()

	origami_name = str(sys.argv[1]).replace(".rev", "")

	print("\nOrigami to reverse engineer:")
	print("%s\n" % origami_name)
	
	print("Algorithm main parameters:")
	print("MU_MIN:\t\t%d bp" % MU_MIN)
	print("SIGMA:\t\t%d bp" % SIGMA)
	print("BETA:\t\t%1.2f\n" % BETA)

	print("Deterministic staple placement?:")
	deterministic = "YES - staples polled in order of ascending staple id"
	if parameters.IS_DETERMINISTIC == False :
		deterministic = "NO - staples polled in random order"
	print("%s\n" % deterministic)

	# --------------------
	# Read raw sequences file, and initialise origami data object to hold sequences, and other information
	# --------------------

	origami = read_rev_file(origami_name, "_assets/", verbose = True)

	if origami == None :
		quit()
	
	# --------------------
	# Reverse engineer origami contact map from raw sequences
	# --------------------

	t = time.time()

	# reverse engineer origami contact map  (note: passed origami dict is also updated by reference)
	contactmap, errormsg = reverse_engineer(origami, MU_MIN, SIGMA, BETA, parameters.IS_DETERMINISTIC, verbose = True)

	elapsed_sec = time.time() - t

	# check if an error occurred that prevented a result being returned
	if errormsg != None :
		print("\n%s\n" % errormsg)
		quit()

	# --------------------
	# Results summary
	# --------------------

	c1, staples_placed = origami["stage4_result"]

	print("*********** FINAL SUMMARY ***********")

	# details of problematic staples (before loopout staples re-composed)
	if staples_placed < origami["num_physical_staples"] :

		print("Staples omitted:")

		if len(origami["staple_ids_with_no_tree"]) > 0 :
			print("Staples for which no viable scaffold routes were found at Stage 0")
			for staple_id in origami["staple_ids_with_no_tree"] :
				print("%d\t5'--%s--3'" % (staple_id, origami["staples53"][staple_id]))

		if len(origami["staple_ids_with_zero_route"]) > 0 :
			print("Staples which collapsed to 0 viable scaffold routes, due to the placement of other staples at Stage 1 or 2")
			for staple_id in origami["staple_ids_with_zero_route"] :
				print("%d\t5'--%s--3'" % (staple_id, origami["staples53"][staple_id]))			

		if len(origami["staple_ids_with_persistent_ambiguous_route"]) > 0 :
			print("Staples omitted because their route remained ambiguous after Stage 2")
			for staple_id in origami["staple_ids_with_persistent_ambiguous_route"] :
				print("%d\t5'--%s--3'" % (staple_id, origami["staples53"][staple_id]))		

		if len(origami["staple_ids_with_overlap_problems"]) > 0 :
			print("Staples omitted because their overlaps with other staples could not be resolved at Stage 3")
			for staple_id in origami["staple_ids_with_overlap_problems"] :
				print("%d\t5'--%s--3'" % (staple_id, origami["staples53"][staple_id]))

		print("")	

	print("Origami: %s" % (origami_name))
	percent_staples_placed = (staples_placed / origami["num_physical_staples"]) * 100
	print("%d of %d staples placed (%2.1f%%)" % (staples_placed, origami["num_physical_staples"], percent_staples_placed))
	print("%d of %d scaffold bases hybridised\n" % (c1, origami["scaffold_nt"]))

	# save the reverse engineered contact map as CSV
	cio.write_origami_contactmap(contactmap, "_assets", "%s.csv" % origami_name, madeby="REVNANO")
	print("Reverse engineered origami contact map saved to --> _assets/%s.csv" % (origami_name))
	
	print("*************************************\n")

	# --------------------
	# Make separate display contact maps, if the derived contact map has 1nt domains
	# --------------------

	make_modified_contact_maps_for_guide_schematic_display(origami_name, "_assets/", verbose = True)

	# --------------------
	# save all reverse engineering data (for future programatic analysis)
	# --------------------

	pickle.dump( (origami, contactmap), open( "_assets/%s.data" % origami_name, "wb" ) )

	print("\n[REVNANO running time %.2f seconds]\n" % elapsed_sec)



