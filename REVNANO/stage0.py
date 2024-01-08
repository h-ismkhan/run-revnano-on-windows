"""REVNANO : Stage 0 : Read origami sequence file and compute a routing tree for each staple

Ben Shirt-Ediss, 2020, 2023
"""

import sys
sys.path.append("../")

import random

from REVNANO import parameters
from REVNANO import srt
from REVNANO import dna




#
#
# CREATE STAPLE ROUTING TREES
#
#


"""
Structure of a staple routing tree data structure:

A group of child nodes, coming off a parent node, is defined like this

tree[(i53, runlen)] = [ history_CHILD, runlen_so_far + runlen, False, subtree ]
							1 				2 					3 		4

1 = the PATH to get to this node, a list [(i53, runlen),(i53, runlen)]
2 = bases on staple placed so far
3 = true if node is pruned (i.e. forms part of a path which cannot reach the total staple length, or is a node hybridising to an already claimed part of the scaffold)
4 = the CHILDREN of this node, a tree dictionary as above
"""


def build_staple_routing_trees(origami, MU_MIN, SIGMA, BETA, revnano_deterministic) :
	""" Two-stage procedure:

	Pass 1: find all definite 1-route staples (for the value of parameters.SIGMA_1ROUTE)
	Construct trees of definite 1-route staples, then place these staples (in footprints dict, and St dict)

	This optimisation step effectivley reduces the length of the scaffold for routing trees of other
	staples to be created in Pass 2.

	Pass 2: construct trees of all other staples (for the passed value of SIGMA =< parameters.SIGMA_1ROUTE) 
	with the constraint that definite 1-route staples cannot be overlapped (subject to tolerance parameter BETA)

	Updates origami["staple_routing_trees"] variable to be a dictionary     
		staple_id: routing_tree (a dictionary of children, where each child points to a dict of its own children)
	"""

	staple_routing_trees = {}

	# --------------------
	# Pass 1
	# --------------------
	
	print("Finding definite 1-route staples... [SIGMA_1ROUTE = %d, MU_MIN = %d]" % (parameters.SIGMA_1ROUTE, MU_MIN), flush=True)
	
	# decide poll order
	poll_order_staple_ids = list(range(0, origami["num_staples"]))
	if revnano_deterministic == False :
		# poll in random order: algorithm no longer strictly deterministic 
		random.shuffle(poll_order_staple_ids)


	for staple_id in poll_order_staple_ids :

		staple_seq = origami["staples53"][staple_id]

		if srt.staple_has_1route_fast_check(origami, staple_id, staple_seq, set(), 0, len(staple_seq), parameters.SIGMA_1ROUTE, MU_MIN, BETA) :
			# yes, it is a definite 1-route staple
			
			# ------ gen routing tree ------
			srt.subtree_count = 0
			try :
				rtree = srt.generate_staple_routing_tree(origami, staple_id, staple_seq, [], set(), 0, len(staple_seq), parameters.SIGMA_1ROUTE, MU_MIN, BETA)
			except ValueError as exp :
				# staple routing tree too large
				staple_routing_trees[staple_id] = None
				origami["staple_ids_with_no_tree"].append(staple_id)
				print("\tFailed. Staple routing tree growing too large (> %d subtrees)" % parameters.MAX_SUBTREES)
				continue
			# -------------------------------

			staple_routing_trees[staple_id] = rtree

			# check placing this staple does not overlap other definite 1-route staples
			# these definite 1-route staples are NOT allowed to clash: if they do, its a terminal error

			srt.prune_blocked_nodes_of_staple_routing_tree(origami, staple_id, rtree, BETA)
			root_pruned = srt.propagate_pruning_in_staple_routing_tree(rtree)

			if rtree == None :
				root_pruned = True 	# if no tree at all, correct the root pruned flag

			if root_pruned :
				errormsg = "REVNANO Error #1 [Stage 0]: Definite 1-route staple %d overlaps or fully nests a definite 1-route staple already placed. Stop.\n" % staple_id
				raise ValueError(errormsg)
			else :
				# routing tree is valid:
				# UPDATE footprints dict with definite 1-route staple
				footprint = [0] * origami["scaffold_nt"]
				num_staple_routes = srt.staple_routing_tree_footprint(origami, rtree, footprint)
				dna.cast_staple_footprint_on_scaffold(origami, staple_id, footprint, num_staple_routes)
					# updates origami["footprints"]

				# UPDATE St dict with definite 1-route staple
				leaf_histories = srt.list_of_valid_leaf_histories(rtree, [])
				for section in leaf_histories[0] :
					i53, runlen = section
					origami["staples_placed"][staple_id].append( [i53, runlen] )

				origami["staple_ids_with_1definite_route"].append(staple_id)
			
			print("--> 1-route: staple %d" % staple_id, flush=True)

	print("[Done]\n", flush=True)

	# --------------------
	# Pass 2
	# --------------------

	print("Generating routing trees for other staples... [SIGMA = %d, MU_MIN = %d]" % (SIGMA, MU_MIN))

	for staple_id, staple_seq in enumerate(origami["staples53"]) :

		if staple_id in origami["staple_ids_with_1definite_route"] :		# already has tree made
			continue

		print("--> Making routing tree for staple %d" % staple_id, flush=True)

		if len(staple_seq) <= parameters.SHORT_STAPLE_LEN :
			# SHORT STAPLE
			print("\tShort staple of %d bases. Using MU_MIN = %d" % (len(staple_seq), MU_MIN-1))

			# ------ gen routing tree ------
			srt.subtree_count = 0
			try :
				rtree = srt.generate_staple_routing_tree(origami, staple_id, staple_seq, [], set(), 0, len(staple_seq), SIGMA, MU_MIN-1, BETA)
			except ValueError as exp :
				# staple routing tree too large
				staple_routing_trees[staple_id] = None
				origami["staple_ids_with_no_tree"].append(staple_id)
				print("\tFailed. Staple routing tree growing too large (> %d subtrees)" % parameters.MAX_SUBTREES)
				continue			
			# -------------------------------

			root_pruned = srt.propagate_pruning_in_staple_routing_tree(rtree)
			
			if rtree == None :
				root_pruned = True 	# if no tree at all, correct the root pruned flag

			if root_pruned :
				staple_routing_trees[staple_id] = None
				origami["staple_ids_with_no_tree"].append(staple_id)
				print("\tFailed. Staple cannot be routed")
			else :
				# root not pruned, valid possible routings exist
				staple_routing_trees[staple_id] = rtree
				print("\t[Success]")

		else :
			# NORMAL LENGTH STAPLE
			
			# ------ gen routing tree ------
			srt.subtree_count = 0
			try :
				rtree = srt.generate_staple_routing_tree(origami, staple_id, staple_seq, [], set(), 0, len(staple_seq), SIGMA, MU_MIN, BETA)
			except ValueError as exp :
				# staple routing tree too large
				staple_routing_trees[staple_id] = None
				origami["staple_ids_with_no_tree"].append(staple_id)
				print("\tFailed. Staple routing tree growing too large (> %d subtrees)" % parameters.MAX_SUBTREES)
				continue
			# ------------------------------

			root_pruned = srt.propagate_pruning_in_staple_routing_tree(rtree)
			
			if rtree == None :
				root_pruned = True 	# if no tree at all, correct the root pruned flag

			# staple may have very short sections, thats why graph cannot be made
			if root_pruned : 
				print("\tFailed, but now retrying at MU_MIN = %d" % (MU_MIN-1))
				
				# ------ gen routing tree ------
				srt.subtree_count = 0
				try :
					rtree = srt.generate_staple_routing_tree(origami, staple_id, staple_seq, [], set(), 0, len(staple_seq), SIGMA, MU_MIN-1, BETA)			
				except ValueError as exp :
					# staple routing tree too large
					staple_routing_trees[staple_id] = None
					origami["staple_ids_with_no_tree"].append(staple_id)
					print("\tFailed. Staple routing tree growing too large (> %d subtrees)" % parameters.MAX_SUBTREES)
					continue					
				# ------------------------------

				root_pruned = srt.propagate_pruning_in_staple_routing_tree(rtree)

				if rtree == None :
					root_pruned = True 	# if no tree at all, correct the root pruned flag

				if root_pruned :
					staple_routing_trees[staple_id] = None
					origami["staple_ids_with_no_tree"].append(staple_id)
					print("\tFailed. Staple cannot be routed")
				else :
					staple_routing_trees[staple_id] = rtree
					print("\t[Success]")
			else :
				# root not pruned, valid routings exist
				staple_routing_trees[staple_id] = rtree

	print("[Done]\n", flush=True)				

	origami["staple_routing_trees"] = staple_routing_trees













#
#
# EXECUTE
#
#

def process(origami, MU_MIN, SIGMA, BETA, revnano_deterministic) :

	print("-----------------------------------------------------------\n")
	print("REVNANO STAGE 0 (Build staple routing trees)               \n")
	print("-----------------------------------------------------------\n")

	# --------------------
	# Build staple routing trees
	# --------------------

	build_staple_routing_trees(origami, MU_MIN, SIGMA, BETA, revnano_deterministic)		
		# updates origami["staple_routing_trees"] and origami["footprints"] with definite 1-route staples

	# --------------------
	# Stats
	# --------------------

	countNoTree, count0, count1, countG1, G1av, G1max = srt.routing_tree_stats(origami)
	origami["stage0_result"] = (countNoTree, count0, count1, countG1, G1av, G1max, 0, 0)

	str_no_tree = ""
	if countNoTree > 0 :
		str_no_tree = "(Staples: %r)" % origami["staple_ids_with_no_tree"]

	print("The %d staples have status:" % origami["num_staples"])
	print("- No tree: \t%d\t%s" % (countNoTree, str_no_tree))
	print("- 1 route: \t%d" % count1)
	print("- >1 routes: \t%d" % countG1)
	print("\tAverage number of routes on >1 route staples: \t%1.2f" % (G1av))
	print("\tMaximum number of routes on >1 route staples: \t%d" % (G1max))

	print("")	

