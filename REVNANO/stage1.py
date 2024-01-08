"""REVNANO : Stage 1 : Find approximate positions for staples by using an iterative tree pruning approach 

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
# PROPAGATE CONSTRAINTS ("SUDOKU STYLE")
#
#

def place_staples_by_constraint_propagation(origami, BETA, revnano_deterministic) :
	""" "Crystalises" the footprint of each staple, by pruning the tree of each staple
	in response to the current footprints of other staples. Iterative process.

	Going into this function, only the footprints of the definite 1-route staples exist
	All other staples exist just as routing trees which are yet to cast a footprint

	Loop until there are no footprint changes:
		For each staple: 	(staples polled in order, or by random)
			- Get staple
			- Update its routing tree, according to current state of footprints of other staples (origami["footprints"])
			- Cast an enlarged footprint on the scaffold for this staple, updating origami["footprints"]
			- Move to next staple

		(by the time end of staples reached, the pruning state of SRTs of
		staples polled first will be **out of date**: but the whole batch of staples is gone
		through again on the next iteration)

	Afterwards, those staples reducing to 1 route are added to the St dict
	"""

	# iterative loop, where staple trees successively prune each other
	claimed_diff = 1
	iteration = 1
	c0 = 0
	c1 = 0

	c0 = dna.count_scaffold_bases_claimed(origami)
	print("%d bases already claimed by definite 1-route staples\n" % c0, flush=True)

	print("Crystallising staple footprints...", flush=True)

	while claimed_diff != 0 :

		print("--> Iteration %d" % iteration)

		# decide staple poll order for this iteration
		poll_order_staple_ids = list(range(0, origami["num_staples"]))
		if revnano_deterministic == False :
			# poll in random order: algorithm no longer strictly deterministic 
			random.shuffle(poll_order_staple_ids)
			
		for staple_id in poll_order_staple_ids :

			if staple_id in origami["staple_ids_with_no_tree"] :
				continue 	# skip if staple had no tree built at Stage 0

			if staple_id in origami["staple_ids_with_zero_route"] :
				continue 	# skip if staple has already collapsed to 0 viable routes at Stage 0

			if staple_id in origami["staple_ids_with_1definite_route"] :
				continue	# skip if staple is a deinite 1-route staple, and has been already placed at Stage 0

			# A. get current routing tree for staple
			rtree = origami["staple_routing_trees"][staple_id]

			# B. prune any new nodes which have become blocked by other staples (uses origami["footprints"])
			srt.prune_blocked_nodes_of_staple_routing_tree(origami, staple_id, rtree, BETA)

			# C. propagate the effect of pruned nodes up and down the tree
			#	 and so only valid routing paths remain
			root_pruned = srt.propagate_pruning_in_staple_routing_tree(rtree)
				# if the root node is pruned, then no viable routing paths exist for this staple on the scaffold now
		
			if root_pruned :
				print("\tStaple %d has collapsed to 0 valid routes" % (staple_id))

				# remove any existing footprint this staple had on the scaffold
				dna.remove_staple_footprint_from_scaffold(origami, staple_id)

				origami["staple_ids_with_zero_route"].append(staple_id)
			else :
				# D. routing tree is still valid: cast enlarged/same sized staple footprint on scaffold
				footprint = [0] * origami["scaffold_nt"]
				num_staple_routes = srt.staple_routing_tree_footprint(origami, rtree, footprint)
				dna.cast_staple_footprint_on_scaffold(origami, staple_id, footprint, num_staple_routes)
					# updates origami["footprints"]

		c1 = dna.count_scaffold_bases_claimed(origami)
		claimed_diff = c1 - c0
		print("\t%d scaffold bases now claimed (+%d since prev)" % (c1, claimed_diff))
		c0 = c1

		iteration += 1

	print("[Done]\n", flush=True)















#
#
# EXECUTE
#
#

def process(origami, BETA, revnano_deterministic) :

	print("-----------------------------------------------------------\n")
	print("REVNANO STAGE 1 (Place staples by propagating constraints) \n")
	print("-----------------------------------------------------------\n")

	# --------------------
	# Propagate constraints and place staples
	# --------------------

	place_staples_by_constraint_propagation(origami, BETA, revnano_deterministic)

	# --------------------
	# Record positions of newly placed single route staples in St dictionary
	# --------------------

	for staple_id, staple_seq in enumerate(origami["staples53"]) :	

		origami["staples_placed"][staple_id] = []
		rtree = origami["staple_routing_trees"][staple_id]
		
		if rtree == None : 	# skip staples that never had a staple routing tree made in the initial stage0
			continue

		leaf_histories = srt.list_of_valid_leaf_histories(rtree, [])

		if len(leaf_histories) == 1 : # single unique route
			for section in leaf_histories[0] :
				i53, runlen = section
				origami["staples_placed"][staple_id].append( [i53, runlen] )

	# --------------------
	# Stats
	# --------------------

	countNoTree, count0, count1, countG1, G1av, G1max = srt.routing_tree_stats(origami)
	c1 = dna.count_scaffold_bases_claimed(origami)
	staples_placed = dna.count_staples_placed(origami)
	origami["stage1_result"] = (countNoTree, count0, count1, countG1, G1av, G1max, c1, staples_placed)

	str_no_tree = ""
	if countNoTree > 0 :
		str_no_tree = "(Staples: %r)" % origami["staple_ids_with_no_tree"]

	str_zero_routes = ""
	if count0 > 0 :
		str_zero_routes = "(Staples: %r)" % origami["staple_ids_with_zero_route"]

	print("The %d staples have status:" % origami["num_staples"])
	print("- No tree: \t%d\t%s" % (countNoTree, str_no_tree))
	print("- 0 routes: \t%d\t%s" % (count0, str_zero_routes))
	print("- 1 route: \t%d" % count1)
	print("- >1 routes: \t%d" % countG1)
	print("\tAverage number of routes on >1 route staples: \t%1.2f" % (G1av))
	print("\tMaximum number of routes on >1 route staples: \t%d" % (G1max))

	print("")

	if count1 == origami["num_staples"] :
		print("\n ** All staples placed on origami ** \n")
	else :
		print("%d staples currently placed on origami\n" % (staples_placed))
