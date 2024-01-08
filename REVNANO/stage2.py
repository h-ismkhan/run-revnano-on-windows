"""REVNANO : Stage 2 : Use shortest path approach to solve problem >1 route staples, after origami is largely solved by Stage 1

(This stage seems necessary particularly for the M13mp18 scaffold, which has many repeats)

Ben Shirt-Ediss, 2020, 2023
"""

import sys
sys.path.append("../")

import networkx as nx 		# for shortest path calculation
import random

from REVNANO import parameters
from REVNANO import srt
from REVNANO import dna





#
#
# ORIGAMI BASE GRAPH
#
#

def initialise_base_graph_B(origami) :
	"""Goes through each staple in St, and creates an approximate base connectivity graph for the origami

	The graph B is a networkX Graph object

	Note: Because the solver de-constructs staples with loop-outs before starting, staple loopoouts
	do not contribute to shortest path distances. However, this omission is expected to be minor.
	"""
	
	St = origami["staples_placed"]
	B = {}

	# 1. make graph of just the scaffold
	for i53 in range(0, origami["scaffold_nt"]) :

		B[i53] = {}	# dictionary for neighbours of this base

		# join to previous base
		if i53 > 0 :
			B[i53][i53-1] = 1
		else :
			# i53 = base at 5 prime end
			if origami["scaffold_type"] == "CIRCULAR" :
				B[i53][origami["scaffold_nt"]-1] = 1 # join 5 prime backwards, to meet 3 prime end and complete loop

		# join to next base
		if i53 < (origami["scaffold_nt"] - 1) :
			B[i53][i53+1] = 1
		else :
			# i53 = base at 3 prime base
			if origami["scaffold_type"] == "CIRCULAR" :
				B[i53][0] = 1  # join 3 prime forwards, to meet 5 prime end and complete loop


	# 2. add staple crossovers, to make shorter range base connections
	for staple_id in St :

		crossover_pair_list = dna.staple_crossover_pairs(origami, staple_id)

		for pair in crossover_pair_list :

			i53_section_3prime, i53_next_section_5prime = pair

			# update B graph - two way link
			B[i53_section_3prime][i53_next_section_5prime] = 1
			B[i53_next_section_5prime][i53_section_3prime] = 1			

	origami["B"] = nx.Graph(B)  # save graph as a NetworkX graph, for shortest path querying




def staple_routing_distance(origami, route) :
	"""Returns the total distance that a staple route has
	when all of its crossovers are routed via the existing origami mesh (and not by the staple itself)
	"""

	route_len = 0
	i53_section_3prime = None

	for section_id, section in enumerate(route) :

		i53_section_5prime, runlen = section

		# with crossover from previous section -- find how long is the route between crossover bases, via the existing origami mesh
		if section_id > 0 :
			sp = nx.shortest_path(origami["B"],i53_section_3prime,i53_section_5prime)
			route_len += len(sp)

		# calc 3 prime base of this section
		i53_section_3prime = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53_section_5prime - runlen + 1)

		# add this dsDNA section length to route length
		route_len += runlen

	return route_len























#
#
# G1 STAPLES METHODS
#
#

def initialise_G1_staples_dictionary(origami) :
	"""Makes origami["G1_staples"], a dictionary of staples that still have more than 1 valid route

	G1_staples dictionary has the form:
	staple_id : [num_routes, routes_list, list of shortest path length for each route, true/false shortest path is unique, true/false added to origami]
	"""

	G1_staples = {}

	for staple_id, _ in enumerate(origami["staples53"]) :

		rtree = origami["staple_routing_trees"][staple_id]
		if rtree == None : # skip staples that never had a staple routing tree made in the initial stage0
			continue

		routes_list = srt.list_of_valid_leaf_histories(rtree, [])
		num_routes = len(routes_list)

		if(num_routes > 1) :			
			G1_staples[staple_id] = [0, [], [], False]  # place holder

	origami["G1_staples"] = G1_staples





def update_G1_staples_dictionary(origami) :
	"""Updates the G1 set of staples
		Deletes staples which now have 0 or 1 routes
	"""

	del_list = []

	for staple_id in origami["G1_staples"]:

		rtree = origami["staple_routing_trees"][staple_id]
		if rtree == None :
			continue

		routes_list = srt.list_of_valid_leaf_histories(rtree, [])
		num_routes = len(routes_list)

		if num_routes == 0 or num_routes == 1 :
			del_list.append(staple_id)

			if num_routes == 0 :
				print("\tStaple %d has collapsed to 0 valid routes" % staple_id)
				origami["staple_ids_with_zero_route"].append(staple_id)
		else : 
			# num_routes > 1
			# find length of all routes
			route_len_list = []
			for route in routes_list :
				route_len_list.append( staple_routing_distance(origami, route) )

			# see if shortest route is well defined
			shortestpath_unique = True
			route_len_list_sorted = sorted(route_len_list)
			if route_len_list_sorted[0] == route_len_list_sorted[1] :
				shortestpath_unique = False

			origami["G1_staples"][staple_id] = [num_routes, routes_list, route_len_list, shortestpath_unique]

	for staple_id in del_list :
		del origami["G1_staples"][staple_id]





def add_1route_G1_staples_to_origami(origami) :
	"""Adds to the origami staples in the G1 set that have a single route (and returns number of single route staples)
		This modifies the shortest path routes through the origami

		These single route staples are then cleared out of G1 the next time it is updated
	"""

	num1 = 0

	for staple_id in origami["G1_staples"] :

		rtree = origami["staple_routing_trees"][staple_id]
		routes_list = srt.list_of_valid_leaf_histories(rtree, [])
		num_routes = len(routes_list)
		
		if num_routes == 1 :   # 1 route staple, not yet added to origami

			# update St -- staples dictionary
			for section in routes_list[0] :
				i53, runlen = section
				origami["staples_placed"][staple_id].append( [i53, runlen] )

			# update B -- base connectivity graph, to include edges for crossovers of this staple
			crossover_pair_list = dna.staple_crossover_pairs(origami, staple_id)
			for pair in crossover_pair_list :
				i53_section_3prime, i53_next_section_5prime = pair
				origami["B"].add_edge(i53_section_3prime, i53_next_section_5prime)

			num1 += 1

	return num1




def select_G1_staple_with_clearest_shortest_route(origami) :
	"""Finds the G1 staple whose shortest path is the longest distance away from the next shortest path

	Returns the id of the staple placed
	Returns -1 if no staple placed, i.e. if all staples in G1_staples are "0", "1" or "G1 with multiple shortest paths" staples

	Casts the footprint of this staple on the scaffold

	The staple with the clearest shortest route is not correct in small % of cases

	Note: If 2 or more staples share the clearest shortest path, the one with the higher staple_id is selected
	Note: Staples with non-unique shortest paths can only be resolved if another staple is set, with blocks one of these shortest paths
	"""

	# 1. find id of staple with the clearest shortest path
	clearest_id = -1
	sp = []
	gdiff = 0   # greatest diff between shortest path and next shortest path, in all G1 staple set

	for staple_id in origami["G1_staples"] :

		num_routes, routes_list, route_len_list, shortestpath_unique = origami["G1_staples"][staple_id]

		# assess this staple, if its a G1 staple with a unique shortest path
		if num_routes > 1 and shortestpath_unique :

			route_len_list_sorted = sorted(route_len_list)
				# in ascending order: index 0=shortest path, 1=second shortest path

			diff = route_len_list_sorted[1]-route_len_list_sorted[0]
			if diff > gdiff :
				clearest_id = staple_id
				sp = routes_list[route_len_list.index(route_len_list_sorted[0])]
				gdiff = diff
	
	if clearest_id == -1 :
		return -1

	# prune all routes in routing tree, except the shortest route
	rtree = origami["staple_routing_trees"][clearest_id]
	srt.prune_all_but_shortest_route_of_staple_routing_tree(rtree, sp)
	srt.propagate_pruning_in_staple_routing_tree(rtree) 		# root cannot be pruned in this operation, as a single path will now exist

	# cast the footprint of the staple on the scaffold
	footprint = [0] * origami["scaffold_nt"]
	srt.staple_routing_tree_footprint(origami, rtree, footprint)
	dna.cast_staple_footprint_on_scaffold(origami, clearest_id, footprint, 1)

	return clearest_id




def place_G1_staples_by_constraint_propagation(origami, BETA, revnano_deterministic) :
	"""Updates all staples in the G1 set, once one of them is forced to a single route via shortest path
	Accounts for all the "ripple" knock-on effects. Once these are dealt with, the next G1 staple can be forced to a single route by shortest path.
	"""

	# identical to Stage 1 loop
	claimed_diff = 1
	c0 = dna.count_scaffold_bases_claimed(origami)
	c1 = 0

	while claimed_diff != 0 :

		# polling deterministic, or non-deterministic
		poll_order_staple_ids = list(origami["G1_staples"].keys())	# in ascending order
		if revnano_deterministic == False :
			random.shuffle(poll_order_staple_ids)

		for staple_id in poll_order_staple_ids :

			# A. get current routing tree for G1 staple
			rtree = origami["staple_routing_trees"][staple_id]

			if rtree == None :
				continue

			# B. prune any new nodes which have become blocked by other staples
			srt.prune_blocked_nodes_of_staple_routing_tree(origami, staple_id, rtree, BETA)

			# C. propagate the effect of pruned nodes up and down the tree
			#	 and so only valid routing paths remain
			root_pruned = srt.propagate_pruning_in_staple_routing_tree(rtree)
		
			if root_pruned :
				# 0 route staple has its footprint removed from scaffold
				dna.remove_staple_footprint_from_scaffold(origami, staple_id)
			else :
				# D. cast enlarged/same sized staple footprint on scaffold
				footprint = [0] * origami["scaffold_nt"]
				num_staple_routes = srt.staple_routing_tree_footprint(origami, rtree, footprint)
				dna.cast_staple_footprint_on_scaffold(origami, staple_id, footprint, num_staple_routes)

		c1 = dna.count_scaffold_bases_claimed(origami)
		claimed_diff = c1 - c0
		c0 = c1




















#
#
# EXECUTE
#
#

def process(origami, BETA, revnano_deterministic) :
	
	print("-----------------------------------------------------------\n")
	print("REVNANO STAGE 2 (Place remaining staples by shortest path) \n")
	print("-----------------------------------------------------------\n")

	# --------------------
	# Skip this stage if all staples already placed in Stage 1
	# --------------------

	_, _, count1, countG1, _, _, _, _ = origami["stage1_result"]

	if countG1 == 0 :
		print("No staples of >1 route exist to be placed.")
		print("(Stage 2 skipped)\n")
		return

	# --------------------
	# Exit with error if not enough staples are placed in Stage 1
	# --------------------

	fraction_staples_placed = count1 / origami["num_staples"]

	if fraction_staples_placed < parameters.GAMMA :
		# unrealiable to calculate shortest paths with so little staples placed
		errormsg = "REVNANO Error #2 [Stage 2]: Only %2.1f%% of staples have a single route on the origami (below the GAMMA limit). Stop.\n" % (fraction_staples_placed*100)
		raise ValueError(errormsg)

	# --------------------
	# Build origami base connectivity graph
	# --------------------

	initialise_G1_staples_dictionary(origami)
	initialise_base_graph_B(origami)

	# --------------------
	# Place staples by shortest path
	# --------------------

	print("%d staples have >1 route and remain to be placed" % len(origami["G1_staples"]))

	print("\nForce-placing remaining staples by shortest path heuristic...", flush=True)
	
	iteration = 1 		# number of staples dijkstras algorithm forced placement of
	num1_total = 0		# number of staples reduced to 1 route

	# 2. reduce G1 list iteratively
	clearest_id = 0
	while clearest_id > -1 :

		# update the status of all staples in the original G1 set (number of routes, shortest path lengths)
		# clear out 0 and 1 route staples
		update_G1_staples_dictionary(origami)

		# select the G1 staple with the clearest shortest path,
		# prune its routing tree to the shortest path, and place its footprint on the scaffold
		clearest_id = select_G1_staple_with_clearest_shortest_route(origami)
		if clearest_id == -1 :  
			continue			# no G1 clearest path staple exists, exit while loop

		print("--> Iteration %d" % iteration)
		print("\tStaple %d forced to its shortest path route" % clearest_id)
		
		# update G1 staple set iteratively, until the ripple effects of force placing
		# the staple above come to an end
		place_G1_staples_by_constraint_propagation(origami, BETA, revnano_deterministic)

		# when done, add new 1-route staples in G1 staple set to the origami
		num1 = add_1route_G1_staples_to_origami(origami)
		print("\tRipple effect: %d staples in total now become single route" % (num1))
		num1_total += num1

		iteration += 1


	# take footprints of still existing G1 staples off the origami scaffold (they are never added to St)
	# these staples are OMITTED from the final design
	for staple_id in origami["G1_staples"] :
		origami["staple_ids_with_persistent_ambiguous_route"].append(staple_id)
		dna.remove_staple_footprint_from_scaffold(origami, staple_id)
	
	print("[Done]\n", flush=True)

	print("%d staples reduced to 1 route\n" % (num1_total))

	# --------------------
	# Stats
	# --------------------

	countNoTree, count0, count1, countG1, G1av, G1max = srt.routing_tree_stats(origami)
	c1 = dna.count_scaffold_bases_claimed(origami)
	staples_placed = dna.count_staples_placed(origami)
	origami["stage2_result"] = (countNoTree, count0, count1, countG1, G1av, G1max, c1, staples_placed)
	
	# --- error if force placement of staples has caused too many staples to collapse to 0 routes --
	# 		(ensures that number of staples placed is always above gamma threshold)
	fraction_staples_placed = count1 / origami["num_staples"]
	if fraction_staples_placed < parameters.GAMMA :
		errormsg = "REVNANO Error #2 [Stage 2]: Only %2.1f%% of staples have a single route on the origami (below the GAMMA limit). Stop.\n" % (fraction_staples_placed*100)
		raise ValueError(errormsg)
	# ----------

	str_no_tree = ""
	if countNoTree > 0 :
		str_no_tree = "(Staples: %r)" % origami["staple_ids_with_no_tree"]

	str_zero_routes = ""
	if count0 > 0 :
		str_zero_routes = "(Staples: %r)" % origami["staple_ids_with_zero_route"]

	str_g1_routes = ""
	if countG1 > 0 :
		str_g1_routes = "(Staples: %r)" % origami["staple_ids_with_persistent_ambiguous_route"]

	print("The %d staples have status:" % origami["num_staples"])
	print("- No tree: \t%d\t%s" % (countNoTree, str_no_tree))
	print("- 0 routes: \t%d\t%s" % (count0, str_zero_routes))
	print("- 1 route: \t%d" % count1)
	print("- >1 routes: \t%d\t%s" % (countG1, str_g1_routes))
	print("\tAverage number of routes on >1 route staples: \t%1.2f" % (G1av))
	print("\tMaximum number of routes on >1 route staples: \t%d" % (G1max))

	print("")

	if count1 == origami["num_staples"] :
		print("\n ** All staples placed on origami ** \n")
	else :
		print("%d staples currently placed on origami\n" % (staples_placed))
