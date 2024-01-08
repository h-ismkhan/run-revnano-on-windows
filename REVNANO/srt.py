""" Staple routing tree operations

Ben Shirt-Ediss, 2020, 2023
"""

import sys
sys.path.append("../")

from REVNANO import parameters
from REVNANO import dna



TREE_HISTORY_INDEX = 0
TREE_RUNLEN_SO_FAR_INDEX = 1
TREE_IS_PRUNED_INDEX = 2





#
#
# FINDING SCAFFOLD MATCH REGIONS FOR A STAPLE FRAGMENT
#
#

def truncated_runlen(origami, i53_fragment_5prime, runlen, i53_of_bases_in_staple_so_far) :

	runlen_trunc = 0
	i53_of_bases_to_add_to_staple = set()

	for i53raw in range(i53_fragment_5prime, i53_fragment_5prime-runlen, -1) :

		i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53raw)	

		if i53 in i53_of_bases_in_staple_so_far :
			# base already on staple, cut off the match before this base
			return runlen_trunc, i53_of_bases_to_add_to_staple
		else :
			# base not seen before
			runlen_trunc += 1
			i53_of_bases_to_add_to_staple.add(i53)

	return runlen_trunc, i53_of_bases_to_add_to_staple




def get_match_candidates(origami, staple_id, subseq53, i53_of_bases_in_staple_so_far, SIGMA, MU_MIN, BETA) :

	match_candidates = {}		# i53 index: (runlen, set of bases in run)

	# 1. find length of consecutive matches from each base along the scaffold, for the staple fragment subseq53
	run = dna.get_run(origami, subseq53)
		# run is a dict. i35 on scaffold: number of consecutive bases complementary to subseq53 (going from 3' to 5' of scaffold)

	if len(run) == 0 :				# no more match candidates
		return match_candidates

	# ...if this far, there are complementary matches for subseq53 on the scaffold

	mu_max = max(run.values()) 				
	mu_threshold = mu_max - SIGMA
	if mu_threshold < MU_MIN :
		mu_threshold = MU_MIN

	# 2. validate the scaffold match regions: each must be long enough
	# (also check that none of the match regions cause staple self-intersection, and if so, truncate them)
	
	for i35, runlen in run.items() :
		# test 0: is match length in the acceptance window?
		if (runlen >= mu_threshold and runlen <= mu_max) :
			# yes

			# test 1: does this candidate intersect a previous part of this same staple?
			# if it does, truncate it. Otherwise, dont truncate
			i53 = origami["scaffold_nt"] - 1 - i35	 # i53 index of staple 5' (this needs to be decremented, to get to staple 3')
			runlen_trunc, i53_of_bases_to_add_to_staple = truncated_runlen(origami, i53, runlen, i53_of_bases_in_staple_so_far)

			# re-test if we still have a potential candidate
			if (runlen_trunc >= mu_threshold and runlen_trunc <= mu_max) :
				pass 		# yes
			else :
				continue	# no

			# test 2: does this candidate want to claim a scaffold region
			# that is already occupied by an already placed definite 1-route staple? 
			# or a region that is fully NESTED by an already definite 1-route staple? 

			if (dna.scaffold_segment_get_claimed_percent(origami, staple_id, i53, runlen_trunc) >= BETA) or (dna.scaffold_segment_fully_nests_a_staple_section(origami, staple_id, i53, runlen_trunc)) :
				continue 	# candidate not valid, don't explore it

			# candidate is fine
			match_candidates[i53] = (runlen_trunc, i53_of_bases_to_add_to_staple)

	return match_candidates



























#
#
# STAPLE ROUTING TREE GENERATION
#
#

def staple_has_1route_fast_check(origami, staple_id, subseq53, i53_of_bases_in_staple_so_far, runlen_so_far, staplelen, SIGMA, MU_MIN, BETA) :
	"""Recursive (self-calling) function. 

	Fast-check procedure to determine if a staple has a single route tree
	This function ONLY follows a single branch of the staple routing tree, thus is fast

	The staple does not need a "perfect" linear 1-route tree: dead-end children at 
	different levels are tolerated. This feature is essential to catch most 1-route trees

	Note: this function does not determine ALL 1-route trees, just a significant proportion of them
	To determine all single route trees entails constructing the whole tree (slow) and then
	determining if a single route exists pruning the whole tree down

	Note: tree leaf nodes are not fed into this function: recursion always stops
	1 level above tree leaf nodes

	Pseudocode:

		assess all children of node subseq53

		if all children have 0 children themselves:
			if exactly 1 child is a valid tree leaf (has cumulative runlen equal to staple length):
				RETURN TRUE  (bottom out - staple has 1 route)
			else :
				RETURN FALSE (bottom out - the staple route is singular, but not long enough)
		else:
			if exactly 1 child has further children AND the other children are not valid tree leaves
				(i.e. the other children do not have cumulative runlen equal to staple length)
					RECURSE (on the 1 child that has further children)
			else :
				RETURN FALSE (tree is probably multi-route)
	"""

	candidates_which_are_valid_tree_leaves = []
	candidates_with_further_children = []

	match_candidates = get_match_candidates(origami, staple_id, subseq53, i53_of_bases_in_staple_so_far, SIGMA, MU_MIN, BETA)

	# 1. get children of current node of singular staple route
	for i53 in match_candidates :

		runlen, i53_of_bases_to_add_to_staple = match_candidates[i53]

		# does this child represent a valid tree leaf? (i.e. it runs to the end of the staple)
		if runlen_so_far + runlen == staplelen :
			candidates_which_are_valid_tree_leaves.append(i53)

		# does this child have any children itself?
		subseq53_CHILD = subseq53[runlen:]
		i53_of_bases_in_staple_so_far_CHILD = i53_of_bases_in_staple_so_far.union(i53_of_bases_to_add_to_staple)

		match_candidates2 = get_match_candidates(origami, staple_id, subseq53_CHILD, i53_of_bases_in_staple_so_far_CHILD, SIGMA, MU_MIN, BETA) 
		if len(match_candidates2) > 0 :
			candidates_with_further_children.append(i53)

	# 2. decide to recurse or bottom out, depending on properties of the children
	if len(candidates_with_further_children) == 0 :
		if len(candidates_which_are_valid_tree_leaves) == 1 :
			return True  # recursion bottom out -- tree is single route
		else :
			return False # recursion bottom out -- there was a single route, but its not long enough
	else :
		if len(candidates_with_further_children) == 1 and \
			len(candidates_which_are_valid_tree_leaves) == 0 :		# note: if candidate has children, by definition it is not a tree leaf
				# no answer yet: recurse to go further down single tree branch

			runlen, i53_of_bases_to_add_to_staple = match_candidates[candidates_with_further_children[0]]
			subseq53_CHILD = subseq53[runlen:]
			i53_of_bases_in_staple_so_far_CHILD = i53_of_bases_in_staple_so_far.union(i53_of_bases_to_add_to_staple)

			return staple_has_1route_fast_check(origami, staple_id, subseq53_CHILD, i53_of_bases_in_staple_so_far_CHILD, runlen_so_far + runlen, staplelen, SIGMA, MU_MIN, BETA)				
		else :	
			return False	# recursion bottom out -- tree looks probable to be multi-route




subtree_count = 0

def generate_staple_routing_tree(origami, staple_id, subseq53, history, i53_of_bases_in_staple_so_far, runlen_so_far, staplelen, SIGMA, MU_MIN, BETA) :
	"""Recursive (self-calling) function. 

	Generates staple routing tree using depth-first recursion

	Specifically, takes a staple subsequence and returns a dictionary of child nodes (match sites for this subsequence)
	Each child node also has a sequence, and links to a dictionary of its own child nodes, and so on
	(This is how a tree is implemented, as a dictionary of nested dictionaries)

	A child nodes dictionary has the form:
	(i53, runlen) : [ path to this node, total staple length to this point, node pruned = True/False, subtree ]

	- The key is the binding site on the scaffold
	- Value is various information about the binding site (tree node), plus the dictionary of its own children

	When all recursions have been done on all subtrees, the full routing tree for a staple is formed

	Note 1: The staple tree itself does NOT have a single major root node, it starts with the children of the root
	(all binding positions of the first section)

	Note 2: Other arguments are also passed to this function which are
	not stored in the tree, but which carry information about the branch being constructed
	"""
	global subtree_count

	# -----------------------
	# 1 call per sub tree
	# if call count exceeds maximum, tree creation raises exception
	subtree_count += 1
	if subtree_count > parameters.MAX_SUBTREES :
		raise ValueError("Staple Routing Tree Too Large")
	# -----------------------

	match_candidates = get_match_candidates(origami, staple_id, subseq53, i53_of_bases_in_staple_so_far, SIGMA, MU_MIN, BETA)

	if len(match_candidates) == 0 : 	# i.e. subseq53 has no valid scaffold matches, there is no tree
		return None 					# this is where recursion bottoms out

	# 2. add each of these potential matches as children coming off the parent
	# (recursion)
	tree = {}

	for i53 in match_candidates :

		runlen, i53_of_bases_to_add_to_staple = match_candidates[i53]

		subseq53_CHILD = subseq53[runlen:]
		history_CHILD = list(history)				# deep copy of history to this point, to give to child node
		history_CHILD.append( (i53, runlen) )   	# history is a list of (i53, runlen) tuples up to the current node, and includes the current node
		i53_of_bases_in_staple_so_far_CHILD = i53_of_bases_in_staple_so_far.union(i53_of_bases_to_add_to_staple)
													# no need to deep copy the set, as union method returns a new set object

		# get subtree coming from this child node
		subtree = generate_staple_routing_tree(origami, staple_id, subseq53_CHILD, history_CHILD, i53_of_bases_in_staple_so_far_CHILD, runlen_so_far + runlen, staplelen, SIGMA, MU_MIN, BETA)

		if subtree == None :
			# this node is a leaf node
			
			if (runlen_so_far + runlen) < staplelen :
				# mark it as PRUNED: the path to the leaf falls SHORT of matching the entire staple sequence length
				tree[(i53, runlen)] = [ history_CHILD, runlen_so_far + runlen, True, subtree ]
			else :
				# valid leaf node
				tree[(i53, runlen)] = [ history_CHILD, runlen_so_far + runlen, False, subtree ]
		else :
			# this node is not a leaf node, there is a subtree coming from it
			tree[(i53, runlen)] = [ history_CHILD, runlen_so_far + runlen, False, subtree ]

	return tree




























#
#
# STAPLE ROUTING TREE PRUNING
#
#

def prune_down_staple_routing_tree(tree) :
	"""Recursive (self-calling) function. 

	The parent of tree is pruned, this method recursively marks all children as pruned
	"""

	if tree == None :
		return

	for child in tree :

		# prune this child node
		tree[child][TREE_IS_PRUNED_INDEX] = True

		# and prune the child's children
		history, runlen_so_far, is_pruned, subtree = tree[child]
		prune_down_staple_routing_tree(subtree)




def propagate_pruning_in_staple_routing_tree(tree) :
	"""Recursive (self-calling) function. 

	Given a tree, marks nodes as pruned if their children are all pruned
	and returns True / False if the parent of tree itself should be pruned

	Thus, propagates pruning UP the staple tree
	Also, propagates DOWN the staple tree: any pruned nodes also have their children marked as pruned

	Important Note: The initial tree supplied (before the recursive calls) should not be None
	"""
	
	if tree == None :
		return False  # an unpruned leaf node

	total_pruned = 0

	for child in tree :

		history, runlen_so_far, is_pruned, subtree = tree[child]

		if is_pruned :
			# child is already marked as pruned
			total_pruned += 1
			
			# Note: the following step is NOT strictly necessary, because trees are only ever queried from the ROOT, not the leaves
			# but I do it for completeness -- to make all pruned nodes appear in red when the graph is rendered

			# make sure that all of its children are, in turn, marked as pruned
			prune_down_staple_routing_tree(subtree)
		else :
			# don't know if child is pruned or not, we need to query its children to know this
			if propagate_pruning_in_staple_routing_tree(subtree) :
				# yes, in the end, the child is pruned (because all of its children are pruned)
				tree[child][TREE_IS_PRUNED_INDEX] = True
				total_pruned += 1				
			else :
				# no, the child is not pruned
				pass

	return total_pruned == len(tree)




def prune_blocked_nodes_of_staple_routing_tree(origami, staple_id, tree, BETA) :
	"""Recursive (self-calling) function. 

	Goes through staple routing tree, and marks nodes which are already covered by other staples as PRUNED
	(origami["footprints"] is used)

	Note: does not marked the children of a pruned node as also pruned: 
	propagate_pruning_in_staple_routing_tree() needs to be called afterwards for this
	"""

	if tree == None :
		return			# stop when leaf nodes reached

	for child in tree :

		i53, runlen = child
		history, runlen_so_far, is_pruned, subtree = tree[child]

		if is_pruned :
			continue
		
		if (dna.scaffold_segment_get_claimed_percent(origami, staple_id, i53, runlen) >= BETA) or (dna.scaffold_segment_fully_nests_a_staple_section(origami, staple_id, i53, runlen)) :
			
			# (1) this node (potential section) is already significantly claimed by other staples
			# or (2) this node claims bases that FULLY NEST the section of an existing staple (not allowed)
			# -- mark it as pruned, then skip it (its children will all be pruned later by calling the propagate_pruning_in_staple_routing_tree() function)
			tree[child][TREE_IS_PRUNED_INDEX] = True
		else :
			# this node is not sufficiently blocked by other staples (...yet)
			# -- see if any of its children are blocked
			prune_blocked_nodes_of_staple_routing_tree(origami, staple_id, subtree, BETA)




def prune_all_but_shortest_route_of_staple_routing_tree(tree, history_of_shortest_route) :
	"""Recursive (self-calling) function. 

	Used in Stage 2

	Marks as pruned child notes which are not on the shortest route history supplied

	Note: does not marked the children of a pruned node as also pruned: propagate_pruning_in_staple_routing_tree() needs to be used for this
	"""

	if tree == None :	# no children of parent
		return

	for child in tree :

		i53, runlen = child
		history, runlen_so_far, is_pruned, subtree = tree[child]

		if (i53, runlen) in history_of_shortest_route :
			# this child is on the shortest route, go ahead and check its children
			prune_all_but_shortest_route_of_staple_routing_tree(subtree, history_of_shortest_route)
		else :
			# child not on the shortest path -- prune it (and its children, grand children...)
			tree[child][TREE_IS_PRUNED_INDEX] = True





























#
#
# FOOTPRINT OF STAPLE ROUTING TREE ON SCAFFOLD
#
#

def staple_routing_tree_footprint(origami, tree, footprint) :
	"""Recursive (self-calling) function. 

	Sets the number of times that ALL valid (non-pruned) routes in staple routing tree "tree" cross each base position on the scaffold

	This is called the "footprint" of the tree

	If a base position is crossed by all possible routes of the staple, then
	the base position can be set as belonging to this staple.

	Returns the current number of non-pruned leaf nodes the staple tree has
	(e.g. the current number of valid ROUTES through the tree)	
	
	footprint is a dict 	i53: number of times crossed by staple routing tree
	footprint is modified by REFERENCE, it is also a return value
	"""
	
	if tree == None :
		return 0

	num_staple_routes = 0

	for child in tree :

		i53, runlen = child
		history, runlen_so_far, is_pruned, subtree = tree[child]

		if is_pruned :
			continue

		if subtree == None :
			# this is a leaf node, mark bases en-route to this node
			for h in history :
				for i53h in range(h[0], h[0]-h[1], -1) :
					i53h = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53h)
					footprint[i53h] += 1

			num_staple_routes += 1

		else :
			# not reached leaf yet
			num_staple_routes += staple_routing_tree_footprint(origami, subtree, footprint)

	return num_staple_routes
































#
#
# STAPLE ROUTING TREE STATS
#
#

def list_of_valid_leaf_histories(tree, leaf_histories) :
	"""Recursive (self-calling) function. 

	Finds all valid (non-pruned) leaf nodes of the staple routing tree "tree"
	Finds the route to each valid leaf node (stored in the leaf node history attribute)
	Returns all histories as a list
	"""

	for child in tree :

		i53, runlen = child
		history, runlen_so_far, is_pruned, subtree = tree[child]

		if is_pruned :
			continue
		elif subtree == None :	# at a non-pruned leaf node (valid route)
			leaf_histories.append( history )
		else :
			leaf_histories = list_of_valid_leaf_histories(subtree, leaf_histories)

	return leaf_histories




def routing_tree_stats(origami) :
	"""Over all staple staple routing trees, returns number of staples with 
	a) no tree, b) 0 valid routes, c) exactly 1 valid route, and c) >1 valid routes
	
	For staples with >1 valid routes, returns average number of routes and max number of routes
	"""

	countNoTree = 0
	count0 = 0
	count1 = 0
	countG1 = 0

	G1cumu = 0
	G1max = 0

	for staple_id, _ in enumerate(origami["staples53"]) :	

		rtree = origami["staple_routing_trees"][staple_id]
		if rtree == None :
			countNoTree += 1
			continue

		num_staple_routes = len(list_of_valid_leaf_histories(rtree, []))
			# this is only valid routes, which are not pruned

		if num_staple_routes == 0 : 	# this only applies to stages 1 and above, not stage 0
			count0 += 1
		elif num_staple_routes == 1 :
			count1 += 1			
		elif num_staple_routes > 1 :
			countG1 += 1
			G1cumu += num_staple_routes
			if num_staple_routes > G1max :
				G1max = num_staple_routes

	G1av = 0
	if countG1 > 0 :	# avoid division by 0 in calc average route number
		G1av = G1cumu/countG1

	return countNoTree, count0, count1, countG1, G1av, G1max




