"""Quantifies the staple split length ambiguity of a particular origami contact map, based on sequence information

Creates an "ambig" file, detailing by how many bases each edge end in the origami graph 
can be shortened (-) or elongated (+) [and no new graph edges added] 
while still ensuring scaffold-staples sequence complementarity

NOTE: it is assumed that the origami is in fully-assembled state, and no ss_domains connect to other ss_domains in the origami graph.

Ben Shirt-Ediss, November 2021
"""

import sys
sys.path.append("../")

from random import randrange
import networkx

from contactmap import iohandler
from contactmap import dna



MIN_STAPLE_SECTION_LENGTH = 5 	# bp
MIN_SS_DOMAIN_LENGTH = 1 		# nt




#
#
# LOCAL VALIDITY OF A GRAPH NODE MOVEMENT
#
#


def get_node_id_of_opposite_edge_end(edge, node_id) :

	e = [x for x in edge if x != node_id]
	return e[0]




def get_node_info(G, node_id) :

	if node_id not in G.nodes() :
		return False	

	node = {}

	node["connects_to"] = {}			# dict of edge types this node connects to
	node["domain_edge"] = None 			# (u,v) of the ds or ss scaffold domain connected to this node
	node["ssDNA_edge"] = None 			# (u,v) of the single-stranded dangle edge, or loopout edge connected to this node
	node["c_node_id"] = None 			# id of node on other side of crossover  (if it exists)
	node["n_node_id"] = None 			# id of node on other side of ghost nick (if it exists)
	node["staple_id"] = None 			# id of the staple the node is part of
	
	if G.nodes[node_id].get("staple_id", -1) != -1 :
		node["staple_id"] = G.nodes[node_id]["staple_id"]

	for neighbour_edge in G.edges(node_id) :
		
		etype = G.edges[neighbour_edge]["type"]
		
		node["connects_to"][etype] = True

		if etype in ["ss_domain", "ds_domain"] :
			node["domain_edge"] = neighbour_edge
		if etype in ["dangle", "loopout"] :
			node["ssDNA_edge"] = neighbour_edge
		if etype == "crossover" :
			node["c_node_id"] = get_node_id_of_opposite_edge_end(neighbour_edge, node_id)
		elif etype == "ghost_nick" :
			node["n_node_id"] = get_node_id_of_opposite_edge_end(neighbour_edge, node_id)

	return node




def num_succesive_WC_complement_bases(seq_top, seq_bottom) :
	"""Compares seq_top and seq_bottom base-for-base and returns
	how many complementary matches there are. The strands can be unequal lengths.
	"""

	complements = 0
	n = min(len(seq_top), len(seq_bottom))

	for i in range(0, n) :
		if dna.bases_are_WC_compliments(seq_top[i], seq_bottom[i]) :
			complements += 1
		else :
			break

	return complements




def is_locally_valid(G, node_change) :
	"""If node_change is locally valid, returns a list of child node_change tuples to also check
	Otherwise, returns False
	"""

	children = []

	node_id, change = node_change

	if change == 0 :
		return False 		# Fail reason: no change is specified! change < 0 or change > 0 is required

	node = get_node_info(G, node_id)

	if node == False :
		return False 		# Fail reason: node to move cannot be found in graph

	if ("dangle" in node["connects_to"]) or ("loopout" in node["connects_to"]) :
		# ----------------------------------------------------
		# A staple dangling end, or a loopout is being resized
		# ----------------------------------------------------

		# note that dangles and staple loopouts are treated exactly the same

		# Test 1. node must be on scaffold and have adjacent ghost nick to enable sliding
		if node["n_node_id"] == None :
			return False	# Fail reason: node is either not on scaffold at all (end of staple dangle), or is on 
							# scaffold, but has no adjacent ghost nick to enable sliding

		domain_edge = node["domain_edge"]
		ssDNA_edge = node["ssDNA_edge"]
		
		# Test 2. node movement must not make connected staple section on the scaffold too short
		minlen = 1
		if G.edges[domain_edge]["runlen"] > MIN_STAPLE_SECTION_LENGTH :
			minlen = MIN_STAPLE_SECTION_LENGTH

		if G.edges[domain_edge]["runlen"] + change < minlen :
			return False	# Fail reason: node movement makes a double stranded staple section too short

		# Test 3. node movement must not consume all of the dangle, or loopout bases, making it disappear. At least 1nt must remain
		if change > G.edges[ssDNA_edge]["runlen"] - 1 :
			return False

		# Test 4. if change > 0, the sequence following the nick must be complementary to sequence on the dangle/loopout
		if change > 0 :		
			across_nick_domain_edge = get_node_info(G, node["n_node_id"])["domain_edge"]
			
			across_nick_scaffold_seq = G.edges[across_nick_domain_edge]["scaffold_sequence53"][::-1]
			ssDNA_seq = G.edges[ssDNA_edge]["staple_sequence53"]

			if node_id == G.edges[domain_edge]["target"] :
				across_nick_scaffold_seq = G.edges[across_nick_domain_edge]["scaffold_sequence53"]
				ssDNA_seq = G.edges[ssDNA_edge]["staple_sequence53"][::-1]

			if change > num_succesive_WC_complement_bases(ssDNA_seq, across_nick_scaffold_seq) :
				return False	# Fail reason: cannot hybridise that much ssDNA to scaffold -- bases on dangle/loopout and scaffold become non-complementary too soon

		# this move is locally valid, prepare the "children" nodes affected
		children.append( (node["n_node_id"],-1*change) )

	elif "ds_domain" in node["connects_to"] :
		# --------------------------------------------------------------------------------------
		# A staple section double helix (not connected to a dangle or loopout) is being re-sized
		# --------------------------------------------------------------------------------------
	
		# Test 1. node must have BOTH crossover AND ghost_nick edges connecting it, to be moveable at all			
		if node["c_node_id"] == None or node["n_node_id"] == None :
			return False 	# Fail reason: node being moved is on a staple section, but does not have crossover AND ghost_nick edge connections

		domain_edge = node["domain_edge"]
		
		# Test 2. node movement must not make staple section too short
		minlen = 1 # if the edge is already too short in the contact map, don't let it go lower than 1bp
		if G.edges[domain_edge]["runlen"] > MIN_STAPLE_SECTION_LENGTH :
			minlen = MIN_STAPLE_SECTION_LENGTH

		if G.edges[domain_edge]["runlen"] + change < minlen :
			return False	# Fail reason: node movement makes a double stranded staple section too short

		# Test 3. if change > 0, sequence following the nick must be complementary to accomodate extra bases coming across the crossover
		"""
		Note: a change < 0 also requires checking sequence complementarity for a crossover. But this is not done here.
		A change > 0 will always mean a change > 0 for a child node. And this is where the change < 0 complementarity for this node is checked.
		"""
		if change > 0 :

			across_nick_domain_edge = get_node_info(G, node["n_node_id"])["domain_edge"]
			across_xover_domain_edge = get_node_info(G, node["c_node_id"])["domain_edge"]

			# 3.1 read the sequences across the nick, and across the crossover, the correct way around!
			# default: node_id moving is the SOURCE node of domain_edge (closest to scaffold 5')
			across_nick_scaffold_seq = G.edges[across_nick_domain_edge]["scaffold_sequence53"][::-1]
			across_xover_staple_seq = G.edges[across_xover_domain_edge]["staple_sequence53"]

			if node_id == G.edges[domain_edge]["target"] :
				# node_id moving is the TARGET node of domain_edge (closest to scaffold 3')
				across_nick_scaffold_seq = G.edges[across_nick_domain_edge]["scaffold_sequence53"]
				across_xover_staple_seq = G.edges[across_xover_domain_edge]["staple_sequence53"][::-1]

			# 3.2 find if sequence complementarity will LET the crossover slide that far
			if change > num_succesive_WC_complement_bases(across_xover_staple_seq, across_nick_scaffold_seq) :
				return False	# Fail reason: crossover cannot slide that far -- staple-scaffold bases become non-complementary too soon

		# this move is locally valid, prepare the "child" nodes affected
		children.append( (node["c_node_id"],-1*change) )
		children.append( (node["n_node_id"],-1*change) )

	elif "ss_domain" in node["connects_to"] :
		# ---------------------------------------------------
		# A single stranded scaffold domain is being re-sized
		# ---------------------------------------------------
		
		domain_edge = node["domain_edge"]

		# Test 1. node must have a ghost_nick edge connecting it,
		# and ds_domain on the other side of the ghost nick, to be moveable at all
		if node["n_node_id"] == None :
			return False 	# Fail reason: no ghost nick connected
		
		across_nick_domain_edge = get_node_info(G, node["n_node_id"])["domain_edge"]

		if G.edges[across_nick_domain_edge]["type"] == "ss_domain" :
			return False 	# Fail reason: ghost nick connects to another ss_domain, not a ds_domain

		# Test 2. node movement must not make single stranded domain too short
		if G.edges[domain_edge]["runlen"] + change < MIN_SS_DOMAIN_LENGTH :
			return False	# Fail reason: node movement gives a single stranded domain which is too short

		# this move is locally valid, prepare the "children" nodes affected
		children.append( (node["n_node_id"],-1*change) )

	return children














#
#
# GLOBAL VALIDITY OF A GRAPH NODE MOVEMENT
#
#


def is_globally_valid(G, node_change, visited) :
	"""** Recursive function **

	Returns True if this node is locally valid, and all of its children (dependent events) 
	are also locally valid, meaning that this node is globally valid

	Returns False if either this node is not locally valid, or if one of its children is not 
	locally valid.

	The argument "visited" is a list, and must be specified
	Lists are passed by reference in Python, and so changes made to visited here persist
	visited holds the order of graph nodes traversed in the depth-first constraints search
	The change to make to the contact map can be deduced from visited

	This function is built on the simpler template is_globally_valid__SIMPLE_VERSION()
	"""

	#print("node ", node_change)

	# BOTTOM OUT 1 (success)
	# stop this branch when node seen before (and thus is locally valid)
	if node_change in visited :
		visited.append(node_change)
		#print("** seen before, locally valid")
		return True

	locally_valid = is_locally_valid(G, node_change)

	# BOTTOM OUT 2 (fail)
	# stop this branch when node is locally invalid
	if not locally_valid :
		#print("** not locally valid")
		visited.append(node_change)
		return False

	# if here, this node is locally valid
	visited.append(node_change)

	# this node becomes GLOBALLY valid, if all its children are locally valid
	# to establish the truth of this statement, you have to RECURSE depth-first down all the subtrees until the children are BOTTOM OUTS (leaf nodes)
	children = locally_valid
	#print("children ", children)

	if len(children) == 1 :
		return is_globally_valid(G, children[0], visited)
	elif len(children) == 2 :
		return is_globally_valid(G, children[0], visited) and is_globally_valid(G, children[1], visited)

	"""
	Notes: 
		In Python, False and(statement 2) does not execute statement 2

		Return values are used to stop branching
		List arguments are used to record search tree structure made
	"""



	



	










#
#
# BASIC RECURSION EXAMPLE
#
#

def is_globally_valid__SIMPLE_VERSION(node, order) :
	"""
	A stripped down skeleton of the recursion loop used
	The order list is passed by reference, so changes to it inside here persist

	Usage:
	order = [] 			# order of nodes traversed
	valid = test(1, order)
	print(valid)
	print(order)
	"""

	print("on node %d" % node)

	# you stop a branch when a node is invalid (in this case, you are invalid if you have value 3)
	# BOTTOM OUT 1
	
	if node == 3 :
		order.append(node)
		return False

	#... or when a child has been seen before
	# BOTTOM OUT 2
	
	if node in order :
		order.append(node)
		return True

	order.append(node)

	# node not seen before, generate 2 children
	c1 = randrange(5)
	c2 = randrange(5)
	print("children are %d %d" % (c1, c2))

	# this node is globally valid, if all of its children are locally valid
	# to establish the truth of this statement, you have to recurse down all the levels of the tree until the leaf nodes are BOTTOM OUTS
	return test(c1, order) and test(c2, order)

















#
#
# MAKE AMBIG DICT FOR ORIGAMI
#
#


def build_ambig_dict(G) :
	"""Returns a dictionary: node_id: most - change, most + change

	Note that a change of +1 at node_id means that this strand can be lengthened by 1 bp/nt at the strand end where node_id is
	Note that a change of -1 at node_id means that this strand can be shortened by 1 bp/nt at the strand end where node_id is
	"""

	ambig = {}

	# go through all nodes of the origami graph, and try to move them - and +
	for node_id in G.nodes() :
		
		changeplus = 1
		while True :
			node_change = (node_id, changeplus)
			if not is_globally_valid(G, node_change, []) :
				changeplus -= 1
				break
			changeplus += 1
		
		changeminus = -1
		while True :
			node_change = (node_id, changeminus)
			if not is_globally_valid(G, node_change, []) :
				changeminus += 1
				break
			changeminus -= 1

		ambig[node_id] = (changeminus, changeplus)

	return ambig

















#
#
# STAPLE SECTION CHANGES TO MAKE IN CONTACT MAP, TO IMPLEMENT A GRAPH NODE CHANGE
#
#


def get_staple_section_changes_for_node_change(G, node_change) :
	# Note 1: the node_change is assumed globally valid

	# Note 2: movement of nodes on ss_domain's don't directly produce a staple change, and so are not handled here

	# Note 3: the code here has many repeats and is not optimal
	# however, I write out all the cases explicitly to ease understanding

	node_id, change = node_change

	node = get_node_info(G, node_id)

	staple_id = node["staple_id"]
	if staple_id == None :
		return {}

	staple_section_changes = {}		# (staple id, section id) : (new length, new base_i53)

	if ("dangle" in node["connects_to"]) or ("loopout" in node["connects_to"]) :
		# ----------------------------------------------------
		# A staple dangling end, or a loopout is being resized
		# (There are 2 staple section changes)
		# ----------------------------------------------------

		domain_edge = node["domain_edge"]
		ssDNA_edge = node["ssDNA_edge"]

		if change < 0 :
			# --- the scaffold-hybrised section of the staple becomes shorter, and the dangle/loopout section becomes longer ---

			change = abs(change)

			if node_id == G.edges[domain_edge]["source"] :		# dangle 5'-3' goes OUT of node

				# A. scaffold-hybrised section of the staple
				length = G.edges[domain_edge]["runlen"] - change
				base_i53 = G.edges[domain_edge]["base_i53"][change:]
				base_i53 = base_i53[::-1]
				# [change:] = list, with CHANGE bases taken off the start

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

				# B. dangle/loopout section of the staple
				length = G.edges[ssDNA_edge]["runlen"] + change
				base_i53 = [-1] * length

				staple_section_changes[ (staple_id, G.edges[ssDNA_edge]["section_id"]) ] = (length, base_i53)

			elif node_id == G.edges[domain_edge]["target"] :	# dangle 5'-3' goes INTO node

				# A. scaffold-hybrised section of the staple
				length = G.edges[domain_edge]["runlen"] - change
				base_i53 = G.edges[domain_edge]["base_i53"][:-change]
				base_i53 = base_i53[::-1]
				# [:-change] = list, with CHANGE bases taken off the end

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)				

				# B. dangle/loopout section of the staple
				length = G.edges[ssDNA_edge]["runlen"] + change
				base_i53 = [-1] * length

				staple_section_changes[ (staple_id, G.edges[ssDNA_edge]["section_id"]) ] = (length, base_i53)

		elif change > 0 :
			# --- the scaffold-hybrised section of the staple becomes longer, and the dangle/loopout section becomes shorter ---

			across_nick_domain_edge = get_node_info(G, node["n_node_id"])["domain_edge"]

			if node_id == G.edges[domain_edge]["source"] :		# dangle 5'-3' goes OUT of node
			
				# A. scaffold-hybrised section of the staple
				length = G.edges[domain_edge]["runlen"] + change
				base_i53 = G.edges[across_nick_domain_edge]["base_i53"][-change:] + G.edges[domain_edge]["base_i53"]
				base_i53 = base_i53[::-1]
				# [-change:] = last CHANGE bases of list

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

				# B. dangle/loopout section of the staple
				length = G.edges[ssDNA_edge]["runlen"] - change
				base_i53 = [-1] * length

				staple_section_changes[ (staple_id, G.edges[ssDNA_edge]["section_id"]) ] = (length, base_i53)

			elif node_id == G.edges[domain_edge]["target"] :	# dangle 5'-3' goes INTO node
				
				# A. scaffold-hybrised section of the staple
				length = G.edges[domain_edge]["runlen"] + change
				base_i53 = G.edges[domain_edge]["base_i53"] + G.edges[across_nick_domain_edge]["base_i53"][0:change]
				base_i53 = base_i53[::-1]
				# [0:change] = first CHANGE bases of list

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

				# B. dangle/loopout section of the staple
				length = G.edges[ssDNA_edge]["runlen"] - change
				base_i53 = [-1] * length

				staple_section_changes[ (staple_id, G.edges[ssDNA_edge]["section_id"]) ] = (length, base_i53)

	elif "ds_domain" in node["connects_to"] :
		# --------------------------------------------------------------------------------------
		# A staple section double helix (not connected to a dangle or loopout) is being re-sized
		# (There is 1 staple section change)
		# --------------------------------------------------------------------------------------

		domain_edge = node["domain_edge"]

		if change < 0 :

			change = abs(change)

			if node_id == G.edges[domain_edge]["source"] :		# crossover 5'-3' direction goes OUT of node

				length = G.edges[domain_edge]["runlen"] - change
				base_i53 = G.edges[domain_edge]["base_i53"][change:]
				base_i53 = base_i53[::-1]
				# [change:] = list, with CHANGE bases taken off the start

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

			elif node_id == G.edges[domain_edge]["target"] :	# crossover 5'-3' direction goes INTO node

				length = G.edges[domain_edge]["runlen"] - change
				base_i53 = G.edges[domain_edge]["base_i53"][:-change]
				base_i53 = base_i53[::-1]
				# [:-change] = list, with CHANGE bases taken off the end

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

		elif change > 0 :

			across_nick_domain_edge = get_node_info(G, node["n_node_id"])["domain_edge"]
	
			if node_id == G.edges[domain_edge]["source"] :		# crossover 5'-3' direction goes OUT of node

				length = G.edges[domain_edge]["runlen"] + change
				base_i53 = G.edges[across_nick_domain_edge]["base_i53"][-change:] + G.edges[domain_edge]["base_i53"]
				base_i53 = base_i53[::-1] 
				# [-change:] = last CHANGE bases of list

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

			elif node_id == G.edges[domain_edge]["target"] :    # crossover 5'-3' direction goes INTO node

				length = G.edges[domain_edge]["runlen"] + change
				base_i53 = G.edges[domain_edge]["base_i53"] + G.edges[across_nick_domain_edge]["base_i53"][0:change]
				base_i53 = base_i53[::-1] 
				# [0:change] = first CHANGE bases of list

				staple_section_changes[ (staple_id, G.edges[domain_edge]["section_id"]) ] = (length, base_i53)

	return staple_section_changes









