"""Utilities for origami contact maps

Ben Shirt-Ediss, October 2021
"""

import sys
sys.path.append("../")


import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from itertools import groupby
import copy

from contactmap import iohandler
from contactmap import dna
from contactmap import contact2dlgraph
from contactmap import contact2ambig





#
#
# HELPER FUNCTIONS
#
#


def _scaffold_force_last_domain_id_to_zero(contactmap) :
	"""(This makes the last scaffold domain have id of zero (same as the first scaffold domain)
	"""

	_, sc, _ = contactmap

	dom_id_3prime = sc[len(sc)-1]["scaffold_domain_id"]

	i53 = len(sc)-1
	while sc[i53]["scaffold_domain_id"] == dom_id_3prime :
		sc[i53]["scaffold_domain_id"] = 0
		i53 -= 1




def _reset_sc_from_st(contactmap) :
	"""Resets the sc part of the contact map from the st part.

	Only the following fields of st are required, for each staple:

	"sequence53", "section_lengths", "base_i53"

	(note: the "num_sections" and "section_scaffold_domain_ids" fields can be empty, because they are re-calculated)
										 
	The following fields of sc are set, for each scaffold base:

	"scaffold_domain_id", "staple_id", "staple_section_id", "staple_base_id"

	Important Note: The scaffold sequence in sc CANNOT be inferred from staples in SINGLE STRANDED regions of the scaffold
	"""

	st, sc, scaffold_circular = contactmap

	# 1. clear sc dictionary, keeping only base_letter for each i53 index
	for i53 in sc :
		sc[i53]["scaffold_domain_id"] = -1
		sc[i53]["staple_id"] = -1
		sc[i53]["staple_section_id"] = -1
		sc[i53]["staple_base_id"] = -1

	# 2. fill staple details in sc, except for "scaffold_domain_id"
	for staple_id in st :

		# break base_i53 list up, so there is a sublist for each staple section
		b53_chunks = _chunk_staple_base_i53_by_staple_section_lengths(contactmap, staple_id)

		st[staple_id]["num_sections"] = len(b53_chunks)

		staple_base_id = 0

		for section_id, base_i53_list_for_section in enumerate(b53_chunks) :

			# go through all bases of staple section,
			for i53 in base_i53_list_for_section :
				
				if i53 == -1 :	# staple base not connected to scaffold
					staple_base_id += 1
					continue

				sc[i53]["staple_id"] = staple_id
				sc[i53]["staple_section_id"] = section_id
				sc[i53]["staple_base_id"] = staple_base_id

				staple_base_id += 1

	# 3. fill in "scaffold_domain_id" in sc, treating scaffold as linear (as in scadnano2contact.py)
	domain_id = 0
	ssid = "%d_%d" % (sc[0]["staple_id"], sc[0]["staple_section_id"])

	for i53 in range(0, len(sc)) :

		ssid_new = "%d_%d" % (sc[i53]["staple_id"], sc[i53]["staple_section_id"])

		if ssid_new != ssid :
			# staple section hybridised has changed
			domain_id += 1
			ssid = ssid_new

		sc[i53]["scaffold_domain_id"] = domain_id	

	# 4. add correction to domain id numbers, if scaffold is circular
	# if scaffold is circular, and the same staple section is hybridised to both ends of the scaffold, 
	# or if both ends of the scaffold are ssDNA,
	# then set the domain id of the last domain of scaffold to the domain id of the first domain (0)
	if scaffold_circular :

		bind3 = (sc[len(sc)-1]["staple_id"], sc[len(sc)-1]["staple_section_id"])
		bind5 = (sc[0]["staple_id"], sc[0]["staple_section_id"])

		if bind3 == bind5 :
			_scaffold_force_last_domain_id_to_zero(contactmap)


	# 5. fill in "section_scaffold_domain_ids" in st (as in scadnano2contact.py)
	for staple_id in st :

		staple = st[staple_id]

		# go through all staple bases in order, and find scaffold domain each hybridises to, from sc
		domain_id_list = []
		for base_i53 in staple["base_i53"] :

			# get scaffold domain id, of this staple base
			domain_id = -1 				# -1 i53's lead to -1 domain ids
			if base_i53 != -1 :
				domain_id = sc[base_i53]["scaffold_domain_id"]

			domain_id_list.append(domain_id)

		# compress groups of same domain ids in list, into a single number
		# https://stackoverflow.com/questions/39340345/how-to-count-consecutive-duplicates-in-a-python-list
		unique_domain_ids = [list(set(group))[0] for _, group in groupby(domain_id_list)]

		staple["section_scaffold_domain_ids"] = tuple(unique_domain_ids)




def _chunk_staple_base_i53_by_staple_section_lengths(contactmap, staple_id) :

	st, _, _ = contactmap

	b53 = st[staple_id]["base_i53"]
	slen = st[staple_id]["section_lengths"]

	# chunk the base_i53 list according to section lengths
	b53_chunks = []; i = 0
	for s in slen :
		b53_chunks.append( b53[i:i+s] )
		i = i+s

	return b53_chunks




def _make_staple_identifier(contactmap, staple_id) :
	"""Hashable unique identifier for a staple
	"""

	st, _, _ = contactmap

	iden = []
	for key in ["sequence53", "num_sections", "section_lengths", "section_scaffold_domain_ids", "base_i53"] :
		item = st[staple_id][key]
		if isinstance(item, list) :
			iden.append( tuple( item ) )
		else :
			iden.append( item )

	return tuple(iden)




def _make_diff(contactmap1, set1, set2) :
	"""Note: the diff can be made the other way around, by switching the contact maps in the arguments
	"""

	st1, _, _ = contactmap1

	# staples in origami1, not in origami2
	diff1 = {}		

	for tup in set1.difference(set2) :	# in set 1, not in set 2

		# get the sequence
		sequence53 = tup[0]

		# match it to a staple id in origami 1 contact map (assumes all staple sequences are unique)
		for staple_id in st1 :
			if st1[staple_id]["sequence53"] == sequence53 :
				break 	# has to be a match

		diff1[staple_id] = tup

	return diff1




MIN_STAPLE_SECTION_LENGTH = 7 	# linear scaffold cannot be placed in rotations where the
								# scaffold split breaks an existing staple section into shorter lengths than this


def _valid_rotation_offsets(offsets53) :

	valid = []  # valid offsets where 3' of linear scaffold can be placed 
						   # (5' of scaffold follows immediately after)
	
	for i in range(1, len(offsets53)) :
		
		first = offsets53[:i] 	# how the scaffold split breaks up the domain into two pieces
		last = offsets53[i:]

		if (len(first) >= MIN_STAPLE_SECTION_LENGTH) and (len(last) >= MIN_STAPLE_SECTION_LENGTH) :
			valid.append(first[-1])

	return valid




def _linear_scaffold_path(G, i53_scaffold_start, i53_scaffold_end) :
	"""Returns list
	[ (edge_type, i53_domain_start, i53_domain_end)... ]

	The i53_scaffold_start must be LESS THAN i53_scaffold_end
	These indexes identify the nodes at the beginning and end of the scaffold and are assumed to exist
	"""

	scaffold_path = []

	i53 = i53_scaffold_start
	prev_etype = None

	while i53 != i53_scaffold_end :
		
		# get all undirected edges connected to this i53 node
		for neighbour_edge in G.edges(i53) :

			etype = G.edges[neighbour_edge]["type"]

			# an ss_domain, ds_domain or ghost_nick edge will be followed from here to the next node on the scaffold
			# (note: etype is not allowed to be prev_etype, to stop back tracking: the graph is undirected)
			if (etype != prev_etype) and (etype in ["ss_domain", "ds_domain", "ghost_nick"]) :
				prev_etype = etype
				i53 = neighbour_edge[1]

				scaffold_path.append( (etype, neighbour_edge[0], neighbour_edge[1] ) )
									# edge_type, i53_domain_start, i53_domain_end
				continue

	return scaffold_path


























#
#
# SCAFFOLD OPS
#
#


def change_scaffold_sequence(contactmap, scaffold_material, staple_material, scaffold53_new) :
	"""Changes the scaffold sequence in a contact map, starting at scaffold 5', and re-derives the complementary staple sequences
	Staple dangles or loopouts that do not contact the scaffold retain their current sequences

	New June 2022: If the scaffold sequence given is longer than the scaffold in the contact map, and the scaffold is CIRCULAR, 
	then an ssDNA scaffold domain is added to the 3' end of the scaffold in the contact map, to accommodate the extra bases
	"""

	# contact map is changed by reference, so no return value
	st, sc, scaffold_circular = contactmap

	# 1. initialise a new staple set. New staples originally have question marks ? for bases
	staple_set = {}
	for staple_id in st :
		staple_set[staple_id] = ["?"] * len(st[staple_id]["sequence53"])

	# 2. assign new sequence to scaffold, and derive a new staple sequences from the new scaffold
	scaffold53_new = list(scaffold53_new)
	for i53 in sc :

		sc[i53]["base_letter"] = scaffold53_new[i53]

		staple_id = sc[i53]["staple_id"]
		staple_base_id = sc[i53]["staple_base_id"]

		if staple_id != -1 :	# check scaffold base is not single stranded, and is mapped to a staple

			staple_set[staple_id][staple_base_id] = \
				dna.complementary_staple_base(sc[i53]["base_letter"], scaffold_material, staple_material, watson_crick = True)[0]

	# 3. copy across new staple bases, set by the new scaffold, to the previous staples
	#	 ? bases remaining in the new staples are not copied across to previous staples, and the previous staples keep these bases the same
	for staple_id in st :
		oldseq53 = list(st[staple_id]["sequence53"])	# a sequence list
		newseq53 = staple_set[staple_id]				# a sequence list

		for i, base in enumerate(newseq53) :
			if base != "?" :
				oldseq53[i] = base

		st[staple_id]["sequence53"] = "".join(oldseq53)

	# 4. if the scaffold sequence being set is longer than the scaffold in the contact map,
	# and the scaffold is LINEAR, then
	# add an extra single stranded domain to the 3' scaffold in the contact map
	extra_bases = len(scaffold53_new) - len(sc)
	if extra_bases > 0 and (not scaffold_circular) :

		# the extra bases to add to 3' of the existing scaffold
		extra53 = scaffold53_new[len(sc):len(sc)+extra_bases]

		final_dom_id = sc[len(sc)-1]["scaffold_domain_id"] + 1

		i53 = len(sc)
		for i in range(0, extra_bases) :
			
			sc[i53] = 		{"base_letter"			: extra53[i],
							 "scaffold_domain_id"	: final_dom_id,
							 "staple_id"			: -1,
							 "staple_section_id"	: -1,
							 "staple_base_id"		: -1}
			i53 += 1

		# note: if necessary, a linear scaffold should be extended
		# before being circularised




def make_scaffold_circular(contactmap) :
	"""Returns a contact map with the original LINEAR scaffold changed to a circular scaffold
	
	Constructs a domain-level graph of the contact map, to check the ends.
	
	Note: Contact maps with 1nt domains cannot have the domain-level graph made, and thus cannot be circularised
	"""
	
	st, sc, scaffold_circular = contactmap

	if scaffold_circular :
		return False 		# already circular

	# see what the ending domains of the scaffold are connected to
	bind3 = (sc[len(sc)-1]["staple_id"], sc[len(sc)-1]["staple_section_id"])
	bind5 = (sc[0]["staple_id"], sc[0]["staple_section_id"])

	# Does the same staple occupy both scaffold ends?
	if bind3[0] == bind5[0] :
		# Maybe...

		if bind3 == (-1,-1) :
			# No: the same staple does NOT occupy both scaffold ends
			# Both ends of the scaffold are ssDNA
			# --> Case 2: circularisation just requires making domain ids at both ends of scaffold the same

			contactmap[2] = True
			_scaffold_force_last_domain_id_to_zero(contactmap)

		else :
			# Yes: The same staple DOES occupy both scaffold ends
			# Does this staple connect the two scaffold ends by a crossover?

			G = contact2dlgraph.origami_dl_graph(contactmap)
			u = 0
			v = len(sc)-1

			if G.has_edge(u, v) and G.edges[u,v]["type"] == "crossover" :
				# Yes. The staple connects the scaffold ends with a crossover
				# --> Case 3: circularisation requires to MERGE the two sections of the crossover staple into one
	
				"""
				How attributes change of staple that links scaffold ends with a crossover

				"sequence53" 						does not change	
				"num_sections"		 				minus 1
				"section_lengths" 				    remove section that hybridises at 3'
													add length of above section, to section that hybridises at 5'
				"section_scaffold_domain_ids" 	    (_reset_sc_from_st() handles updating this)
				"base_i53"							does not change	(just the crossover is removed)
				"""
				
				staple_id = bind5[0] 	# stame staple id is hybridised at both scaffold ends

				staple_section_id_3prime = bind3[1]
				staple_section_id_5prime = bind5[1]

				# delete section of the staple bound to scaffold 3', and add it's length to preceeding staple section, bound to scaffold 5'
				st[staple_id]["num_sections"] -= 1

				slen = st[staple_id]["section_lengths"][staple_section_id_3prime]
				del st[staple_id]["section_lengths"][staple_section_id_3prime]
				st[staple_id]["section_lengths"][staple_section_id_5prime] += slen

				contactmap[2] = True
				_reset_sc_from_st(contactmap)

			elif G.has_edge(u, v) and G.edges[u,v]["type"] == "loopout" :

				# a staple loopout across scaffold ends PREVENTS a scaffold being circularised
				# reason: a loopout is a single crossover edge in the DL graph, just like a standard crossover
				# 1) it cannot be "joined" like a crossover and 2) a nick edge and a loopout edge cannot
				# co-exit across the scaffold ends in the DL graph
				return False

			else :
				# No: The staple does not have a crossover across the scaffold ends
				# It just has different sections at the scaffold ends.
				# --> Case 1: circularisation is trivial. Just requires ends joining with a nick (i.e. simply changing circular flag to true)

				contactmap[2] = True

	else :
		# No. The same staple does not occupy scaffold ends. Scaffold ends are "distinct domains".
		# --> Case 1: circularisation is trivial. Just requires ends joining with a nick (i.e. simply changing circular flag to true)
		# Examples: Start domain-End domain is;  ssDNA-staple 1 section, staple 4 section0-staple 5 section2

		contactmap[2] = True		

	return True




def make_scaffold_linear(contactmap, i53_scaffold_3prime) :
	"""The i53 position of the new physical scaffold 3' is supplied. 
	The physical scaffold 5' end follows at the next position.
	A circular scaffold can be nicked at any position, but some nick points will cause 1bp domains
	"""

	st, sc, scaffold_circular = contactmap

	if scaffold_circular == False:
		return False 		# already linear

	# from here, we know scaffold is circular...

	i53_3prime = len(sc)-1

	# calculate other side of scaffold nick point   --->3' 5'---->
	i53_scaffold_5prime = i53_scaffold_3prime + 1
	if i53_scaffold_3prime == i53_3prime :
		i53_scaffold_5prime = 0

	# ------------
	# 1. see if the new scaffold nick will break a hybridised staple section
	bind3 = (sc[i53_scaffold_3prime]["staple_id"], sc[i53_scaffold_3prime]["staple_section_id"])
	bind5 = (sc[i53_scaffold_5prime]["staple_id"], sc[i53_scaffold_5prime]["staple_section_id"])

	if (bind3 == bind5) and (bind3 != (-1,-1)) :
		# yes, the hybridised staple section needs dividing in two

		"""
		How attributes change of staple that is split, by the new scaffold nick

		"sequence53" 						does not change	
		"num_sections"		 				add 1
		"section_lengths" 				    add new section, caused by split
		"section_scaffold_domain_ids" 	    (_reset_sc_from_st() handles updating this)
		"base_i53"							does not change	
		"""

		staple_id = bind5[0] 	# this is the staple to split
		section_id = bind5[1]

		# find exactly how i53_scaffold_3prime splits this staple section

		b53 = _chunk_staple_base_i53_by_staple_section_lengths(contactmap, staple_id)[section_id]	# the scaffold bases hybridised to, in the 5'-3' STAPLE direction
		b35 = b53[::-1]							# list the bases in the direction of the 5'-3' SCAFFOLD direction
		idx = b35.index(i53_scaffold_3prime)	# index of the 3' scaffold end
		sec1 = b35[idx+1:][::-1] 				# how the scaffold nick point divides the original staple section..
		sec2 = b35[:idx+1][::-1]				# going 5'-3' on the staple, sec1 is followed by sec2

		# update attributes of the staple split

		st[staple_id]["num_sections"] += 1

		st[staple_id]["section_lengths"][section_id] = len(sec2) 		# modify length of existing section
		st[staple_id]["section_lengths"].insert(section_id, len(sec1))	# add length of new section before it		
		
		# Note: this change above is not valid, until all the i53 indexes below are changed
		# (its not valid, because its a crossover added on a straight helix, which is not allowed
		# but when the index numbers change, its a crossover across the scaffold ends, which is allowed)

		if len(sec1) == 1 or len(sec2) == 1 :
			print("\nWarning: Placing the scaffold nick at this location creates a 1 base pair domain. " + \
			"The contact map is still produced, but it cannot be converted to a domain-level graph.\n")

	# ------------
	# 2. calculate how all i53 indexes in the contact map need to change, 
	# to implement the new PHYSICAL 5' and 3' ends of the linear scaffold 

	i53_old = list( range(0, i53_3prime+1) )
	i53_new = i53_old[-i53_scaffold_5prime:] + i53_old[:-i53_scaffold_5prime]

	i53_replace = {}	# old i53 index: new i53 index
	for i, i53 in enumerate(i53_old) :
		i53_replace[i53] = i53_new[i]

	# ------------
	# 3. change all i53 indexes in st and sc
	for staple_id in st :
		for i, i53 in enumerate(st[staple_id]["base_i53"]) :
			if i53 != -1 :
				st[staple_id]["base_i53"][i] = i53_replace[i53]

	sc2 = {}
	for i53 in sc :
		sc2[ i53_replace[i53] ] = sc[i53]
	contactmap[1] = sc2

	# ------------
	# 4. mark contact map as linear
	contactmap[2] = False

	# ------------
	# 5. infer the sc part of contact map from the st part
	_reset_sc_from_st(contactmap)

	return True



def apply_scaffold_rotation_offset(contactmap, i53_scaffold_3prime, scaffold_material, staple_material, scaffold53 = None) :

	st, sc, scaffold_circular = contactmap

	"""
	Technical note: 
	
	Circular scaffolds have their sequences rotated to allow ALL rotations to be accessed. The other option of 
	converting a circular contact map to a linear one of desired rotation, and then back to a circular one is not feasible, since the
	rotations that cause 1bp domains in the linear contact map cause problems.
	"""

	if scaffold53 != None :	
		# change scaffold sequence to the one specified
		change_scaffold_sequence(contactmap, scaffold_material, staple_material, scaffold53)
	else :
		# get existing scaffold sequence in contact map
		scaffold53 = "".join([sc[i53]["base_letter"] for i53 in sorted(sc.keys())])


	if scaffold_circular :
		# CIRCULAR SCAFFOLD

		# rotate scaffold sequence itself and re-apply it
	
		# calculate other side of scaffold nick point   --->3' 5'---->
		i53_scaffold_5prime = i53_scaffold_3prime + 1
		if i53_scaffold_3prime == len(sc)-1 :
			i53_scaffold_5prime = 0

		scaffold53_rotated = scaffold53[-i53_scaffold_5prime:] + scaffold53[:-i53_scaffold_5prime]

		change_scaffold_sequence(contactmap, scaffold_material, staple_material, scaffold53_rotated)

	else :
		# LINEAR SCAFFOLD

		# make linear scaffold circular
		success = make_scaffold_circular(contactmap)
		if not success :
			return False  			# (can fail if e.g. loopout is across scaffold ends)

		# then nick it in the new position (new 3' and 5' ends appear; scaffold sequence changes when listed 5' to 3'; staple sequences remain the same)
		make_scaffold_linear(contactmap, i53_scaffold_3prime)

		# finally apply the original scaffold sequence to the scaffold in its new position
		# (this gives new staple sequences which pin the scaffold into this rotation)
		change_scaffold_sequence(contactmap, scaffold_material, staple_material, scaffold53)

	return True



def get_thermodynamically_valid_scaffold_rotation_offsets(contactmap) :
	"""For the current origami design, returns list of offsets (base positions of 3' end) of all the scaffold rotations that are valid thermodynamically
	- For circular scaffolds, all rotation offsets are possible
	- For linear scaffolds, only those rotation offsets that do no create very short staple sections are permissable
	"""

	st, sc, scaffold_circular = contactmap

	if scaffold_circular == True :
		# each position along the circular scaffold is a valid rotation offset
		return list( range(0, len(sc)) )

	else :
		# linear scaffold  (the 5' and 3' scaffold ends do exist)
		
		# 1. check if rotation of scaffold is POSSIBLE at all

		G = contact2dlgraph.origami_dl_graph(contactmap)

		i53_3prime = len(sc)-1		

		# a crossover MUST link scaffold ends, in order for scaffold to rotate
		can_rotate = (G.has_edge(0, i53_3prime) and G.edges[0, i53_3prime]["type"] == "crossover")

		if can_rotate == False :
			return [i53_3prime]		# the only rotation possible is just where the scaffold is now

		# 2. work out valid rotation offsets of scaffold, for this origami design
		# note that rotation offsets are *relative to* the original contact map, supplied to this function

		# the scaffold is made into a linear path, and all double stranded domains along that path
		# have indexes marked where they can be split

		scaffold_path = _linear_scaffold_path(G, 0, i53_3prime)

		rotation_offsets = set()
		for edge in scaffold_path :
			etype, i53_domain_start, i53_domain_end = edge

			if etype == "ds_domain" :
				offsets53 = list( range(i53_domain_start, i53_domain_end+1) )
				rotation_offsets.update( _valid_rotation_offsets(offsets53) )

		# 3. Remember that a staple section bridges the first and last ds_domain's on the scaffold
		# the scaffold split point can ALSO be at multiple locations under this staple section

		offsets53 = list( range( scaffold_path[-1][1], scaffold_path[-1][2]+1 )) + \
					list( range( scaffold_path[0][1], scaffold_path[0][2]+1 )) 
					 
		rotation_offsets.update( _valid_rotation_offsets(offsets53) )

		# 4. add the rotation offset present in the original design (which may be illegal according to MIN_STAPLE_SECTION_LENGTH, but its always included anyway)
		rotation_offsets.add(i53_3prime)

		return sorted(list(rotation_offsets))	






























#
#
# 1nt SCAFFOLD DOMAIN DETECTION AND REMOVAL
#		(A contactmap with 1nt domains CANNOT be converted to a domain-level graph and cannot have a guide schematic generated)
# 		(Wireframe origamis can have 1nt scaffold domains at junctions)

def scaffold_1nt_domain_locations(contactmap) :
	""" Returns the i53 indexes of all 1nt domains on the scaffold
	This takes into account if the scaffold is circular
	"""

	st, sc, scaffold_circular = contactmap

	# get i53 indexes which serve as single base pair domains

	scaffold_nt = len(sc)
	i53_start_domain = []   	# bases which start and end domains are nodes on the scaffold graph
	i53_end_domain = []

	for i53 in sorted(sc.keys()) :		# go from 5' to 3' base on scaffold

		if i53 == 0 :						# vertex for first base
			i53_start_domain.append(i53)

		if i53 == scaffold_nt - 1 :			# vertex for last base
			i53_end_domain.append(i53)
			break

		domain_id = sc[i53]["scaffold_domain_id"]
		next_domain_id = sc[i53+1]["scaffold_domain_id"]

		if domain_id != next_domain_id :
			i53_end_domain.append(i53)
			i53_start_domain.append(i53+1)

	I = set.intersection(set(i53_start_domain), set(i53_end_domain))	
		# I is a set containing i53 indexes of all 1nt/1bp scaffold domains

	# circular scaffold correction::
	# if there is a 1nt domain on the first or last base of the scaffold, and the scaffold is circular, 
	# do and extra check of adjacent base around the loop
	if scaffold_circular :

		last = scaffold_nt-1

		if (0 in I) and (sc[0]["scaffold_domain_id"] == sc[last]["scaffold_domain_id"]):
			I.remove(0)	# not really a 1nt domain at virtual beginning of scaffold

		if (last in I) and (sc[0]["scaffold_domain_id"] == sc[last]["scaffold_domain_id"]):
			I.remove(last)	# not really a 1nt domain at virtual end of scaffold

	return I




def remove_staples_to_remove_1nt_scaffold_domains(contactmap) :
	""" Removing staples to remedy 1nt scaffold domains means that calculations of junction ambiguity are still correct
	BUT removing staples can deform crucial parts of the structure in the guide schematic
	"""

	I = scaffold_1nt_domain_locations(contactmap)

	_, sc, _ = contactmap
	scaffold_nt = len(sc)

	contactmapR = copy.deepcopy(contactmap) 	# this "guide schematic display" contact map will have staples removed, so no 1bp or 1nt domains exist on the scaffold
	staple_ids_removed = []

	for i53 in I :

		if sc[i53]["staple_id"] != -1 :
			# the 1nt scaffold domain at i53 is hybridised to a staple, believe it or not....
			raise ValueError("Error. 1nt scaffold domain hybridised to a staple found (i53 = %d).\n" % i53)

		# the 1nt scaffold domain is not hybridised
		# remove a neighbour staple, one side of the 1nt unhybridised scaffold domain, to get rid of the 1nt unhybridised domain
		if i53 == scaffold_nt - 1 :
			# last scaffold base
			# remove the staple hybridised to the previous scaffold base	
			neighbour_staple_id = sc[i53-1]["staple_id"]
		else :
			# not last scaffold base
			# remove the staple hybridised to the following scaffold base
			neighbour_staple_id = sc[i53+1]["staple_id"]

		staple_ids_removed.append(neighbour_staple_id)

	delete_staples(contactmapR, staple_ids_removed)

	return contactmapR, staple_ids_removed










def delete_1nt_scaffold_domains(contactmap) :
	""" Note: deleting 1nt scaffold domains will allow an origami a DL graph, and thus allow it to be displayed
	but it also changes the specification of the origami 
	"""

	# deleting unhybridised 1nt scaffold domains does change
	# staple["section_scaffold_domain_ids"], staple["base_i53"]   (both tuples)
	#
	# and it does NOT change
	# staple["sequence53"], staple["num_sections"], staple["section_lengths"]


	# 1. get base indexes of 1nt domains
	I = scaffold_1nt_domain_locations(contactmap)

	st, sc, scaffold_circular = contactmap

	# 2. work out mappings from old (pre-deletions) to new (post-deletions) values
	i53_change_map = { -1: -1 } 		# old i53 index -> new i53 index
	dom_change_map = { -1: -1 } 		# old scaffold dom id -> new scaffold dom id

		# -1 bases and domains on staples are not connected to scaffold, and remain not-connected during this procedure
		# hence I leave -1 mapping to itself

	deletions_so_far = 0
	for i53 in sorted(sc.keys()) :
			
		if i53 in I :

			if sc[i53]["staple_id"] != -1 :
				# the 1nt scaffold domain at i53 is hybridised to a staple, believe it or not....
				raise ValueError("Error. 1nt scaffold domain hybridised to a staple found (i53 = %d).\n" % i53)

			# skip over this i53

			deletions_so_far += 1 		# the i53 and dom_id of the deleted scaffold base wont be used by any staple, because it is not hybridised
		else :
			i53_change_map[i53] = i53 - deletions_so_far

			# get domain id of the current scaffold base
			dom_id = sc[i53]["scaffold_domain_id"]

			# update the domain id, only if the domain is NOT domain 0
			if dom_id == 0 :
				dom_change_map[dom_id] = 0
			else :
				dom_change_map[dom_id] = dom_id - deletions_so_far

			# note :
			# If a 1nt domain is the first scaffold base (domain 0 of scaffold), then this else loop is only entered at domain 1
			# BUT this else can be entered by domain id 0, if the scaffold is circular, and the virtual scaffold nick point falls inside the same domain
			# in this case, domain 0 should NOT be updated, as doing so sets it to -1 * deletions_so_far


	# 3. make the staples part of new contact map
	st2 = {}
	for staple_id in st :

		staple = st[staple_id]

		# most the staple info does not change
		st2[staple_id] = staple

		# update staple base i53 indexes (old replaced by new; -1 replaced by -1)
		st2[staple_id]["base_i53"] = tuple( [ i53_change_map[i53] for i53 in staple["base_i53"] ] )

		# update staple domain ids
		st2[staple_id]["section_scaffold_domain_ids"] = tuple( [ dom_change_map[dom_id] for dom_id in staple["section_scaffold_domain_ids"] ] ) 

		# other attributes of staple (sequence, number of sections, section lengths) are left unchanged

	# 4. make the scaffold part of new contact map
	sc2 = {}
	for i53 in sorted(sc.keys()) :

		if i53 not in i53_change_map :		# this scaffold base was deleted, and is not carried through to sc2
			continue

		# update i53 of scaffold base, if necessary
		i53_new = i53_change_map[i53]
		sc2[i53_new] = sc[i53]

		# update domain id of scaffold base
		dom_id = sc[i53]["scaffold_domain_id"]
		sc2[i53_new]["scaffold_domain_id"] = dom_change_map[dom_id]

		# other attributes of scaffold base (letter, staple bound to) are left unchanged

	contactmapD = [st2, sc2, scaffold_circular]

	return contactmapD

	


































#
#
# STAPLE OPS
#
#


def delete_staples(contactmap, staple_id_list) :
	"""Deletes staples in contact map and updates it by reference
	Returns true if all staples specified are deleted.
	Returns false if one or more of the staples did not exist to be deleted.
	
	Note: The staple id's are not renumbered - there will be gaps where staples were deleted
	"""

	st, _, _ = contactmap

	all_deleted = True

	for staple_id in staple_id_list :	
		if staple_id not in st :
			all_deleted = False
			continue

		del st[staple_id]

	_reset_sc_from_st(contactmap)

	return all_deleted





def add_staple_end_dangles(contactmap, staple_id, dangle5, dangle3) :

	# note: dangle sequences are specified 5' to 3' direction. No dangle is specified by passing the empty string
	# note: contactmap is changed by reference

	st, _, _ = contactmap

	staple = st[staple_id]

	staple["sequence53"] = dangle5 + staple["sequence53"] + dangle3
	
	if len(dangle5) > 0 :
		staple["num_sections"] += 1

		section_lengths = list( staple["section_lengths"] )
		section_lengths.insert(0, len(dangle5))
		staple["section_lengths"] = tuple( section_lengths )

		section_scaffold_domain_ids = list( staple["section_scaffold_domain_ids"] )
		section_scaffold_domain_ids.insert(0, -1)
		staple["section_scaffold_domain_ids"] = tuple( section_scaffold_domain_ids )

		base_i53 = list( staple["base_i53"] )
		base_i53 = ([-1] * len(dangle5)) + base_i53
		staple["base_i53"] = tuple( base_i53 )

	if len(dangle3) > 0 :
		staple["num_sections"] += 1

		section_lengths = list( staple["section_lengths"] )
		section_lengths.append(len(dangle3))
		staple["section_lengths"] = tuple( section_lengths )		

		section_scaffold_domain_ids = list( staple["section_scaffold_domain_ids"] )
		section_scaffold_domain_ids.append(-1)
		staple["section_scaffold_domain_ids"] = tuple( section_scaffold_domain_ids )

		base_i53 = list( staple["base_i53"] )
		base_i53 = base_i53 + ([-1] * len(dangle3))
		staple["base_i53"] = tuple( base_i53 )

	# note: adding sections to a staple at 5', shifts base indexes
	# therefore, the sc part of the contact map needs regenerating from the st part
	_reset_sc_from_st(contactmap)





def add_staple_end_dangles_from_file(staple_dangles_filename, contactmap) :
	"""Given a text file, with one line per staple in the format
	old_sequence new_sequence

	searches for staple with old_sequence, and replaces it with new_sequence, 
	where new_sequence is longer than old_sequence, contains it, and involves adding a beginning and/or end dangle to old_sequence
	"""

	dangles = iohandler.read_staple_dangles_file(staple_dangles_filename)

	if isinstance(dangles, str) :
		# error occurred reading the dangles file
		print(dangles)
		return False

	st, _, _ = contactmap

	for old_sequence in dangles :

		# find FIRST staple_id in contactmap with old_sequence
		# and add 5' and/or 3' dangles to it
		for staple_id in st :
			if st[staple_id]["sequence53"] == old_sequence :
				dangle5, dangle3 = dangles[old_sequence]
				add_staple_end_dangles(contactmap, staple_id, dangle5, dangle3)

	return True




def recompose_staple_with_loopouts(contactmap, staple_id, extra_staples_list, loopout_seq_list) :

	st, _, _ = contactmap

	staple = st[staple_id]

	# update attributes of staple (which is currently a stub)
	sequence53 = staple["sequence53"]
	num_sections = staple["num_sections"]
	section_lengths = list( staple["section_lengths"] )
	section_scaffold_domain_ids = list( staple["section_scaffold_domain_ids"] )
	base_i53 = list( staple["base_i53"] )

	# add loopout, add substaple, in pairs, ending with substaple LLLLLSSSSSLLLLLLSSSSS
	for i, loopout_seq in enumerate(loopout_seq_list) :
		
		# add loopout
		sequence53 += loopout_seq
		num_sections += 1 	# a loopout is always 1 single staple section
		section_lengths.append(len(loopout_seq))
		section_scaffold_domain_ids.append(-1)
		base_i53 += ([-1] * len(loopout_seq))

		# add substaple
		substaple = st[extra_staples_list[i]]
		
		sequence53 += substaple["sequence53"]
		num_sections += substaple["num_sections"]
		section_lengths += list(substaple["section_lengths"])
		section_scaffold_domain_ids += list(substaple["section_scaffold_domain_ids"])
		base_i53 += list(substaple["base_i53"])

	staple["sequence53"] = sequence53
	staple["num_sections"] = num_sections
	staple["section_lengths"] = tuple(section_lengths)
	staple["section_scaffold_domain_ids"] = tuple(section_scaffold_domain_ids)
	staple["base_i53"] = tuple(base_i53)

	# delete the virtual substaples, which are now all in the recomposed loopout staple
	delete_staples(contactmap, extra_staples_list)

	# note: the above function also regenerates the sc part of the contact map, as staples have been deleted















#
#
# MOVE AMBIGUOUS CROSSOVERS, LOOPOUTS, DANGLES
#
#

def move_scaffold_base(contactmap, i53, change) :

	# construct the domain-level origami graph
	G = contact2dlgraph.origami_dl_graph(contactmap)

	node_change = (i53, change)
	visited = []

	# check proposed move is legal, from structural and sequence matching perspective
	if contact2ambig.is_globally_valid(G, node_change, visited)	:

		# 1. get set of unique local moves to make, that implement the proposed move
		local = set(visited)

		# 2. get the staple section changes necessary for each local move, and put them all together
		staple_section_changes = {}
		for node_change in local :
			ssc = contact2ambig.get_staple_section_changes_for_node_change(G, node_change)
			staple_section_changes.update(ssc)

		# 3. update st of the contact map with the staple section changes
		#	 (note: staple sections are not being deleted or added: only re-sized)
		st, _, _ = contactmap
		for (staple_id, section_id) in staple_section_changes :
	
			new_length, new_base_i53 = staple_section_changes[ (staple_id, section_id) ]

			# update the base_i53 list for staple_id
			b53_chunks = _chunk_staple_base_i53_by_staple_section_lengths(contactmap, staple_id)
			b53_chunks[section_id] = new_base_i53
			st[staple_id]["base_i53"] = [item for sublist in b53_chunks for item in sublist]	# flatten list of lists into list

			# update the section lengths for staple_id
			st[staple_id]["section_lengths"][section_id] = new_length

		# 4. reset the sc part of the contact map from the freshly updated st part
		_reset_sc_from_st(contactmap)
		
		return True

	else :
		return False






	



















#
#
# SEQUENCE OPS
#
#

def change_DNA_to_RNA(contactmap) :
	"""All DNA T's in contact map are changed for RNA U's
	"""

	# contact map is changed by reference, so no return value
	st, sc, scaffold_circular = contactmap

	for staple_id in st :
		st[staple_id]["sequence53"] = dna.to_material(st[staple_id]["sequence53"], "RNA")

	for i53 in sc :
		sc[i53]["base_letter"] = dna.to_material(sc[i53]["base_letter"], "RNA")




def polyT_all_unhybridised_unknown(contactmap) :
	"""All unhybrised scaffold and staple bases, which have an unknown base letter (like ?, X, N)
		are changed to polyT"""

	# contact map is changed by reference, so no return value
	st, sc, _ = contactmap

	st_subs = 0
	sc_subs = 0

	for staple_id in st :
		seq = list(st[staple_id]["sequence53"])
		bid = st[staple_id]["base_i53"]

		for i, i53 in enumerate(bid) :
			if i53 == -1 and (seq[i] not in {"A", "C", "G", "T", "U"}) :
				seq[i] = "T"
				st_subs += 1

		st[staple_id]["sequence53"] = "".join(seq)

	for i53 in sc :
		if sc[i53]["staple_id"] == -1 and (sc[i53]["base_letter"] not in {"A", "C", "G", "T", "U"}) :
			sc[i53]["base_letter"] = "T"
			sc_subs += 1

	return st_subs, sc_subs




def polyT_flush(contactmap) :
	"""All scaffold and staple bases in contact map are changed to T -- complete T flush
	"""

	st, sc, _ = contactmap

	for staple_id in st :
		st[staple_id]["sequence53"] = "T" * len(st[staple_id]["sequence53"])

	for i53 in sc :
		sc[i53]["base_letter"] = "T"




def contactmap_sequence_check(contactmap, staple_material, scaffold_material) :

	st, sc, _ = contactmap

	staple_material_correct = True
	for staple_id in st :
		staple53 = st[staple_id]["sequence53"]

		if not dna.is_valid_sequence(staple53, staple_material) :
			staple_material_correct = False
			break
	
	scaffold53 = "".join([sc[i53]["base_letter"] for i53 in sorted(sc.keys())])
	unique_bases = set(scaffold53)

	scaffold_material_correct = dna.is_valid_base_set(unique_bases, scaffold_material)
	scaffold_seq_heterogeneous = len(unique_bases) > 1

	return staple_material_correct, scaffold_material_correct, scaffold_seq_heterogeneous























#
#
# ORIGAMI STATS
#
#


def scaffold_and_staple_stats(contactmap) :
	"""Returns stats about the scaffold and staples
	"""

	stats = {}

	st, sc, scaffold_circular = contactmap

	# total scaffold bases not hybridised in design
	# and GC content of scaffold
	GC = 0
	GCfraction = 0
	sc_nh = 0
	for i53 in sc :
		if sc[i53]["staple_id"] == -1 :
			sc_nh += 1
		if sc[i53]["base_letter"] in ["G", "g", "C", "c"] :
			GC += 1

	GCfraction = GC / len(sc)

	# total staple bases not hybridised in design
	st_nh = 0
	for staple_id in st :
		st_nh += st[staple_id]["base_i53"].count(-1)

	# other staple stats
	avlen = 0 				# average, min and max staple length
	minlen = 1000000
	maxlen = 0
	maxsec = 0 				# maximum NUMBER of sections a staple has

	splits = {1: 0, 2: 0, 3: 0, 4: 0} 	# number of STAPLES with 1 split, 2 splits, 3 splits, 3+ splits

	split_size_distribution = {} 		# split size (nt): frequency   -- over all staples
	total_splits = 0 					# total staple sections hybridising the scaffold

	with_loopouts = [] 	# id's of staples with loopout crossovers
	with_dangles = []	# id's of staples with staple single stranded overhangs

	staples_le5nt_sections = set() 	# id's of staples with 1 or more hybridising section of 5nt or less
	staples_le3nt_sections = set()

	for staple_id in st :

		# min/max/average staple length
		slen = len(st[staple_id]["sequence53"])
		if slen < minlen :
			minlen = slen
		if slen > maxlen :
			maxlen = slen
		avlen += slen

		# min/max staple sections
		nsec = st[staple_id]["num_sections"]
		if nsec > maxsec :
			maxsec = nsec

		# count staples in each 1,2,3,3+ split category
		if nsec < 4 :
			splits[nsec] += 1
		else :
			splits[4] += 1

		# add this staple split sizes to split distribution
		# IF the staple split HYBRIDISES with the scaffold
		for i, section_length in enumerate( st[staple_id]["section_lengths"] ) :
			if st[staple_id]["section_scaffold_domain_ids"][i] != -1 :
				# i.e. staple section does hyb scaffold
				if section_length in split_size_distribution :
					split_size_distribution[section_length] += 1
				else :
					split_size_distribution[section_length] = 1
				total_splits += 1

				# is this staple a revnano 'problem staple' i.e. has scaffold binding sections
				# of 5nt or less?
				if section_length <= 5 :
					staples_le5nt_sections.add(staple_id)
				if section_length <= 3 :
					staples_le3nt_sections.add(staple_id)

		# count staples with dangles or loopouts
		sid = st[staple_id]["section_scaffold_domain_ids"]
		if (sid[0] == -1) or (sid[-1] == -1) :
			with_dangles.append(staple_id)

		if len(sid) > 2 and (-1 in sid[1:-1]) :
			with_loopouts.append(staple_id)

	avlen = avlen / len(st)

	circ = "LINEAR"
	if scaffold_circular : circ = "CIRCULAR"


	stats["scaffold type"] = circ
	stats["scaffold length"] = len(sc)
	stats["scaffold GC fraction"] = GCfraction
	stats["scaffold bases single stranded"] = sc_nh
	stats["scaffold-staple base pairs"] = len(sc)-sc_nh
	
	stats["staple count"] = len(st)
	stats["staples with section"] = splits 	# number of staples with 1, 2, 3, 3+ sections
	stats["staple max sections"] = maxsec
	stats["staple section hybridising scaffold count"] = total_splits
	stats["staple section length distribution"] = split_size_distribution
	stats["staple min length"] = minlen
	stats["staple max length"] = maxlen
	stats["staple mean length"] = avlen
	stats["staple bases single stranded"] = st_nh
	stats["staples with single stranded dangles count"] = len(with_dangles)
	stats["staple ids with single stranded dangles"] = with_dangles
	stats["staples with single stranded loopouts count"] = len(with_loopouts)
	stats["staple ids with single stranded loopouts"] = with_loopouts

	stats["num staples with one or more hyb section <= 5nt"] = len(staples_le5nt_sections)
	stats["num staples with one or more hyb section <= 3nt"] = len(staples_le3nt_sections)

	return stats



def display_scaffold_and_staple_stats(contactmap) :

	stats = scaffold_and_staple_stats(contactmap)

	print("")
	print("- %s scaffold" % stats["scaffold type"])

	print("")
	print("- %d\t\tScaffold length (nt)" % stats["scaffold length"])
	print("- %1.2f\t\tScaffold GC fraction" % stats["scaffold GC fraction"])

	print("- %d\t\tScaffold bases single-stranded, not hybridised to staples" % stats["scaffold bases single stranded"])
	print("- %d\t\tScaffold-staple base pairs" % stats["scaffold-staple base pairs"])
	
	print("")
	print("- %d\t\tStaple count" % stats["staple count"])
	print("- %d/%d/%d/%d\t1/2/3/3+ section staples (breakdown of staple count)" % (stats["staples with section"][1], \
									stats["staples with section"][2],stats["staples with section"][3],stats["staples with section"][4]))
	print("- %d\t\tMost sections on a staple" % stats["staple max sections"])

	print("\n:::::: Staple Split Distribution :::::: (splits which hyb scaffold, only) ")
	print("Size (nt)\tPercent Fraction")
	for section_length in sorted(stats["staple section length distribution"].keys()) :
		print("%d\t\t%3.3f" % (section_length, (stats["staple section length distribution"][section_length]/stats["staple section hybridising scaffold count"])*100 ))
	print(":::::::::::::::::::::::::::::::::::::::\n")

	print("- %d --> %d\tStaple length range (nt --> nt)" % (stats["staple min length"], stats["staple max length"]))
	print("- %2.2f\t\tAverage staple length (nt)" % stats["staple mean length"])	
	print("- %d\t\tStaple bases not hybridised to scaffold (in overhangs/loopouts)" % stats["staple bases single stranded"])

	if stats["num staples with one or more hyb section <= 5nt"] > 0 :
		print("- %d\t\tStaples with hybridising sections of 5nt or less" % stats["num staples with one or more hyb section <= 5nt"])

	if stats["num staples with one or more hyb section <= 3nt"] > 0 :
		print("- %d\t\tStaples with hybridising sections of 3nt or less" % stats["num staples with one or more hyb section <= 3nt"])

	if stats["staples with single stranded dangles count"] > 0 or stats["staples with single stranded loopouts count"] > 0 :
		print("- %d\t\tStaples with dangles: %r" % (stats["staples with single stranded dangles count"], stats["staple ids with single stranded dangles"]))
		print("- %d\t\tStaples with loopout crossovers: %r" % (stats["staples with single stranded loopouts count"], stats["staple ids with single stranded loopouts"]))
		
	print("")



def repeat_distribution(scaffold53, scaffold_type = "LINEAR") :
	"""Returns dictionary: kmer length: repetitions past first occurence found

	Circular scaffolds are properly accounted for.

	Note: DBS sequences of order k don't have repeats more than k bases long
	and also will have a characteristic distribution of repeats of length less than k
	"""

	MAX_KMER_SIZE = 25

	U = {} 	# kmer length: number of times this sized kmer was found unique (before being repeated)
	R = {}  # kmer length: number of times this sized kmer was repeated, after the first find (the most interesting data)
	kmers_found = set()

	scaffold_len = len(scaffold53)

	# slide window of each kmer size over seq
	# if kmer already encountered, repeats = repeats + 1
	# add kmer to set 

	# 0101010 		kmer size 3
	# 010
	#  101 
	#   010  		1 rep
	#    101    	2 rep
	#     010   	3 rep = 3 reps of kmer size 3 in total (for linear sequence)
	#      100
	#		001     ...also 3 reps of kmer size 3 when sequence circular

	# repetitions = 0 when all subsequences are distinct

	for kmer_size in range(1,MAX_KMER_SIZE+1) :

		#print("Counting size %d repeats " % kmer_size)
		unique = 0
		repeats = 0

		# where kmer sliding window starts and ends
		i53r_start = 0
		i53r_end = scaffold_len - kmer_size
		if scaffold_type == "CIRCULAR" :
			i53r_end = scaffold_len - 1

		# sliding window of kmer
		for i53r in range(i53r_start, i53r_end+1):	# from i53r_start to i53r_end inclusive
			
			# make the next kmer
			kmer = ""
			i53rr = i53r
			while len(kmer) < kmer_size :
				i53 = dna.linear_circular(scaffold_type, scaffold_len, i53rr)
				kmer += scaffold53[i53]
				i53rr += 1

			# see if its a repeat
			if kmer in kmers_found :
				repeats += 1
			else :
				unique += 1

			kmers_found.add(kmer)

		U[kmer_size] = unique
		R[kmer_size] = repeats


	# table of repeats data
	print("kmer size k\t\tRepeats")
	for kmer_size in sorted(R.keys()) :
		print("%d\t\t\t%d" % (kmer_size, R[kmer_size]))
	print("")
	print("Repeats = All kmer's of size k, collectively taken together, were found to have this many repeats in total. A repeat is an occurence after the first occurence.")
	print("")
	print("The DBS order is the lowest kmer size for which repeats = 0")
	print("")

	# graph of repeats data
	plt.figure()
	
	kmer_size = sorted(R.keys())
	repeats = [R[k] for k in kmer_size]
	if 0 in repeats :
		DBS_order = str( kmer_size[repeats.index(0)] )
	else :
		DBS_order = "> %d" % (MAX_KMER_SIZE)

	plt.subplot(1,2,1)
	plt.plot(kmer_size, repeats, color="xkcd:black")
	plt.scatter(kmer_size, repeats, color="xkcd:black") #, marker=".", markersize=10)
	plt.xlabel("k-mer size (nt)")
	plt.ylabel("Times repeated after first occurrence")
	plt.title("")

	plt.subplot(1,2,2)
	plt.semilogy(kmer_size, repeats, marker=".", color="xkcd:black")
	plt.xlabel("k-mer size (nt)")
	plt.title("LOG SCALE")

	plt.suptitle("Scaffold Sequence k-mer Repeats (Scaffold is DBS order %s)" % DBS_order)
	plt.show()























#
#
# VERIFY CONTACT MAP
#
#

def verify_contact_map(contactmap) :

	st, sc, scaffold_circular = contactmap

	# ------------------------------------------------------------------
	# Part 1: Validity tests which MUST all be passed
	# ------------------------------------------------------------------

	# 0. Verify that all elements of st agree in size
	verify0 = True 	
	for staple_id in st :

		a = len(st[staple_id]["sequence53"])
		b = len(st[staple_id]["base_i53"])
		c = sum(st[staple_id]["section_lengths"])
		d = len(st[staple_id]["section_lengths"])
		e = len(st[staple_id]["section_scaffold_domain_ids"])
		f = st[staple_id]["num_sections"]

		# sequence length same as base list length, same as sum section lengths
		# section lengths list same size as domain id list, same as num_sections
		if not ((a == b) and (a == c) and (d == e) and (d == f)) :
			verify0 = False
			break

	# 1. Verify that no two staple bases hybridise to the same scaffold base
	i53_paired = set()

	verify1 = True
	for staple_id in st :
		
		base_i53 = st[staple_id]["base_i53"]

		for i53 in base_i53 :
			if i53 != -1 :
				if i53 in i53_paired :
					verify1 = False
					break
				else :
					i53_paired.add(i53)

	# 2. Verify that the two parts of the contact map, st and sc, are CONSISTENT (i.e. they express the same connectivity information toward each other)
	contactmap2 = ( copy.deepcopy(st), copy.deepcopy(sc), scaffold_circular )
	_reset_sc_from_st(contactmap2)

	_, sc2, _ = contactmap2
		
	verify2 = True
	for i53 in sc :
		for key in sc[i53] :
			if sc[i53][key] != sc2[i53][key] :
				verify2 = False
				break

	# 3. Verify scaffold-staple SEQUENCE COMPLEMENTARITY, by generating staple sequences from sc
	verify3 = True
	for i53 in sc :
		scaffold_base = sc[i53]["base_letter"]

		staple_id = sc[i53]["staple_id"]
		staple_base_id = sc[i53]["staple_base_id"]

		if staple_id != -1 :		# skip over scaffold bases not hybridised to a staple
			staple_base = st[staple_id]["sequence53"][staple_base_id]

			if not dna.bases_are_WC_compliments(scaffold_base, staple_base) :
				verify3 = False
				break
	
	print("----------- Verification Results -------------")
	print("%r\tTest 0: information in st part is of correct size" % verify0)
	print("%r\tTest 1: no two staple bases hybridise to the same scaffold base" % verify1)
	print("%r\tTest 2: st and sc halves of contact map are consistent" % verify2)
	print("%r\tTest 3: staple-scaffold watson-crick complementarity exists" % verify3)

	# ------------------------------------------------------------------
	# Part 2: Warnings, which are advice and are not critical
	# ------------------------------------------------------------------

	sc_bases = set()
	for i53 in sc :
		sc_bases.add( sc[i53]["base_letter"] )

	st_bases = set()
	for staple_id in st :
		st_bases.update( list(st[staple_id]["sequence53"]) )

	small4 = {}			# staple index : ([section len list], [scaffold domain id list])
	small1 = {}

	for staple_id in st :
		# record staples with small sections hybridised to scaffold
		lsec = st[staple_id]["section_lengths"]
		sid = st[staple_id]["section_scaffold_domain_ids"]
		for i, l in enumerate(lsec) :
			if (l > 1 and l < 5) and (sid[i] != -1) :
				small4[staple_id] = (lsec, sid) 		# its a very small staple section, bound to the scaffold (possibly an error)
			if (l == 1) and (sid[i] != -1) :
				small1[staple_id] = (lsec, sid)		# its a 1nt staple section bound to the scaffold -- nearly impossible

	# warning message if scaffold homogeneous sequence
	if len(sc_bases) == 1 :
		print("\n*Warning* Scaffold is HOMOGENEOUS sequence of '%s'" % (list(sc_bases)[0]))

	# warning message if scaffold or staples contain unknown base
	diff = sc_bases.difference({"A","C","G","T","U"})
	if len(diff) > 0 :
		print("\n*Warning* Scaffold contains unknown bases: %r" % (diff))
	diff = st_bases.difference({"A","C","G","T","U"})
	if len(diff) > 0 :
		print("\n*Warning* Staples contain unknown bases: %r" % (diff))

	# warning message if some HYBRIDISING staple sections are very small (and list staples)
	if len(small4) > 0 :
		print("\n*Warning* Some staples have very short (1nt< length <5nt) sections which hybridise the scaffold. This may indicate an error in the contact map.")
		print("\nStaple\tSection Lengths (nt)\tHybridises to scaffold domain ids")
		for staple_id in small4 :
			print("%d\t%r\t\t%r" % (staple_id, small4[staple_id][0], small4[staple_id][1]))

	# further warning message if some HYBRIDISING staple sections are 1nt
	if len(small1) > 0 :
		print("\n*Warning* Some staples apparently have 1nt sections which hybridise the scaffold. This may indicate an error in the contact map.")
		for staple_id in small1 :
			print("%d\t%r\t\t%r" % (staple_id, small1[staple_id][0], small1[staple_id][1]))

	print("----------------------------------------------")

	return verify1 and verify2 and verify3

	


























#
#
# EQUIVALENCE AND HAMMING DISTANCE OF CONTACT MAPS
#
#



def renumber_staples_by_sequence_matching(master_contactmap, contactmap) :
	""" This function simply makes sure that a staple in contactmap which has the same sequence as a staple in master_contactmap
	is given the same staple ID number as the staple has in master_contactmap

	(The same sequence staple may be routed differently in the two contact maps, this does not matter)

	Conditions:
	- All staple sequences in master_contactmap must be unique
	- All staple sequences in contactmap must be unique

	Returns the modified contactmap, or None if there was an error
	"""

	st0, sc0, circular0 = master_contactmap
	st1, sc1, circular1 = contactmap
	
	# all staple sequences in master_contactmap must be unique
	seq0 = set()
	for staple_id in st0 :
		seq0.add(st0[staple_id]["sequence53"])

	if len(seq0) != len(st0) :
		return None

	# all staple sequences in contactmap must be unique
	seq1 = set()
	for staple_id in st1 :
		seq1.add(st1[staple_id]["sequence53"])

	if len(seq1) != len(st1) :
		return None

	# find staple partners
	HIGH_STAPLE_ID = 10000
	st2 = {}

	for staple_id1 in st1 :
		staple1 = st1[staple_id1]

		# find sequence match partner in ground truth contact map
		partner = False
		for staple_id0 in st0 :
			staple0 = st0[staple_id0]
			if staple1["sequence53"] == staple0["sequence53"] :
				# partner found
				st2[staple_id0] = staple1  # store staple 1 with updated id
				partner = True
				break

		# if no partner in ground truth contact map, assign it a high ID number
		if partner == False :
			st2[HIGH_STAPLE_ID] = staple1
			HIGH_STAPLE_ID += 1

	contactmap2 = (st2, copy.deepcopy(sc1), circular1)

	_reset_sc_from_st(contactmap2)

	return contactmap2




def check_equivalence(origami_name1, contactmap1, origami_name2, contactmap2) :
	"""Returns True if two contact maps are equivalent

	Equivalence means: 
		- both have same scaffold type LINEAR/CIRCULAR
		- both have same listed 5' to 3' scaffold sequence
		- both have same scaffold circularity
		- both have same st part, although ordering of staples (staple_ids) is allowed to differ

	Useful for a) verifying scadnano2contact is equivalent to oxdna2contact for scadnano designs
	and b) verifying that revnano outputs a contact map equivalent to the maximal 3' contact map
	"""

	st1, sc1, scaffold_circular1 = contactmap1
	st2, sc2, scaffold_circular2 = contactmap2

	if scaffold_circular1 != scaffold_circular2 :
		print("NOT equivalent.")
		print("One scaffold is LINEAR, one is CIRCULAR.")
		return False

	seq53_1 = "".join([sc1[i53]["base_letter"] for i53 in sorted(sc1.keys())])
	seq53_2 = "".join([sc2[i53]["base_letter"] for i53 in sorted(sc2.keys())])

	if seq53_1 != seq53_2 :
		print("NOT equivalent:")
		print("Scaffolds have different sequences and/or lengths when listed 5' to 3'.")	
		return False	

	# compare staple connections and sequences 
	# (staple id's are not compared, these don't have to match, two designs can be identical with differently labelled staples)
	set1 = set(); set2 = set()
	for staple_id in st1 : set1.add(_make_staple_identifier(contactmap1, staple_id))
	for staple_id in st2 : set2.add(_make_staple_identifier(contactmap2, staple_id))

	if set1 != set2 :
		print("NOT equivalent:")
		print("Staple information is not identical in both contact maps.")
		print("Saving file _assets/diff.csv detailing the differences.\n")

		# make CSV file of the difference
		# staples in orgiami1, not in origami2
		diff1 = _make_diff(contactmap1, set1, set2)
		# staples in origami2, not in origami1
		diff2 = _make_diff(contactmap2, set2, set1)
		iohandler.write_contactmap_diff(origami_name1, diff1, origami_name2, diff2, "_assets/diff.csv")

		return False			

	print("** EQUIVALENT **\n")
	return True






def hamming_distance(contactmap1, contactmap2) :
	""" Returns the hamming distance: the number of scaffold base targets which need to be changed in contactmap2 to get to contactmap1

	contactmap1 = original contact map
	contactmap2 = contact map to compare to

	Returns None if the scaffolds are different sequences, and/or different circularity
	(cannot return False, as this is same as returning integer 0 for no hamming distance)

	Calculation simply goes through the scaffold strand and counts how many staple base targets 
	are different on contactmap2, as compared to contactmap1

	Important: This means that equivalent staples in both contact maps are ASSUMED to have EQUAL staple id numbers

	The maximum hamming distance is the scaffold length
	"""

	st1, sc1, scaffold_circular1 = contactmap1
	st2, sc2, scaffold_circular2 = contactmap2

	if scaffold_circular1 != scaffold_circular2 :
		print("One scaffold is LINEAR, one is CIRCULAR.")
		return None

	seq53_1 = "".join([sc1[i53]["base_letter"] for i53 in sorted(sc1.keys())])
	seq53_2 = "".join([sc2[i53]["base_letter"] for i53 in sorted(sc2.keys())])

	if seq53_1 != seq53_2 :
		print("Scaffolds have different sequences and/or lengths when listed 5' to 3'.")	
		return None

	# calc hamming distance
	hamming = 0

	for i53 in sc2 :	# iteration order not important

		# does this scaffold base hybridise with the same staple base as in the master contact map?
		same_staple = (sc2[i53]["staple_id"] == sc1[i53]["staple_id"])
		same_staple_base = (sc2[i53]["staple_base_id"] == sc1[i53]["staple_base_id"])

		if not (same_staple and same_staple_base) :
			hamming += 1

	return hamming




























#
#
# DISPLAY CONTACT MAP AS BASE PAIR MATRIX
#
#


def display_bp_matrix(origami_name, contactmap) :
	"""Note: only for small origamis, e.g. less than 1000nt scaffold. Or else you cannot see data in the matrix.
	"""

	st, sc, scaffold_circular = contactmap

	if len(sc) > 1000 :
		print("Error: Origami too large to display contact map as matrix.\n")
		return 

	VALUE_NOT_PAIRED = 0
	VALUE_NEVER_PAIRED = 0.5
	VALUE_PAIRED = 1
	COLOUR_MAP = "binary"

	# -----------------------
	# 1. Data Matrix
	# -----------------------

	# note: rows (y-axis) of data matrix are scaffold bases
	#		columns (x-axis) of matrix are staple bases

	# x-axis of matrix
	where_is_on_x = {} 	# staple_id, staple_base_id : index on x-axis  (first index is 0)

	x = 0
	x_ticks = []		# want x-ticks only when a new staple starts
	standard_row = []

	for staple_id in sorted(st.keys()) :

		x_ticks.append(x)

		for staple_base_id in range(0, len(st[staple_id]["sequence53"])) :
			
			base_i53 = st[staple_id]["base_i53"][staple_base_id]
			if base_i53 == -1 :
				standard_row.append(VALUE_NEVER_PAIRED)
			else :
				standard_row.append(VALUE_NOT_PAIRED)

			where_is_on_x[(staple_id, staple_base_id)] = x	# used later, when building connection matrix

			x += 1

	x_tick_labels = [str(staple_id) for staple_id in sorted(st.keys())]

	# y-axis of matrix
	y_ticks = [0, len(sc)-1]
	y_tick_labels = [str(y) for y in y_ticks]

	# fill in each row of data matrix
	datamatrix = []

	for i53 in sorted(sc.keys()) :
		staple_id = sc[i53]["staple_id"]
		staple_base_id = sc[i53]["staple_base_id"]
		if staple_id == -1 :
			datamatrix.append( [VALUE_NEVER_PAIRED] * len(standard_row) )
		else :
			# add the standard row, where staple bases that are never paired are marked
			datamatrix.append( copy.deepcopy(standard_row) )
			# change the single staple base that IS paired to this scaffold base
			datamatrix[i53][ where_is_on_x[(staple_id, staple_base_id)] ] = VALUE_PAIRED
	
	# -------------------------
	# 2. Heatmap of data matrix
	# -------------------------

	plt.figure()
	ax = plt.subplot(1,1,1)
	
	# plot
	plt.imshow(datamatrix, interpolation='nearest', cmap=COLOUR_MAP)

	# axis tick labels
	plt.xticks(x_ticks)
	ax.set_xticklabels(x_tick_labels)
	plt.yticks(y_ticks)
	ax.set_yticklabels(y_tick_labels)
	ax.tick_params(labelsize = 3.0) # font size of tick labels

	# axis labels and title
	plt.xlabel("Staples 5'--->3'")
	plt.ylabel("Scaffold 3'<---5'")
	plt.title("%s Contact Map" % (origami_name))

	plt.show()
	#plt.savefig("cmap.png", bbox_inches='tight', dpi=300)
	plt.close()








































#
#
# CHECK STAPLE SEQUENCES MATCH A SUPPLIED CSV LIST
#
#

def staples_match(st, staples_in) :


	# note: the ORDER of staples in both collections is NOT SIGNIFICANT
	# only 1) the same number of staples must exist and 2) the same sequences must exist

	only_in_st = []

	for staple_id in st :
		seq = st[staple_id]["sequence53"].upper()

		try :
			# remove (first occurence of) staple sequence from CSV list
			staples_in.remove(seq)
		except(ValueError) :
			# staple in contact map is not in the remaining list, and cannot be removed from remaining list
			only_in_st.append( seq )

	if len(staples_in) == 0 :
		return "All match!"
	
	if len(only_in_st) > 0 or len(staples_in) > 0 :
		s1 = "Only in contact map: %r" % only_in_st
		s2 = "Only in file: %r" % staples_in

		return "MISMATCH: %s; %s" % (s1, s2)

































#
#
# EXECUTE
#
#

if __name__ == "__main__" :

	print("\n---------------------")
	print("Contact Map Utilities")
	print("---------------------\n")
	
	if len(sys.argv) < 3 :
		print("Common usage: python3 contactutils <command> <origami name>\n")
		quit()

	cmd = sys.argv[1]
	origami_name = sys.argv[2]

	# ---------------------------------------------------------------------------------
	# ** special case not involving contact map: conversion of a revnano file
	# ** usage: python3 contactutils.py rev2revstar myfile.rev
	if cmd == "rev2revstar" :

		# revnano input file, to an auto-starred version of revnano input file, with all polyT islands on staples marked
		# - assumes revnano input file is VALID
		# - RE-WRITES the inputted REVNANO file

		print("Converting REVNANO input file, to a version where all polyT staple regions are * delimited")
		iohandler.REVNANO_auto_star("_assets/%s" % (origami_name))
		print("[Done]\n")
		quit()
	# ---------------------------------------------------------------------------------


	if not os.path.exists("_assets/%s.csv" % (origami_name)) :
		print("Error: contact map file '%s.csv' does not exist in the _assets/ directory.\n" % (origami_name))
		quit()		

	contactmap = iohandler.read_origami_contactmap("_assets/%s.csv" % (origami_name))		

	if isinstance(contactmap, str) :
		print("Error: %s\n" % (contactmap))
		quit()		

	print("Contact map: %s\n" % (origami_name))


	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py checkstaples origami_name staples.csv
	# ---------------------------------------------------------------------------------
	if cmd == "checkstaples" :

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py checkstaples <origami_name> <staples csv file>\n")
			quit()

		staples_filename = "_assets/%s" % (sys.argv[3])

		if not os.path.exists(staples_filename) :
			print("Error: staples CSV file '%s' does not exist.\n" % (staples_filename))
			quit()

		
		staples_in = iohandler.read_staples_csv(staples_filename)

		st, _, _ = contactmap

		print("Checking staple sequences in contact map match those in CSV file...")

		print( staples_match(st, staples_in) )

		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py scaffold origami_name seq.txt dna dna
	# ---------------------------------------------------------------------------------
	elif cmd == "scaffold" :

		if len(sys.argv) != 6 :
			print("Usage: python3 contactutils.py scaffold <origami_name> <sequence file> <scaffold material> <staple material>\n")
			quit()


		seq_filename = "_assets/%s" % (sys.argv[3])

		if not os.path.exists(seq_filename) :
			print("Error: scaffold sequence file '%s' does not exist.\n" % (seq_filename))
			quit()

		# new scaffold sequence has to fit
		scaffold53_new = iohandler.read_sequence(seq_filename)

		_, sc, scaffold_circular = contactmap
		# if scaffold is circular, an exact length sequence is needed
		if scaffold_circular :
			if len(scaffold53_new) != len(sc) :					
				print("Error. The new scaffold sequence must be exactly %d nt. The current scaffold is circular.\n" % (len(sc)))
				quit()
		else :
			# scaffold linear
			if len(scaffold53_new) < len(sc) :		
				print("Error. The new scaffold sequence must be %d nt or longer.\n" % (len(sc)))
				quit()

		scaffold_material = sys.argv[4].upper()
		staple_material = sys.argv[5].upper()
		
		if (scaffold_material not in ["DNA", "RNA"]) or (staple_material not in ["DNA", "RNA"]) :
			print("Error. Scaffold and staple materials must be specified as 'DNA' or 'RNA'.\n")
			quit()

		if not dna.is_valid_sequence(scaffold53_new, scaffold_material) :
			print("Error. New scaffold sequence supplied is not pure %s.\n" % (scaffold_material))
			quit()
				
		print("Applying new scaffold sequence...")

		change_scaffold_sequence(contactmap, scaffold_material, staple_material, scaffold53_new)

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")		
		print("[Done]\n")
		
	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py circ origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "circ" :

		print("Making scaffold circular...")

		success = make_scaffold_circular(contactmap)

		if not success :
			print("Error. This linear scaffold could not be circularised. It could already be circular. Or, a staple may bridge the scaffold ends with a >1nt loopout crossover (not permitted).\n")
			quit()

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py linear origami_name 20
	# ---------------------------------------------------------------------------------
	elif cmd == "linear" :

		print("Making scaffold linear...")

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py linear origami_name <scaffold base number of 3 prime>\n")
			quit()

		try :
			i53_scaffold_3prime = int(sys.argv[3])
		except ValueError:
			# cannot convert value to an int
			print("Error. Supply an integer for the base number of the scaffold 3' end.\n")
			quit()

		_, sc, _ = contactmap

		if (i53_scaffold_3prime) < 0 or (i53_scaffold_3prime > len(sc)-1) :
			print("Error. Base number of the scaffold 3' end is off the beginning or end of the scaffold.\n")
			quit()
			
		success = make_scaffold_linear(contactmap, i53_scaffold_3prime)

		if not success :
			print("Error. Could not make linear.\n")
			quit()

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")
		print("[Done]\n")		

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py thermo origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "thermo" :

		_, _, scaffold_circular = contactmap
		R = get_thermodynamically_valid_scaffold_rotation_offsets(contactmap)

		if scaffold_circular :
			print("Scaffold is CIRCULAR\n")
			print("All base indexes are thermodynamically valid rotation offsets.")
		else :
			print("Scaffold is LINEAR\n")
			print("The following rotation offsets are thermodynamically valid. None lead to staple sections smaller than %d base pairs. The rotation offset of the current design is also included.\n" % (MIN_STAPLE_SECTION_LENGTH))
			print(R)
			if len(R) == 1 :
				print("\nNote: This scaffold cannot be rotated. The scaffold has only one valid rotation offset - its current position.")
		
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py rotate origami_name 34 dna dna
	# ---------------------------------------------------------------------------------
	elif cmd == "rotate" :

		if len(sys.argv) != 6 :
			print("Usage: python3 contactutils.py rotate <origami_name> <rotation offset> <scaffold material> <staple material>\n")
			quit()

		try :
			i53_scaffold_3prime = int(sys.argv[3])
		except ValueError:
			# cannot convert value to an int
			print("Error. Supply an integer for the rotation offset.\n")
			quit()

		_, sc, _ = contactmap

		if (i53_scaffold_3prime) < 0 or (i53_scaffold_3prime > len(sc)-1) :
			print("Error. Rotation offset is off the beginning or end of the scaffold.\n")
			quit()

		# check material supplied is consistent with sequences in contactmap
		scaffold_material = sys.argv[4].upper()
		staple_material = sys.argv[5].upper()
		staple_material_correct, scaffold_material_correct, _ = contactmap_sequence_check(contactmap, staple_material, scaffold_material)

		if not scaffold_material_correct :
			print("Error: Scaffold sequence in contactmap is not pure %s.\n" % (scaffold_material))
			quit()

		if not staple_material_correct :
			print("Error: Staple sequences in contactmap are not all pure %s.\n" % (staple_material))
			quit()

		print("Rotating origami scaffold to offset %d...." % (i53_scaffold_3prime))

		success = apply_scaffold_rotation_offset(contactmap, i53_scaffold_3prime, scaffold_material, staple_material)

		if not success :
			print("Error. Applying this rotation offset failed.\n")
			quit()

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py del 12,34,32,10
	# ---------------------------------------------------------------------------------
	elif cmd == "del" :

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py del origami_name <staple id list>\n")
			quit()

		try :
			data = sys.argv[3].split(",")
			data = [int(d) for d in data]
		except ValueError:
			# cannot convert value to an int
			print("Error: The list of staple ids should be supplied as a comma separated list with no spaces, e.g. 23,12,244\n")
			quit()

		if min(data) < 0 :
			print("Error: Staple ids cannot be negative.\n")
			quit()

		print("Deleting staples %r...." % (data))

		success = delete_staples(contactmap, data)

		if not success :
			print("Error: Some staples specified do not exist in the origami contact map.\n")
			quit()

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py move origami_name 169,-1
	# ---------------------------------------------------------------------------------
	elif cmd == "move" :

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py move origami_name <scaffold base number,change>\n")
			quit()

		error = False

		try :
			data = sys.argv[3].split(",")
			data = [int(d) for d in data]
		except ValueError:
			# cannot convert value to an int
			error = True

		if len(data) != 2 or error:
			print("Error: Scaffold base to move should be supplied as two numbers separated by comma (no spaces): <scaffold base number,change>\n")
			quit()

		i53, change = data

		if i53 < 0 :
			print("Error: Scaffold base cannot be negative.\n")
			quit()

		print("Moving scaffold base %d by %d...." % (i53, change))

		success = move_scaffold_base(contactmap, i53, change)

		if not success :
			print("Error: This move is invalid.\n")
			quit()

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")	
		print("[Done]\n")		

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py rna origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "rna" :
		
		print("Converting to RNA...")
		
		change_DNA_to_RNA(contactmap)

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")	
		print("[Done]\n")	
		
	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py tlost origami_name	
	# ---------------------------------------------------------------------------------
	elif cmd == "tlost" :

		print("Converting unhybrised and unknown bases on scaffold and staples to polyT...")
		
		st_subs, sc_subs = polyT_all_unhybridised_unknown(contactmap)

		print("%d bases change in staples" % (st_subs))
		print("%d bases change on scaffold" % (sc_subs))

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")		
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py tflush origami_name	
	# ---------------------------------------------------------------------------------
	elif cmd == "tflush" :

		print("Converting ALL scaffold and staple bases to T...")

		polyT_flush(contactmap)

		# RE-WRITE the current contact map
		iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py stats origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "stats" :

		print("Origami statistics:")

		display_scaffold_and_staple_stats(contactmap)

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py repeats origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "repeats" :
		
		_, sc, scaffold_circular = contactmap
		scaffold53 = "".join([sc[i53]["base_letter"] for i53 in sorted(sc.keys())])
		scaffold_type = "LINEAR"
		if scaffold_circular : scaffold_type = "CIRCULAR"

		print("Scaffold repeats distribution:")

		repeat_distribution(scaffold53, scaffold_type = scaffold_type)

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py verify origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "verify" :

		verify_contact_map(contactmap)

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py diff origami_name1 origami_name2
	# ---------------------------------------------------------------------------------
	elif cmd == "diff" :

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py diff <origami_name1> <origami_name2>\n")
			quit()

		origami_name2 = sys.argv[3]
		contactmap2 = iohandler.read_origami_contactmap("_assets/%s.csv" % (origami_name2))		

		if isinstance(contactmap2, str) :
			print("Error: %s\n" % (contactmap2))
			quit()		

		print("Contact map: %s" % (origami_name2))
		print("\nChecking if equivalent...\n")

		check_equivalence(origami_name, contactmap, origami_name2, contactmap2)

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py matrix origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "matrix" :		

		display_bp_matrix(origami_name, contactmap)

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py fasta origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "fasta" :

		print("Writing scaffold and staple sequences as FASTA file...")
		iohandler.write_FASTA(contactmap, "_assets", "%s.fasta" % (origami_name))
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py rev origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "rev" :

		# contact map to revnano input file

		print("Generating REVNANO input file from contact map...")
		iohandler.write_REVNANO(contactmap, "_assets", "%s.rev" % (origami_name))
		print("[Done]\n")

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py del1nt origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "del1nt" :

		try :
			print("Deleting 1nt scaffold domains from contact map...")
			
			contactmapD = delete_1nt_scaffold_domains(contactmap)

			iohandler.write_origami_contactmap(contactmapD, "_assets", "%sD.csv" % (origami_name), madeby="")

			print("---> modified contact map saved to: _assets/%sD.csv" % origami_name)

			print("[Done]\n")

		except ValueError as exp :
			print(exp)

	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py rs1nt origami_name
	# ---------------------------------------------------------------------------------
	elif cmd == "rs1nt" :

		try :
			print("Removing some staples from contact map so all scaffold domains larger than 1nt...")
			
			contactmapR, staple_ids_removed = remove_staples_to_remove_1nt_scaffold_domains(contactmap)

			print("The following staples were removed: %r" % staple_ids_removed)
			
			iohandler.write_origami_contactmap(contactmapR, "_assets", "%sR.csv" % (origami_name), madeby="")

			print("---> modified contact map saved to: _assets/%sR.csv" % origami_name)

			print("[Done]\n")

		except ValueError as exp :
			print(exp)


	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py hamming origami_original origami_reconstructed
	# ---------------------------------------------------------------------------------
	elif cmd == "hamming" :

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py hamming <origami_original> <origami_reconstructed>\n")
			quit()

		origami_reconstructed = sys.argv[3]
		contactmap_reconstructed = iohandler.read_origami_contactmap("_assets/%s.csv" % (origami_reconstructed))		

		if isinstance(contactmap_reconstructed, str) :
			print("Error: %s\n" % (contactmap_reconstructed))
			quit()		

		print("Calculating Hamming Distance...")

		h = hamming_distance(contactmap, contactmap_reconstructed)

		if h == None :
			print("Error. The two contact maps cannot be compared.\n")
			quit()

		print("\nH = %d base changes\n" % h)

		print("The reconstructed origami needs the hybridisation target of %d scaffold bases changing, in order to match the original origami precisely.\n" % h)


	# ---------------------------------------------------------------------------------
	# usage: python3 contactutils.py dangles origami_name staple_dangles_filename
	# ---------------------------------------------------------------------------------
	elif cmd == "dangles" :

		if len(sys.argv) != 4 :
			print("Usage: python3 contactutils.py dangles <origami_name> <staple_dangles_filename>\n")
			quit()

		staple_dangles_filename = "_assets/%s" % (sys.argv[3])

		print("Giving staples dangle extensions...")

		success = add_staple_end_dangles_from_file(staple_dangles_filename, contactmap)

		if success :
			# RE-WRITE the current contact map
			iohandler.write_origami_contactmap(contactmap, "_assets", "%s.csv" % (origami_name), madeby="")	
			print("[Done]\n")


	# ---------------------------------------------------------------------------------
	# UNRECOGNISED
	# ---------------------------------------------------------------------------------
	else :
		print("Unrecognised command '%s'.\n" % (cmd))
