"""Transforms an origami contact map into a domain level origami graph in NetworkX

Note: 1nt domains are NOT ALLOWED. All domains must be 2bp minimum.

Ben Shirt-Ediss, September-October 2021
"""

import sys
sys.path.append("../")

import networkx

from contactmap import contactutils


CROSSOVER_SCHEMATIC_LENGTH = 3

SS_DOMAIN_MAX_SCHEMATIC_LEN = 200 			# ssDNA scaffold domains over this length get CAPPED to this display length
 											# this stops large ssDNA scaffold loops from distorting the display of the domain-level graph









#
#
# SCAFFOLD GRAPH (ONLY)
#
#


def cap_ss_domain_schematic_len(nt) :

	if nt > SS_DOMAIN_MAX_SCHEMATIC_LEN :
		return SS_DOMAIN_MAX_SCHEMATIC_LEN

	return nt




def scaffold_dl_graph(contactmap) :
	"""Graph of the scaffold only. This is the initial domain level graph, before staples are added.
	Edges are ssDNA domains, or ghost nicks. The graph is either a line or circle of connected nodes.
	There are no staple crossovers yet (and none of the scaffold nicks are scaffold crossovers either).
	"""

	# Warning: if start and end i53 lists contain the same i53 number, it means that some domains are a SINGLE BASE PAIR (start and end at same i53)
	# This is not allowed in the domain-level graph representation, because the nodes
	# at each side of a domain edge need to have different i53 indexes. This means the minimum edge length is 2 bases. 

	I = contactutils.scaffold_1nt_domain_locations(contactmap)
	if len(I) > 0 :
		print("Error. Could not build the domain-level origami graph.\n1nt scaffold domains found at scaffold bases %r.\n1nt domains cannot be represented in the origami domain-level graph.\n" % I)
		print("Try to remove the 1nt domains by e.g:")
		print("python3 contactutils.py del1nt origami_name")
		print("python3 contactutils.py rs1nt origami_name\n")
		
		quit()

	_, sc, scaffold_circular = contactmap

	# 0. get scaffold sequence
	sequence53 = ""
	for i53 in sorted(sc.keys()) :
		sequence53 += sc[i53]["base_letter"]

	# 1. get i53 indexes where domains start and finish on the scaffold
	scaffold_nt = len(sc)
	i53_start_domain = []   	# bases which start and end domains are nodes on the scaffold graph
	i53_end_domain = []
	i53_all_domain = []

	for i53 in sorted(sc.keys()) :		# go from 5' to 3' base on scaffold

		if i53 == 0 :						# vertex for first base
			i53_start_domain.append(i53)
			i53_all_domain.append(i53)

		if i53 == scaffold_nt - 1 :			# vertex for last base
			i53_end_domain.append(i53)
			i53_all_domain.append(i53)
			break

		# if next base has different domain, make a node pair
		domain_id = sc[i53]["scaffold_domain_id"]
		next_domain_id = sc[i53+1]["scaffold_domain_id"]

		if domain_id != next_domain_id :
			i53_end_domain.append(i53)
			i53_start_domain.append(i53+1)
			i53_all_domain.append(i53)
			i53_all_domain.append(i53+1)

	# 2. construct scaffold graph
	G = networkx.Graph()

	# add nodes, and label them
	for i53 in i53_all_domain :
		G.add_node(i53)

	# mark physical 5' and 3' ends of scaffold (these nodes may get deleted below, if scaffold is circular)
	G.nodes[0]["type"] = "5prime"
	G.nodes[i53_all_domain[-1]]["type"] = "3prime"

	# add edges (scaffold backbone) connecting the nodes
	# all edges point in 5' to 3' direction
	for i, i53s in enumerate(i53_start_domain) :
		
		i53e = i53_end_domain[i]

		# edge representing domain
		G.add_edge(i53s, i53e)
		G.edges[i53s, i53e]["type"] = "ss_domain"
		G.edges[i53s, i53e]["scaffold_domain_id"] = sc[i53s]["scaffold_domain_id"]
		G.edges[i53s, i53e]["runlen"] = i53e - i53s + 1
		G.edges[i53s, i53e]["schematiclen"] = cap_ss_domain_schematic_len(i53e - i53s + 1 - 1)
		G.edges[i53s, i53e]["source"] = i53s 		# i53 of 5' of scaffold domain
		G.edges[i53s, i53e]["target"] = i53e 		# i53 of 3' of scaffold domain
		G.edges[i53s, i53e]["scaffold_sequence53"] = sequence53[i53s:i53e+1]
		G.edges[i53s, i53e]["base_i53"] = list(range(i53s, i53e+1))

		# ghost nick, from this domain to the next (if it exists)
		if i < len(i53_start_domain) - 1 :
			i53s_next = i53_start_domain[i+1]
			G.add_edge(i53e, i53s_next)
			G.edges[i53e, i53s_next]["type"] = "ghost_nick"
			G.edges[i53e, i53s_next]["scaffold_domain_id"] = -1
			G.edges[i53e, i53s_next]["runlen"] = 0
			G.edges[i53e, i53s_next]["schematiclen"] = 1
			G.edges[i53e, i53s_next]["source"] = i53e
			G.edges[i53e, i53s_next]["target"] = i53s_next
			
	# 3. correct graph for case of circular scaffold
	if scaffold_circular :

		last_domain_id = sc[scaffold_nt-1]["scaffold_domain_id"]	# the last domain id on the scaffold. The first is 0.

		if last_domain_id != 0 :
			# Case 1: The domains at the "virtual" 5' and 3' ends of the circular scaffold are distinct, independent domains
			# Just add a ghost nick between end scaffold domain, and start scaffold domain
			G.add_edge(scaffold_nt-1, 0)
			G.edges[scaffold_nt-1, 0]["type"] = "ghost_nick"
			G.edges[scaffold_nt-1, 0]["scaffold_domain_id"] = -1
			G.edges[scaffold_nt-1, 0]["runlen"] = 0
			G.edges[scaffold_nt-1, 0]["schematiclen"] = 1
			G.edges[scaffold_nt-1, 0]["source"] = scaffold_nt-1
			G.edges[scaffold_nt-1, 0]["target"] = 0

		else :
			# Case 2: The domains at the "virtual" 5' and 3' ends of the circular scaffold are the SAME domain
			# The first domain on the scaffold should be removed, and the last domain joined to the ghost nick
			# that followed the first domain

			# Important Note: this will not work for pathological case of just 1 domain in the graph
			# i.e. if the origami is a single helix duplex. This duplex cannot be circularised.
			if G.number_of_edges() == 1 :
				print("Error. The origami is marked as CIRCULAR scaffold in the contact map, but the origami only has 1 domain. A circular scaffold domain-level graph cannot be constructed.\n")
				quit()

			# save attributes of last scaffold domain, then delete it
			i53s = i53_start_domain[-1]
			i53e = i53_end_domain[-1]
			i53s_newedge = i53s
			runlen = G.edges[i53s, i53e]["runlen"]
			seq53 = G.edges[i53s, i53e]["scaffold_sequence53"]
			b53 = G.edges[i53s, i53e]["base_i53"]
			G.remove_edge(i53s, i53e)
			G.remove_node(i53e)

			# save attributes of first scaffold domain, then delete it
			i53s = i53_start_domain[0]
			i53e = i53_end_domain[0]
			i53e_newedge = i53e			
			runlen += G.edges[i53s, i53e]["runlen"]
			seq53 += G.edges[i53s, i53e]["scaffold_sequence53"]
			b53 += G.edges[i53s, i53e]["base_i53"]
			G.remove_edge(i53s, i53e)
			G.remove_node(i53s)

			# put in a new single domain, replacing the two deleted domains
			G.add_edge(i53s_newedge, i53e_newedge)
			G.edges[i53s_newedge, i53e_newedge]["type"] = "ss_domain"	
			G.edges[i53s_newedge, i53e_newedge]["scaffold_domain_id"] = sc[i53s_newedge]["scaffold_domain_id"]
			G.edges[i53s_newedge, i53e_newedge]["runlen"] = runlen
			G.edges[i53s_newedge, i53e_newedge]["schematiclen"] = cap_ss_domain_schematic_len(runlen - 1)
			G.edges[i53s_newedge, i53e_newedge]["source"] = i53s_newedge
			G.edges[i53s_newedge, i53e_newedge]["target"] = i53e_newedge
			G.edges[i53s_newedge, i53e_newedge]["scaffold_sequence53"] = seq53
			G.edges[i53s_newedge, i53e_newedge]["base_i53"] = b53

	return G









#
#
# FULL ORIGAMI GRAPH
#
#

def get_staple_update_info(contactmap) :
	"""Returns modifications to make to the bare scaffold graph, when all staples are hybridised
	These modifications update the scaffold graph to the fully folded origami graph
	"""

	# modifications

	single_domain_staples = {}	# (i53_5prime, i53_3prime) : staple_id

	crossovers = {}				# (i53_3prime, i53_5prime) : (type, staple_id, section_id_3prime, section_id_5prime, dom_id_3prime, dom_id_5prime, runlen)	

	staple_dangles = {}			# (i53, i53) : (type, staple_id, section_id, dom_id_joinedto, runlen)


	"""
	Notes: 
		(1) a crossover is the 3' of one staple section joined to
		the 5' of a following staple section 
			
		--hyb to scaffold-->3' ====crossover=== 5'--hyb to scaffold-->
		dom_id_3prime   i53_3prime         i53_5prime  dom_id_5prime

		(2) a loopout is the same, but the crossover is a run of bases not hybridised to the scaffold
	"""

	st, _, _ = contactmap

	# go through all staples
	for staple_id in st :

		seclen = st[staple_id]["section_lengths"]	# staple section lengths
		base_i53 = st[staple_id]["base_i53"]		# i53 indexes of all staple bases

		if len(seclen) == 1 :
			single_domain_staples[(base_i53[0], base_i53[-1])] = staple_id
			continue

		# go through all sections of current staple

		# break base_i53 into seclen chunks
		# https://www.codegrepper.com/code-examples/typescript/how+to+break+a+list+unequal+size+chunks+in+python
		it = iter(base_i53)
		sections = [[next(it) for _ in range(size)] for size in seclen]
		
		i53_3prime = None
		dom_id_3prime = None
		section_id_3prime = None
		loopout_nt = -1
		section_id_loopout = -1

		for section_id, sec in enumerate(sections) :
 
			if section_id == 0 and (-1 in sec) :
				# the beginning staple section is dangling ssDNA
				# add loopout node, and join it to next sec (which is assumed to exist, or else staple would be completely disconnected)
				sec_next = sections[section_id+1]
				i53_5prime = sec_next[0]
				dom_id_joinedto = st[staple_id]["section_scaffold_domain_ids"][section_id+1]
				
				staple_dangles[(i53_5prime + 10000000, i53_5prime)] = ("beginning_dangle", staple_id, 0, dom_id_joinedto, len(sec))

			elif section_id == (len(sections) - 1) and (-1 in sec) :
				# the ending staple section is dangling ssDNA
				# add loopout node, and join to to last sec
				sec_prev = sections[section_id-1]
				i53_3prime = sec_prev[-1]
				dom_id_joinedto = st[staple_id]["section_scaffold_domain_ids"][section_id-1]

				staple_dangles[(i53_3prime, i53_3prime + 10000000)] = ("ending_dangle", staple_id, section_id, dom_id_joinedto, len(sec))

			elif (-1 in sec) :
				# this is an internal staple loopout crossover
				loopout_nt = len(sec)
				section_id_loopout = section_id

			else :
				# this staple section is hybridised to the scaffold
				i53_5prime = sec[0]
				dom_id_5prime = st[staple_id]["section_scaffold_domain_ids"][section_id] 	# scaffold domain id section is bound to
				section_id_5prime = section_id

				# join 3 prime of last staple section hybridising the scaffold
				# (if it exists) to the 5' of this section, via a crossover or via a loopout
				if i53_3prime != None :

					if loopout_nt == -1 :
						# standard crossover
						crossovers[(i53_3prime, i53_5prime)] = ("crossover", staple_id, section_id_3prime, section_id_5prime, -1, dom_id_3prime, dom_id_5prime, 0)
					else :
						# staple loopout
						crossovers[(i53_3prime, i53_5prime)] = ("loopout", staple_id, section_id_3prime, section_id_5prime, section_id_loopout, dom_id_3prime, dom_id_5prime, loopout_nt)

				i53_3prime = sec[-1] 				# the 5' section now becomes the previous 3' section
				dom_id_3prime = dom_id_5prime
				section_id_3prime = section_id
				loopout_nt = -1

	return single_domain_staples, crossovers, staple_dangles




def get_staple_section_sequence53(contactmap, staple_id, section_id) :

	st, _, _ = contactmap

	seq53 = st[staple_id]["sequence53"]
	slen = st[staple_id]["section_lengths"]

	# chunk the sequence according to section lengths
	seq53_chunks = []; i = 0
	for s in slen :
		seq53_chunks.append( seq53[i:i+s] )
		i = i+s

	return seq53_chunks[section_id]


										 

def origami_dl_graph(contactmap) :

	G = scaffold_dl_graph(contactmap)

	single_domain_staples, crossovers, staple_dangles = get_staple_update_info(contactmap)	

	# 1. apply single domain staples
	for edge in single_domain_staples :
		
		u,v = edge
		G.edges[u, v]["type"] = "ds_domain"

		# update edge
		staple_id = single_domain_staples[edge]
		G.edges[u, v]["staple_id"] = staple_id
		G.edges[u, v]["section_id"] = 0
		G.edges[u, v]["staple_sequence53"] = get_staple_section_sequence53(contactmap, staple_id, 0)
		# update nodes
		G.nodes[u]["staple_id"] = staple_id
		G.nodes[v]["staple_id"] = staple_id

	# 2. apply crossovers and loopouts
	for crossover in crossovers :

		i53_3prime, i53_5prime = crossover

		xtype, staple_id, section_id_3prime, section_id_5prime, section_id_loopout, dom_id_3prime, dom_id_5prime, runlen = crossovers[crossover]

		# add crossover (or loopout) edge
		G.add_edge(i53_3prime, i53_5prime)	
		G.edges[i53_3prime, i53_5prime]["type"] = xtype
		G.edges[i53_3prime, i53_5prime]["scaffold_domain_id"] = -1		# not a part of the scaffold
		G.edges[i53_3prime, i53_5prime]["runlen"] = runlen

		if xtype == "crossover" :
			G.edges[i53_3prime, i53_5prime]["schematiclen"] = CROSSOVER_SCHEMATIC_LENGTH
		else : # loopout: these also have staple sequence added
			G.edges[i53_3prime, i53_5prime]["schematiclen"] = runlen + 1
			G.edges[i53_3prime, i53_5prime]["section_id"] = section_id_loopout
			G.edges[i53_3prime, i53_5prime]["staple_sequence53"] = get_staple_section_sequence53(contactmap, staple_id, section_id_loopout)

		G.edges[i53_3prime, i53_5prime]["source"] = i53_3prime 			# in direction of staple
		G.edges[i53_3prime, i53_5prime]["target"] = i53_5prime
		G.edges[i53_3prime, i53_5prime]["staple_id"] = staple_id
	
		# make scaffold domains joined by the crossover double stranded, by setting type attribute of edges
		edge_3prime = [(u,v) for u,v,a in G.edges(data=True) if a["scaffold_domain_id"] == dom_id_3prime]
		u, v = edge_3prime[0]
		G.edges[u, v]["type"] = "ds_domain"
		G.edges[u, v]["staple_id"] = staple_id
		G.edges[u, v]["section_id"] = section_id_3prime
		G.edges[u, v]["staple_sequence53"] = get_staple_section_sequence53(contactmap, staple_id, section_id_3prime)

		G.nodes[u]["staple_id"] = staple_id
		G.nodes[v]["staple_id"] = staple_id

		edge_5prime = [(u,v) for u,v,a in G.edges(data=True) if a["scaffold_domain_id"] == dom_id_5prime]
		u, v = edge_5prime[0]
		G.edges[u, v]["type"] = "ds_domain"
		G.edges[u, v]["staple_id"] = staple_id
		G.edges[u, v]["section_id"] = section_id_5prime
		G.edges[u, v]["staple_sequence53"] = get_staple_section_sequence53(contactmap, staple_id, section_id_5prime)		

		G.nodes[u]["staple_id"] = staple_id
		G.nodes[v]["staple_id"] = staple_id

	# 3. add extra graph nodes and edges for staple dangling ends
	for dangle in staple_dangles :

		i53_1, i53_2 = dangle
		dtype, staple_id, section_id, dom_id_joinedto, runlen = staple_dangles[dangle]

		# add new graph nodes at beginning or end of dangle (these are free-floating, not connected to scaffold)
		if dtype == "beginning_dangle" :
			G.add_node(i53_1)
			G.nodes[i53_1]["type"] = "dangle"
			G.nodes[i53_1]["staple_id"] = staple_id
		else :
			G.add_node(i53_2)
			G.nodes[i53_2]["type"] = "dangle"
			G.nodes[i53_2]["staple_id"] = staple_id

		# add the dangling edge
		G.add_edge(i53_1, i53_2)
		G.edges[i53_1, i53_2]["type"] = "dangle"
		G.edges[i53_1, i53_2]["scaffold_domain_id"] = -1 		# not a part of the scaffold
		G.edges[i53_1, i53_2]["runlen"] = runlen
		G.edges[i53_1, i53_2]["schematiclen"] = runlen
		G.edges[i53_1, i53_2]["source"] = i53_1 				# in direction of staple
		G.edges[i53_1, i53_2]["target"] = i53_2
		G.edges[i53_1, i53_2]["staple_id"] = staple_id
		G.edges[i53_1, i53_2]["section_id"] = section_id
		G.edges[i53_1, i53_2]["staple_sequence53"] = get_staple_section_sequence53(contactmap, staple_id, section_id)

		# update the scaffold domain on the origami, that the dangling edge joins to 
		# (note: this is particularly necessary if a staple has 2 sections, 1 joined to a scaffold domain, and 1 a dangle)
		# (in this case, these is no crossover. So, adding the dangle must fully set the properties of the scaffold domain.)
		edge_joined = [(u,v) for u,v,a in G.edges(data=True) if a["scaffold_domain_id"] == dom_id_joinedto]
		u, v = edge_joined[0]
		G.edges[u, v]["type"] = "ds_domain"
		G.edges[u, v]["staple_id"] = staple_id
		G.nodes[u]["staple_id"] = staple_id
		G.nodes[v]["staple_id"] = staple_id

		if dtype == "beginning_dangle" :
			G.edges[u, v]["section_id"] = section_id+1 	# the scaffold domain joined is after the dangle
		else :
			G.edges[u, v]["section_id"] = section_id-1	# the scaffold domain joined is before the dangle

		G.edges[u, v]["staple_sequence53"] = get_staple_section_sequence53(contactmap, staple_id, G.edges[u, v]["section_id"])

	return G











#
#
# SCHEMATIC LEN CORRECT
#
#

def schematiclen_corrected(G) :
	"""Sets all edges with schematic len of 0, to a small value
	schematiclen = 0 causes problems with the Kamad-Kawai springer embedder (division by 0 error)

	An edge length of 0 can arise e.g. for a 1nt staple dangle (??)
	"""

	corrected = False

	for u,v,a in G.edges(data=True) :
		if a["schematiclen"] < 1 :
			a["schematiclen"] = 0.01
			corrected = True

	return corrected




