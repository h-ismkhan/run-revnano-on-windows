"""Converts a scadnano .sc output file to a contact map CSV file

Note: It is assumed the origami design in scadnano has a LINEAR scaffold.

Note: For this converter to succeed: 
	- all insertions and deletions must not be on the end 5' or 3' bases of staples
	- where insertions and deletions are at bases where a staple is hybridised with the scaffold, 
		they should comes in PAIRS (one on the scaffold, one of matching length on the staple)
		(insertions and deletions on ssDNA parts of the scaffold are allowed to be singular)

Ben Shirt-Ediss, March 2021 / June 2022
"""

import sys
sys.path.append("../")

import os
from itertools import groupby

from contactmap import iohandler
from contactmap import contactutils
from contactmap import oxdna2contact



#
#
# GET CORE SCAFFOLD INFORMATION
#
#


def insertions_to_dict(insertions_list) :
	"""Converts an insertions list to a friendlier dictionary
	"""

	insertions_dict = {}
	for i in insertions_list :
		insertions_dict[ i[0] ] = i[1]

	return insertions_dict




def get_i53map(scaffold_domains) :
	"""Returns dictionary: i53 on scaffold --> (helix, offset) on the design grid
	"""

	i53map = {}
	i53 = 0		# start at scaffold 5'

	deletions_set = set()		# set of (helix, offset) tuples, which are deletions on the scaffold
	insertions_set = set()		# set of (helix, offset) tuples, which are insertion points on the scaffold

	# go through scaffold helixes, from 5' of scaffold to 3' of scaffold
	for scaffold_domain in scaffold_domains :
		
		# is scaffold domain a LOOPOUT?
		if "loopout" in scaffold_domain :
			# these loopout scaffold bases are not connected to anything
			N = scaffold_domain["loopout"]   # number of bases to loop out
			for _ in range(0, N) :
				i53map[i53] = -1
				i53 += 1
			continue

		# scaffold domain is a HELIX (which may or may not be paired with staples)

		insertions = {}
		if "insertions" in scaffold_domain : insertions = insertions_to_dict(scaffold_domain["insertions"])

		deletions = []
		if "deletions" in scaffold_domain : deletions = scaffold_domain["deletions"]		

		# go the correct way through the offsets, depending on which way the
		# scaffold domain helix is pointing in the scadnano design
		offset_range = list( range(scaffold_domain["start"], scaffold_domain["end"]) )	# start to end-1	
		if scaffold_domain["forward"] == False :		
			offset_range = offset_range[::-1]											# reverse range: from end - 1 to start		

		for offset in offset_range :

			helix_and_offset_tuple = (scaffold_domain["helix"], offset)

			# deletions are skipped (these offsets do not physically exist)
			if offset in deletions :
				deletions_set.add(helix_and_offset_tuple)
				continue

			i53map[i53] = helix_and_offset_tuple
			i53 += 1

			# if offset was an insertion point, now add the extra insertion bases after it
			if offset in insertions :
				insertions_set.add(helix_and_offset_tuple)
				N = insertions[offset]	# number of extra bases to insert
				for _ in range(0, N) :
					i53map[i53] = -1	# this scaffold base does not map to a helix, offset in the design
					i53 += 1					

	return i53map, deletions_set, insertions_set




def print_i53map(i53map) :

	for i53 in sorted(i53map.keys()) :
		print(i53, i53map[i53])




def invert_i53map(i53map) :
	"""Returns dictionary: (helix, offset) on the design grid --> i53 on scaffold
	"""

	i53map_inverse = {}

	for i53 in i53map :

		helix_and_offset_tuple = i53map[i53]
		if helix_and_offset_tuple != -1 :
			i53map_inverse[helix_and_offset_tuple] = i53

	return i53map_inverse




def get_scaffold_info(scdata) :
	"""Returns scaffold sequence, the inverse i53 map and all deletions 
		and insertion points in the scaffold
	"""

	for strand in scdata["strands"] :
		if "is_scaffold" in strand :
			
			# make the i53 map inverse
			i53map, deletions_set, insertions_set = get_i53map(strand["domains"])

			i53map_inverse = invert_i53map(i53map)

			# if scaffold sequence is not specified, make it NNNNNN....				
			scaffold53 = ""
			if not ("sequence" in strand) :	
				scaffold53 = "N" * (len(i53map))
			else :
				scaffold53 = strand["sequence"]

			return scaffold53, i53map_inverse, deletions_set, insertions_set

	return None, None, None, None 	# no scaffold strand


























#
#
# 1. STAPLES TO SCAFFOLD MAPPING
# 		and 2. SCAFFOLD TO STAPLES MAPPING
#
# Note: 2 can be inferred from 1, if scaffold length and sequence is known
# 1 can be inferred from 2, if details of all ssDNA staple sections are known
# Both mappings are necessary for a complete contact map description


def record_scaffold_to_staple_binding(sc, i53, scaffold53, staple_id, staple_section_id, staple_base_id) :

	if i53 != -1 : 
		sc[i53] = 	{	"base_letter"			: scaffold53[i53],
						"scaffold_domain_id"	: None,
						"staple_id"				: staple_id,
						"staple_section_id"		: staple_section_id,
						"staple_base_id"		: staple_base_id
					}




def get_mappings(scdata) :
	"""Map all staples to the scaffold (to make st mapping). 
	When a staple base is mapped to the scaffold, the reverse also happens: 
	the scaffold base is mapped to the staple (to make the sc mapping).
	"""

	scaffold53, i53map_inverse, deletions_set, insertions_set = get_scaffold_info(scdata)

	if scaffold53 == None :		# no scaffold = not an origami
		return None, None

	st = {}	  # id of staple : dict of how staple binds to scaffold
	sc = {}	  # i53 of scaffold base : dict of which staple base this scaffold base binds to

	# make empty sc mapping
	for i53 in range(0, len(scaffold53)) :
		record_scaffold_to_staple_binding(sc, i53, scaffold53, -1, -1, -1)

	# go through all staples
	staple_id = 0
	for strand in scdata["strands"] :

		# ignore scaffold strand
		if "is_scaffold" in strand :
			continue	

		# attribute to compute for this staple
		base_i53 = []

		# go through domains of the staple

		# NOTE: "domains" of a strand in scadnano are simply runs of bases on the same row ("helix") of the design grid
		# only when the row changes, or when there is a loopout, does the domain change.
		# However, scadnano strand domains not always correspond to staple "sections" in the sense of a domain-level graph.
		# For instance, a run of staple bases initially off the scaffold, and running onto the scaffold is a single strand "domain"
		# but it is two separate staple sections: a dangle staple section and a hybridised staple section. 
		# Therefore, I first compute base_i53, and then I compute sections of the staple from this.

		for staple_domain in strand["domains"] :

			# is this staple section a special "loopout"? i.e. an interior ssDNA part of a staple
			if "loopout" in staple_domain :
				# yes: mark all bases in the loopout as not hybridising
				N = staple_domain["loopout"]	# number of bases to loop out
				base_i53 += [-1] * N
				continue

			# trace staple domain left to right, or right to left depending on its orientation
			offset_range = list( range(staple_domain["start"], staple_domain["end"]) )	# staple 5' to 3' goes left to right in scadnano main view
			if staple_domain["forward"] == False :											
				offset_range = offset_range[::-1]											# staple 5' to 3' goes right to left in scadnano main view

			# loop through all offsets, going from 5' to 3' in the staple domain
			for offset in offset_range :

				helix_and_offset_tuple = (staple_domain["helix"], offset)			

				# --------- deletion case ---------
				# if staple base is a deletion, skip it
				if helix_and_offset_tuple in deletions_set :
					continue
				# --------- deletion case ---------	

				# record i53 where staple base binds scaffold
				# (or -1, if base does not bind scaffold)
				i53 = i53map_inverse.get(helix_and_offset_tuple, -1)

				# --------- insertion case ---------
				# if base is an insertion, add in the extra bases BEFORE the insertion is reached
				if helix_and_offset_tuple in insertions_set :

					try :
						insertions = insertions_to_dict(staple_domain["insertions"])
					except KeyError:
						print("\nError: There is an base insertion point on the scaffold at helix number %d, offset %d. However, the staple at this point does not have a matching base insertion point. Therefore the scadnano design has an error.\n" % (helix_and_offset_tuple))
						quit()

					N = insertions[offset] # number of bases to insert

					for n in range(i53 + N, i53, -1) :
						base_i53.append(n)
				# --------- insertion case ---------

				base_i53.append(i53)
				

		# ignore free floating staples, not connected to scaffold in any way
		if sum(base_i53) == -len(base_i53) :
			continue

		# 1. update st information
		sequence53 = "N" * len(base_i53)		# a scadnano design may or may not have staple sequences. If not, write them as NNNNN....
		if "sequence" in strand : sequence53 = strand["sequence"]

		num_sections, section_lengths = oxdna2contact.get_staple_section_info(base_i53, False, len(scaffold53))

		st[staple_id] = { 	"sequence53" 					:	sequence53,
							"num_sections"					:	num_sections,
							"section_lengths"				:	section_lengths,
							"section_scaffold_domain_ids"	:	None,				# not set here, because they are only defined once all staples are placed.
							"base_i53"						:   tuple(base_i53)
						}
		
		# 2. update sc information
		for staple_base_id, i53 in enumerate(base_i53) :
			staple_section_id = oxdna2contact.base_in_which_staple_section(section_lengths, staple_base_id)
			record_scaffold_to_staple_binding(sc, i53, scaffold53, staple_id, staple_section_id, staple_base_id)				

		staple_id += 1

	return st, sc




def update_mappings(st, sc) :
	"""Associates a scaffold domain number to all scaffold bases, and all staple sections

	Assumes scaffold is LINEAR

	The passed dictionaries are updated by reference - there is no returned value.
	Domain ids are needed for a domain level graph, but are not strictly necessary for a bare contact map
	"""

	# 1. assign scaffold domain numbers to scaffold bases, now all staples are placed
	# Again, assumes a LINEAR scaffold. No correction for circular scaffold is applied.
	domain_id = 0
	ssid = "%d_%d" % (sc[0]["staple_id"], sc[0]["staple_section_id"])

	for i53 in range(0, len(sc)) :

		ssid_new = "%d_%d" % (sc[i53]["staple_id"], sc[i53]["staple_section_id"])

		if ssid_new != ssid :
			# staple section hybridised has changed
			domain_id += 1
			ssid = ssid_new

		sc[i53]["scaffold_domain_id"] = domain_id

	# 2. assign scaffold domain numbers to staple sections
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































#
#
# EXECUTE
#
#


if __name__ == "__main__" :
	
	print("\n--------------------------------")
	print("Scadnano --> Origami Contact Map")
	print("--------------------------------\n")

	if len(sys.argv) != 2 :
		print("Usage: python3 scadnano2contact.py <origami_name>\n")
		quit()	

	origami_name = sys.argv[1]
		
	sc_filename = "%s.sc" % (origami_name)
	if not os.path.exists("_assets/%s" % (sc_filename)) :
		print("Error: scadnano file '%s' does not exist in the _assets/ directory.\n" % (sc_filename))
		quit()

	print("scadnano design file:   %s\n" % (sc_filename))

	scdata = iohandler.read_scadnano_json("_assets/%s" % (sc_filename))

	if isinstance(scdata, str) :
		print("Error: %s\n" % scdata)
		quit()

	print("Converting to origami contact map....")	

	st, sc = get_mappings(scdata)

	if st == None :
		print("Error. No scaffold strand has been assigned in the scadnano design.\n")
		quit()

	update_mappings(st, sc)

	contactmap = [st, sc, False]	# always have a linear scaffold, from scadnano

	print("[Done]")		

	print("Verifying origami contact map.....")
	success = contactutils.verify_contact_map(contactmap)
	if success :
		print("[Success]")
	else :
		print("[Failed]\n")		# write contact map anyway

	cm_filename = "%s.csv" % (origami_name)
	print("Writing contact map to _assets/%s....." % (cm_filename))
	iohandler.write_origami_contactmap(contactmap, "_assets", cm_filename, madeby="scadnano2contact.py")
	print("[Done]\n")			

	