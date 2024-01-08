"""Builds origami contact map from the origami oxDNA topology (.top) and config (.oxdna) files

Note 1: the config file should be the one exported by a design tool, without relaxation. 
It should detail initial nucleotide positions only, and not have trajectory information embedded.

Note 2: properly handles both LINEAR and CIRCULAR scaffolds.

Ben Shirt-Ediss, September 2021, November 2021
"""

import sys
sys.path.append("../")

import math
import numpy
from itertools import groupby
import bisect
import os.path
import warnings

from contactmap import iohandler
from contactmap import contactutils



# ----- Consensus parameters for base pair detection ----------

# Derived from oxdna2param: see docs for oxdna2param

EUCLIDEAN_CUTOFF = 1.75  	# oxDNA simulation length units
ANGLE_MAX = 8.75 			# degrees

# -------------------------------------------------------------





#
#
# oxDNA TOPOLOGY FILE OPERATIONS
#
#

def get_strand_5prime_3prime(strand_id_target, top) :
	"""Given a strand id, returns the nucleotide id of its 5' nucleotide, 3' nucleotide
	and whether the strand is circular or not

	For robustivity, minimal assumptions are made:
	Strands in the topology file can be broken into multiple non-consecutive segments
	and be connected from 3 to 5 prime as the file advances, or the other way around (both ways are common)

	Note: well-formedness of the topology file is assumed
	All linear strands are assumed to have 3' and 5' end nucleotides defined
	"""

	strand_exists = False
	id_5prime = None
	id_3prime = None
	strand_circular = False

	nuc = [] 		# list of all nucleotide ids on the target strand

	for nucleotide_id in sorted(top.keys()) :
		strand_id, base, id3prime, id5prime = top[nucleotide_id]

		if strand_id == strand_id_target :

			strand_exists = True

			if id3prime == -1 :
				id_3prime = nucleotide_id
			if id5prime == -1 :
				id_5prime = nucleotide_id 		# note: strand being a single nucleotide is allowed

			nuc.append(nucleotide_id)

	# if the strand exists, but has no definitive 5' or 3' end, then it must be a circular strand
	if strand_exists and (id_3prime == None) and (id_5prime == None):

		strand_circular = True

		# find first occurring, and last occuring nucleotides of this strand in the topology file
		# these muat be joined if the strand is circular, even if the strand is defined in 
		# multiple non-consecutive segments
		first_nucleotide_id = min(nuc)
		last_nucleotide_id = max(nuc)

		_, _, _, id5 = top[first_nucleotide_id]
		if id5 == last_nucleotide_id :
			# 5' points to the last nucleotide in the loop:
			# decreasing line numbers go toward 5' and increasing line numbers go toward 3'. 
			# Hence scaffold is defined 5' to 3' as line number advances
			id_5prime = first_nucleotide_id
			id_3prime = last_nucleotide_id
		else :
			# 3' must instead point to the last nucleotide in the loop:
			# decreasing line numbers go toward 3' and increasing line numbers go toward 5'. 
			# Hence scaffold is defined 3' to 5' as line number advances
			id_5prime = last_nucleotide_id
			id_3prime = first_nucleotide_id

	return id_5prime, id_3prime, strand_circular




def walk_strand53(id_5prime, id_3prime, top) :
	"""Walks a strand from 5' to 3' and returns the nucleotide ids, and sequence of strand (as lists)
	"""

	n = []	# nucleotide ids
	s = []	# sequence

	nucleotide_id = id_5prime

	# start at 5' nucleotide and follow 3' pointers, until 3' reached
	while True :
		_, base, id3, _ = top[nucleotide_id]
		n.append(nucleotide_id)
		s.append(base)

		if nucleotide_id == id_3prime :
			break

		nucleotide_id = id3

	return n, s




def oxdna_topology_file_to_useful_format(top) :
	"""Converts the dictionary of the oxDNA topology file into a more useful dictionary
	"""

	topology = {}

	# applies to all strands
	topology["strand_nucleotide_ids_53"] 	= {}	# strand id : [list of nucleotide ids, from 5' to 3']
	topology["strand_sequence_53"] 			= {}	# strand id : [list of base sequence, 5' to 3']
	topology["nucleotide_id_to_strand_id"] 	= {}	# nucleotide id : strand id
	topology["strand_circular"] 			= {}	# strand id : True or False
	
	# applies only to scaffold strand
	topology["scaffold_strand_id"]			= -1	# strand id of scaffold strand (the longest strand in the topology file)
	topology["scaffold_nucleotide_ids_53"]	= []	# list of scaffold nucleotide ids, going 5' to 3' on scaffold
	topology["scaffold_sequence_53"]		= []	# list of scaffold sequence
	topology["scaffold_i53"]				= []	# list of i53 indexes of above scaffold nucleotides

	L = {}		# strand id : length  (for finding scaffold strand)

	for nucleotide_id in sorted(top.keys()) :
		strand_id, _, _, _ = top[nucleotide_id]

		topology["nucleotide_id_to_strand_id"][nucleotide_id] = strand_id

		# has this strand been encountered before?
		if topology["strand_nucleotide_ids_53"].get(strand_id, -1) == -1 :
			# no, add it

			id_5prime, id_3prime, strand_circular = get_strand_5prime_3prime(strand_id, top)
			
			n, s = walk_strand53(id_5prime, id_3prime, top)

			topology["strand_nucleotide_ids_53"][strand_id] = n
			topology["strand_sequence_53"][strand_id] = s
			topology["strand_circular"][strand_id] = strand_circular
			L[strand_id] = len(s)

	# find scaffold strand id (the longest strand)
	# https://stackoverflow.com/questions/268272/getting-key-with-maximum-value-in-dictionary
	# Note: if there are >1 max length strands, the first is returned
	topology["scaffold_strand_id"] = max(L, key=L.get)
	topology["scaffold_nucleotide_ids_53"] = topology["strand_nucleotide_ids_53"][ topology["scaffold_strand_id"] ]
	topology["scaffold_i53"] = list( range(0, len(topology["scaffold_nucleotide_ids_53"])) )
	topology["scaffold_sequence_53"] = topology["strand_sequence_53"][ topology["scaffold_strand_id"] ]

	return topology































#
#
# oxDNA CONFIG FILE OPERATIONS
#
#

def xyz_staple_nucleotides(topology, config) :
	"""Make x, y and z dictionaries of the positions of staple nucleotides on the oxDNA config file
	
	Each dictionary is of format:
		x/y/z float coordinate: [list of staple nucleotide ids with that x/y/z float coordinate]

	The keys of these dictionaries can be sorted to quickly find which staple nucleotides occupy
	different box regions of x,y,z space
	"""

	xdict = {}; ydict = {}; zdict = {}

	for nucleotide_id in config :
		if nucleotide_id not in topology["scaffold_nucleotide_ids_53"] :
			
			# its a staple nucleotide, get its position
			x,y,z = config[nucleotide_id][0:3]

			# round float coordinates, to get a controllable dictionary key
			# https://stackoverflow.com/questions/23721230/float-values-as-dictionary-key
			x = round(x, 10)
			y = round(y, 10)
			z = round(z, 10)

			if xdict.get(x, -1) == -1 :
				xdict[x] = [nucleotide_id]
			else :
				xdict[x].append(nucleotide_id)

			if ydict.get(y, -1) == -1 :
				ydict[y] = [nucleotide_id]
			else :
				ydict[y].append(nucleotide_id)

			if zdict.get(z, -1) == -1 :
				zdict[z] = [nucleotide_id]
			else :
				zdict[z].append(nucleotide_id)

	return xdict, ydict, zdict




def data_in_range(mylist, minval, maxval) :
	"""Returns all numbers in mylist that are within closed interval [minval,maxval]
	(Note that a closed interval means that the bounds are themselves included)

	mylist is a list of DISTINCT floating point values, ordered ascending
	"""

	if minval > maxval :
		return []

	minval_idx = bisect.bisect(mylist, minval) - 1
	maxval_idx = bisect.bisect(mylist, maxval) - 1

	"""
	simple bisect example:

	#       0 1 2   3 4 5 
	keys = [1,2,3.3,6,7,8]
	index = bisect.bisect(keys, r) - 1

	- keys is an sorted ascending list of floats
	- a pair of list indexes contain the floating point number r
	- returns the lower bound index of this pair (which may be r itself)
	- IF r >= last item, the lower bound index is the last index
	- IF r < first item, -1 is returned
	"""	

	# if lower bound indexes are the same then both bounds lie in the same interval between two data points of the list
	# and don't encapsulate any data points of the list
	# (this also catches case that both bounds are higher than the maximum data point, or both are too smaller than the minimum data point)
	if minval_idx == maxval_idx :
		return []

	if minval_idx == -1 or mylist[minval_idx] != minval :
		# range is minval_idx+1 to maxval_idx
		return mylist[minval_idx+1:maxval_idx+1]
	else :
		# range is minval_idx to maxval_idx  (minval is actually the lower bound)
		return mylist[minval_idx:maxval_idx+1]




def build_scaffold_nucleotide_neighbour_index(topology, config, detection_sphere_radius) :

	N = {}		# scaffold nucleotide id: [staple nucleotide id]
				# 							^ list of potential staple nucleotides that can base-pair the scaffold nucleotide

	xdict, ydict, zdict = xyz_staple_nucleotides(topology, config)
	
	xlist = sorted(xdict.keys()); ylist = sorted(ydict.keys()); zlist = sorted(zdict.keys())
	
	# go through all scaffold nucleotides
	for scaffold_nucleotide_id in topology["scaffold_nucleotide_ids_53"] :

		# details of current scaffold nucleotide
		r = config[scaffold_nucleotide_id][0:3]
		b = config[scaffold_nucleotide_id][3:6]
		x,y,z = r

		# get staple nucleotides in bounding cube, 
		# where current scaffold nucleotide lays at the centre of the cube
		in_x_range = set()
		for x_ in data_in_range(xlist, x-detection_sphere_radius, x+detection_sphere_radius) :	
			in_x_range.update( xdict[ round(x_, 10) ] )

		in_y_range = set()
		for y_ in data_in_range(ylist, y-detection_sphere_radius, y+detection_sphere_radius) :
			in_y_range.update( ydict[ round(y_, 10) ] )

		in_z_range = set()
		for z_ in data_in_range(zlist, z-detection_sphere_radius, z+detection_sphere_radius) :
			in_z_range.update( zdict[ round(z_, 10) ] )

		N[scaffold_nucleotide_id] = [staple_nucleotide_id for staple_nucleotide_id in set.intersection(in_x_range, in_y_range, in_z_range)]

	return N

























#
#
# GEOMETRIC BASE PAIR DETECTION
#
#


def vector_between_points_3d(u, v) :
	"""Vector pointing from u to v
	"""

	return ( (v[0]-u[0]), (v[1]-u[1]), (v[2]-u[2]) )




def euclidean_distance_3d(u, v) :

	# simply the magnitude of the vector connecting the two points
	return math.sqrt( (v[0]-u[0])**2 + (v[1]-u[1])**2 + (v[2]-u[2])**2 )




def vector_angle_degrees(p, q) :

	# cosine of the radian angle between vectors. Angle calc by vector dot products divided by multiplication of vector magnitudes
	# -1 <= cos_angle <= 1
	# the cos() of which angle is -1: 180 degrees (pi radians)
	# the cos() of which angle is 1:  0 degrees (0 radians)	
	cos_angle = numpy.dot(p, q) / (numpy.linalg.norm(p) * numpy.linalg.norm(q))

	# sometimes numerical imprecision can cause cos_angle to fall out of range -1 <= cos_angle <= 1
	# and arccos does not like that, so correct here
	if abs(cos_angle) > 1 :
		if cos_angle > 1 :
			return 0.0
		else :
			return 180.0

	# arc cos gives radian angle (only defined for -1 <= cos_angle <= 1)
	angle_rad = numpy.arccos(cos_angle)	

	# which is then converted to degrees
	return math.degrees(angle_rad)




def is_candidate_base_pair(sc_r, sc_b, st_r, st_b, euclidean_cutoff, angle_max) :

	if euclidean_distance_3d(sc_r, st_r) <= euclidean_cutoff :

		sc_st_vector = vector_between_points_3d(sc_r, st_r)
		st_sc_vector = vector_between_points_3d(st_r, sc_r)

		if (vector_angle_degrees(sc_b, sc_st_vector) <= angle_max) and \
				(vector_angle_degrees(st_b, st_sc_vector) <= angle_max) :
			return True

	return False




def fast_detect_base_pairings(topology, config) :
	"""Similar layout to stats_base_pairings() in oxdna2param
	However, for speed, now assumed is that all base pairs detected are monogamous, 
	and loops quit when FIRST base pair match is found
	"""

	N = build_scaffold_nucleotide_neighbour_index(topology, config, EUCLIDEAN_CUTOFF)

	bases_paired = {}
	bases_paired["scaffold_to_staple"] = {} 	# scaffold nucleotide id: staple nucleotide id paired to
	bases_paired["staple_to_scaffold"] = {} 	# staple nucleotide id: scaffold nucleotide id paired to

	# each time a scaffold base pairs uniquely to a staple base, record the staple base
	for scaffold_nucleotide_id in N :

		sc_r = config[scaffold_nucleotide_id][0:3]
		sc_b = config[scaffold_nucleotide_id][3:6]

		if sc_b == (0,0,0) :	# skip over scaffold nucleotides with no backbone vector 
			continue			# (i.e. scaffold nucleotides in scadnano scaffold loopouts)

		for staple_nucleotide_id in N[scaffold_nucleotide_id] :

			st_r = config[staple_nucleotide_id][0:3]
			st_b = config[staple_nucleotide_id][3:6]

			if st_b == (0,0,0) :	# skip over staple nucleotides with no backbone vector
				continue			# (i.e. staple nucleotides in scadnano staple loopouts)

			if is_candidate_base_pair(sc_r, sc_b, st_r, st_b, EUCLIDEAN_CUTOFF, ANGLE_MAX) :

				# raise error and QUIT if duplicates found
				if (scaffold_nucleotide_id in bases_paired["scaffold_to_staple"]) or \
					(staple_nucleotide_id in bases_paired["staple_to_scaffold"]) :
					return False

				# staple nucleotide pairing to scaffold found
				bases_paired["scaffold_to_staple"][scaffold_nucleotide_id] = staple_nucleotide_id
				bases_paired["staple_to_scaffold"][staple_nucleotide_id] = scaffold_nucleotide_id

				# move directly to next scaffold base for speed
				break

	return bases_paired


























#
#
# BUILD CONTACT MAP
#
#


def get_staple_section_info(base_i53, scaffold_circular, scaffold_len) :
	"""Returns number of staple sections, and a semicolon delimited string of section lengths,
	given the i53 index for each base on the staple
	"""

	# base_i53, are the i53 indexes bound to on the scaffold, as you traverse the staple in 5' to 3' direction
	i53_previous = None

	# find length of staple sections
	s = 1 			# length of current section, in bases
	seclen = []

	for i, i53 in enumerate(base_i53) :
		if i < len(base_i53) - 1: 					# go up to penultimate i53
			pair = (i53, base_i53[i+1])

			a = sum([v == -1 for v in pair]) == 1  	# any pair containing a single -1 is a border
			b = abs(pair[0] - pair[1]) > 1 			# any pair with an i53 distance greater than 1 is a border

													# ...UNLESS the scaffold is circular, 
													# and the i53 change corresponds to going around to other end of scaffold
			if scaffold_circular :
				if (pair == (0, scaffold_len-1)) or (pair == (scaffold_len-1, 0)) :
					b = False

			if a or b:
				# border between sections found
				seclen.append(s)
				s = 1
			else :
				s += 1
		else :
			# at last element in list: add length of final staple section
			seclen.append(s)

	# format section length list, into a semicolon delimited string
	return len(seclen), seclen




def base_in_which_staple_section(section_lengths, staple_base_id) :
	"""Given the id of a base on the staple (5' is index 0), returns the id of the staple
	section this base is in (5' is index 0)
	"""

	S = []
	for section_id, l in enumerate(section_lengths) :
		S.extend([section_id] * l)
		
	return S[staple_base_id]




def build_contactmap(topology, bases_paired) :

	# ------------------------------------
	# Part 1: staples to scaffold
	# ------------------------------------
	
	st = {}		# id of staple : dict of how staple binds to scaffold

	staple_id = 0  	

	oxdna_strand_id_2_staple_id = {}		# Note: staple_id is not synonymous with the strand_id in the oxDNA topology file, 
											# hence this mapping between the two

	scaffold_circular = topology["strand_circular"][topology["scaffold_strand_id"]]
	scaffold_len = len(topology["scaffold_sequence_53"])

	for strand_id in sorted(topology["strand_nucleotide_ids_53"].keys()) :

		if strand_id == topology["scaffold_strand_id"] :
			continue

		# work out attributes for this staple
		
		# 1. sequence53
		sequence53 = "".join(topology["strand_sequence_53"][strand_id])
		
		# 2. base_i53
		base_i53 = []
		# go through all bases on staple, and find which scaffold bases they bind to
		for staple_nucleotide_id in topology["strand_nucleotide_ids_53"][strand_id] :
			# corresponding scaffold nucleotide, paired to this staple nucleotide
			scaffold_nucleotide_id = bases_paired["staple_to_scaffold"].get(staple_nucleotide_id, None)

			if scaffold_nucleotide_id == None :
				# this staple base is not paired with a scaffold base
				base_i53.append(-1)
			else :
				idx = topology["scaffold_nucleotide_ids_53"].index(scaffold_nucleotide_id)
				base_i53.append( topology["scaffold_i53"][idx] )
		
		# add staple
		# (don't add staple if all of its nucleotides are unpaired with scaffold -- such a staple is not part of the origami)
		if sum(base_i53) != -len(base_i53) :

			# 3. section_lengths, num_sections
			num_sections, section_lengths = get_staple_section_info(base_i53, scaffold_circular, scaffold_len)

			st[staple_id] = { 	"sequence53" 					:	sequence53,
								"num_sections"					:	num_sections,
								"section_lengths"				:   section_lengths,
								"section_scaffold_domain_ids"	:	None,				# to set in update_mappings()
								"base_i53"						:   tuple(base_i53)
							}

			oxdna_strand_id_2_staple_id[strand_id] = staple_id

			staple_id += 1

	# ------------------------------------
	# Part 2: scaffold to staples
	# ------------------------------------

	sc = {}		# i53 of scaffold base : dict of which staple base this scaffold base binds to

	for idx, scaffold_nucleotide_id in enumerate(topology["scaffold_nucleotide_ids_53"]) :
		
		# get i53 index of scaffold nucleotide
		i53 = topology["scaffold_i53"][idx]

		# get base letter of scaffold nucleotide
		base_letter = topology["scaffold_sequence_53"][idx]

		# find the staple base this scaffold nucleotide is hybridised with
		staple_nucleotide_id = bases_paired["scaffold_to_staple"].get(scaffold_nucleotide_id, None)

		staple_id = -1
		staple_section_id = -1
		staple_base_id = -1

		if staple_nucleotide_id != None :
			# yes, this scaffold nucleotide IS hybridised to a staple nucleotide

			# get the strand id, and the staple id, of the staple that this staple nucleotide is on
			strand_id = topology["nucleotide_id_to_strand_id"][staple_nucleotide_id]
			staple_id = oxdna_strand_id_2_staple_id[strand_id]

			# find index on the staple, of the base
			staple_base_id = topology["strand_nucleotide_ids_53"][strand_id].index(staple_nucleotide_id)

			# find what staple section the staple base is in
			staple_section_id = base_in_which_staple_section(st[staple_id]["section_lengths"], staple_base_id)
			
		# add scaffold base
		sc[i53] = 	{	"base_letter"			: base_letter,
						"scaffold_domain_id"	: None,
						"staple_id"				: staple_id,
						"staple_section_id"		: staple_section_id,
						"staple_base_id"		: staple_base_id
					}


	# ----------------------------------------------------
	# Part 3: number scaffold domains in sc, and update st
	# ----------------------------------------------------

	# 3a. number scaffold domains in sc
	domain_id = 0
	ssid = "%d_%d" % (sc[0]["staple_id"], sc[0]["staple_section_id"])

	for i53 in range(0, len(sc)) :

		ssid_new = "%d_%d" % (sc[i53]["staple_id"], sc[i53]["staple_section_id"])

		if ssid_new != ssid :
			# staple section hybridised has changed
			domain_id += 1
			ssid = ssid_new

		sc[i53]["scaffold_domain_id"] = domain_id	

	# we also need to consider a subtle case with circular scaffolds:
	# All scaffold end cases:
	# SCENARIO 																		RESOLUTION
	# 1. linear scaffold, physical scaffold ends are NOT spanned by same staple 	-- scaffold ends have different domain ids
	# 2. linear scaffold, physical scaffold ends spanned by the same staple 		-- scaffold ends have different domain ids
	# 3. circular scaffold, virtual ends are NOT spanned by same staple 			-- scaffold ends have different domain ids
	# 4. circular scaffold, virtual ends spanned by same staple 					-- scaffold ends have SAME domain ids
			# note: if beginning and end of the scaffold are single stranded, this is staple id -1
			# thus case 4 above applies also in the single stranded case

	# we must correct for scenario 4 above if it occurs

	if scaffold_circular :
		# if the staple section hybridising the last scaffold base at 3'
		# is the SAME staple section hybridising the first scaffold base at 5'
		# then make the first and last domains on the scaffold the same domain

		ssid5 = "%d_%d" % (sc[0]["staple_id"], sc[0]["staple_section_id"])
		ssid3 = "%d_%d" % (sc[scaffold_len-1]["staple_id"], sc[scaffold_len-1]["staple_section_id"])
		if ssid5 == ssid3 :
			# change domain id at 3' end of scaffold to the id of the first scaffold domain
			scaffold_domain_id3 = sc[scaffold_len-1]["scaffold_domain_id"]		
			i53 = scaffold_len-1
			while sc[i53]["scaffold_domain_id"] == scaffold_domain_id3 :
				sc[i53]["scaffold_domain_id"] = 0
				i53 -= 1


	# 3b. update st
	for staple_id in st : 

		# go through all staple bases in order, and record what scaffold domain id each base binds to
		domain_id_list = []
		for base_i53 in st[staple_id]["base_i53"] :

			# get scaffold domain id, of this staple base
			domain_id = -1 				# -1 i53's lead to -1 domain ids
			if base_i53 != -1 :
				domain_id = sc[base_i53]["scaffold_domain_id"]

			domain_id_list.append(domain_id)

		unique_domain_ids = [list(set(group))[0] for _, group in groupby(domain_id_list)]

		st[staple_id]["section_scaffold_domain_ids"] = unique_domain_ids

		
	contactmap = [st, sc, scaffold_circular]
								
	return contactmap
























#
#
# EXECUTE
#
#

if __name__ == "__main__" :

	print("\n----------------------------------------------")
	print("oxDNA --> Origami Contact Map")
	print("----------------------------------------------\n")

	if len(sys.argv) != 2 :
		print("Usage: python3 oxdna2contact.py <origami_name>\n")
		quit()

	origami_name = sys.argv[1]
		
	top_filename = "%s.top" % (origami_name)
	config_filename = "%s.dat" % (origami_name)
	config_filename2 = "%s.oxdna" % (origami_name)

	if os.path.exists("_assets/%s" % (top_filename)) :
		if not os.path.exists("_assets/%s" % (config_filename)) :
			if os.path.exists("_assets/%s" % (config_filename2)) :
				config_filename = config_filename2
			else :
				print("Error: oxDNA config file '%s.dat' (or '%s.oxdna') does not exist in the _assets/ directory.\n" % (origami_name, origami_name))
				quit()
	else :
		print("Error: oxDNA topology file '%s' does not exist in the _assets/ directory.\n" % (top_filename))
		quit()

	print("oxDNA topology file:   %s" % (top_filename))
	print("oxDNA config file:     %s" % (config_filename))
	print("")

	top = iohandler.read_oxdna_topology_file("_assets/%s" % (top_filename))
	config = iohandler.read_oxdna_config_file("_assets/%s" % (config_filename))

	print("Converting oxDNA topology file.....")
	topology = oxdna_topology_file_to_useful_format(top)
	print("[Done]")

	print("Finding base pairs in oxDNA config file.....")
	print("------ Base-pair Detection Parameters --------")
	print("%2.2f\tEUCLIDEAN_CUTOFF (oxDNA distance units)" % (EUCLIDEAN_CUTOFF))
	print("%2.2f\tANGLE_MAX (degrees)" % (ANGLE_MAX))
	print("----------------------------------------------")
	bases_paired = fast_detect_base_pairings(topology, config)

	if bases_paired == False :
		print("Failed. Duplicate base pairs detected. Get suitable EUCLIDEAN_CUTOFF and ANGLE_MAX base pair detection parameters for these oxDNA files by running oxdna2param.\n")
		quit()

	print("%d\tScaffold-staple base pairs found" % (len(bases_paired["scaffold_to_staple"])))

	oxview_filename = "_assets/%s.js" % (origami_name) 		# javascript to check the bases paired in oxview
	bp_filename = "_assets/%s.bp" % (origami_name)
	iohandler.write_bases_paired(bases_paired, oxview_filename, bp_filename)
	print("[Done]")

	print("Building origami contact map.....")
	contactmap = build_contactmap(topology, bases_paired)
	print("[Done]")

	print("Verifying origami contact map.....")
	success = contactutils.verify_contact_map(contactmap)
	if success :
		print("[Success]")
	else :
		print("[Failed]") 	# save contact map anyway
	
	cm_filename = "%s.csv" % (origami_name)
	print("Writing contact map to _assets/%s....." % (cm_filename))
	iohandler.write_origami_contactmap(contactmap, "_assets", cm_filename, madeby = "oxdna2contact.py")
	print("[Done]")
	print("")








