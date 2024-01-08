"""REVNANO functions for manipulating dna/rna sequences, 
querying origami scaffold, and miscellaneous tasks

Ben Shirt-Ediss, 2020, 2023
"""



#
#
# RAW DNA/RNA SEQUENCE OPERATIONS
#
#


def star_delim_valid(seq) :

	if "*" not in seq :
		return True

	if seq.count("*") % 2 == 1 :
		return False 	# odd number of * delim in string not allowed

	if "**" in seq :
		return False 	# empty base sequence between delimiters not allowed

	return True



def sequence_type(seq) :
	"""Given a string of uppercase letters, returns if its a nucleic acid sequence, or text

	The string is allowed to contain an even number of * delimiters, where at least
	1 base must be present between each star pair

	Returns:
	"NUCLEIC ACID" if its a valid nucleic acid sequence (which may contain * delimiters or not)
	"UT MISTAKE" if its a nucleic acid sequence containing both RNA and DNA bases by mistake 
	(Some origami sequence listings have this error. Sequences are either DNA or RNA.)
	"* DELIMITER MISTAKE" if and odd number of * delimiters are present, or if a 
	delimiter pair does not contain at least one base letter
	"TEXT" if it contains letters other than the 5 nucleic acid DNA/RNA bases
	"""
	
	if not star_delim_valid(seq) :
		return "* DELIMITER MISTAKE"

	dna = [s in ["A", "C", "G", "T", "*"] for s in seq]
	is_dna = (sum(dna) == len(seq))

	rna = [s in ["A", "C", "G", "U", "*"] for s in seq]
	is_rna = (sum(rna) == len(seq))

	mistake = [s in ["A", "C", "G", "U", "T", "*"] for s in seq]
	is_mistake = (sum(mistake) == len(seq))

	if is_dna or is_rna :
		return "NUCLEIC ACID"
	elif is_mistake : 		# important to check AFTER is_dna or is_rna
		return "UT MISTAKE"
	else :
		return "TEXT"



def dangles_loopouts(seq53) :
	""" Given a staple sequence with * delimiters, returns

	DDDDDSSSSLLLLSSSSSSLLLLSSSSSDDDDD

	- 5 dangle
	- 3 dangle
	- midlist [("substaple", seq), ("loopout", seq), ("substaple", seq)]

	Note: seq53 is assumed to contain none, or an even number of delimiters where each pair encloses a non-empty base sequence
	"""

	dangle5 = None
	dangle3 = None
	midlist = []

	chunks = []
	for c in seq53.split("*") :
		if c.strip() != "" :
			chunks.append(c.strip())

	if seq53.startswith("*") :
		dangle5 = chunks.pop(0)

	if seq53.endswith("*") :
		dangle3 = chunks.pop()

	# remaining chunks alternate between substaple and loopout
	
	for seq in chunks :

		seqtype = "substaple"
		if len(midlist) % 2 == 1 :
			seqtype = "loopout"

		midlist.append( (seqtype, seq) )

	return dangle5, dangle3, midlist





def bases_complementary(b1, b2) :
	"""Returns 1 if two bases are complementary For DNA-DNA, RNA-RNA or RNA-DNA/DNA-RNA hybrid Watson-Crick pairings

	Returns 0 if bases not complementary
	GU wobble pairs NOT COUNTED

	RNA - 	RNA
	A	-	U *
	C	-	G *
	G	-	C *
	U	-	A *

	DNA	-	DNA  (canonical base pairing)
	A	-	T *
	C	-	G *
	G	-	C *
	T	-	A *

	RNA	-	DNA
	A	-	T *
	C	-	G *
	G	-	C *
	U	-	A *

	Overall canonical Watson-Crick base pairs:

	b1		b2
	G 		C
	C 		G
	A 		T/U
	T 		A
	U 		A

	"""

	if b1=="G" and b2=="C" : return 1
	if b1=="C" and b2=="G" : return 1
	if b1=="A" and b2=="T" : return 1
	if b1=="A" and b2=="U" : return 1
	if b1=="T" and b2=="A" : return 1
	if b1=="U" and b2=="A" : return 1

	return 0	





























#
#
# SCAFFOLD OPERATIONS
#
#


def linear_circular(scaffold_type, scaffold_len, i53) :
	"""Allows a scaffold sequence -- specified as a linear string in the 5' to 3' direction -- to be used in a linear or a circular way

	Given i53, returns the effective i53 index to use to access the linear string of scaffold

	When CALCULATING the base index of any base that is not the section 5 prime,
	this function MUST be used to make circular scaffolds safe

	
	"""

	if scaffold_type == "LINEAR" :

		if (i53 < 0) or (i53 > (scaffold_len - 1)) :
			return -1	# index goes off end of linear scaffold!
		else :
			return i53

	elif scaffold_type == "CIRCULAR" :

		if i53 >= 0 :   # positive index, which can go off right end of scaffold53
			return i53 % scaffold_len
		else :			# negative index, which can go off left end of scaffold53
			idx_circular = scaffold_len - (abs(i53) % scaffold_len)
			if (abs(i53) % scaffold_len) == 0 :
				idx_circular = 0			
			return idx_circular
	else :
		print("linear_circular not given valid scaffold type")
		quit()





def get_run(origami, subseq53) :
	"""From EVERY base position on the scaffold, returns how many consecutive bases staple subsequence subseq53 
	is complementary for (the "run" of consecutive complementary bases)

	Returned dictionary (run) has form:
		i35: number of consecutive bases complementary to subseq53 (going from i35 towards 5' of scaffold)

	If the origami has a circular scaffold, this is accounted for (the beginning of the scaffold sequence is joined to the end)
	"""

	run = {}

	if subseq53 == "" : 
		return run

	# for each base index along scaffold, find what size the complementary run with subseq53 is

	scaffold35 = origami["scaffold53"][::-1]

	# i moves the staple window along the scaffold (from scaf 3' to 5')
	# j moves along bases within the window (from staple 5' to 3')

	for i in range(0, len(scaffold35)) :		# i is an i35 index

		run[i] = 0
		b_sc = scaffold35[i]

		# at the current scaffold base,
		# start at 5 prime of staple, go to 3' of staple, and see how many bases are complementary in a row
		# test all bases in staple sequence (in staple 5' to 3' direction), at the current scaffold base
		for j in range(0, len(subseq53)) :		# j is an addition to the i35 index

			b_st = subseq53[j]

			if bases_complementary(b_sc, b_st) :
				run[i] += 1
				
				# get index of next scaffold base in window
				# taking into account scaffold may be linear or circular
				i35 = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i + j + 1)
					# note: i+j+1 is an i35 index
					# but linear_circular can take either an i53 index, or i35 index

				if i35 == -1 :
					break  # end of linear scaffold reached -  next base does not exist, matching run is over at scaffold base i

				b_sc = scaffold35[i35]
			else :
				break # bases not complementary, matching run is over at i

	
	# i35 is base index from scaffold 3 prime
	# i53 is base index from scaffold 5 prime
	# to convert between them:
	# i53 = len(scaffold53) - 1 - i35	
	
	# run is indexed by the i35 index

	return run





def scaffold_pairs(origami) :
	"""All pairs of bases joined by the scaffold backbone
	"""
	
	scaffold_pair_list = []
	for i53, base in enumerate(origami["scaffold53"]) :

		if i53 < len(origami["scaffold53"]) - 1:
			scaffold_pair_list.append( (i53, i53+1) )
		else :
			if origami["scaffold_type"] == "CIRCULAR" :
				scaffold_pair_list.append( (i53, 0) )

	return scaffold_pair_list





def staple_crossover_pairs(origami, staple_id) :
	"""Returns the scaffold base index pairs where the staple specified makes crossovers

	Staple is traversed from 5 prime to 3 prime of staple

	A list of tuples is returned, a tuple for each crossover
	[(3 prime of previous staple section, 5 prime of next staple section), ...]
	"""

	St = origami["staples_placed"]

	crossover_pair_list = []

	# get staple sections
	all_sections = St[staple_id]		
	num_sections = len(all_sections)

	if all_sections == [] :  # staple not placed yet
		return []

	if num_sections == 1 :	# no crossover exists on the staple
		return []

	for i, section in enumerate(all_sections) :

		if i == (num_sections - 1) :
			continue	# don't count the last section - it has no crossover to a next section

		# add crossover pair
		
		# calc i53 of 3 prime
		i53_section_5prime, runlen = section
		i53_section_3prime = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - runlen + 1))
		
		# get i53 of 5 prime of following section
		i53_next_section_5prime, _ = all_sections[i+1]

		crossover_pair_list.append( (i53_section_3prime, i53_next_section_5prime) )
		
	return crossover_pair_list





def staple_section_pairs(origami, staple_id, section_id) :
	"""Returns pairs of bases joined on the scaffold by this staple section
		Also returns length of section, i53 of 5 prime, and i53 of 3 prime
	"""
	
	section_pair_list = []

	section = origami["staples_placed"][staple_id][section_id]
	
	i53_section_5prime, runlen, _ = section	
	i53_section_3prime = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - runlen + 1))

	current_pair = []
	r = runlen
	while r > 0 :
		i53 = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - r + 1))

		if len(current_pair) == 0 :
			current_pair.append(i53)
		else :
			current_pair.append(i53)
			section_pair_list.append( tuple(current_pair) )
			current_pair = [i53]

		r -= 1

	return section_pair_list, runlen, i53_section_5prime, i53_section_3prime





def cast_staple_footprint_on_scaffold(origami, staple_id, footprint, num_staple_routes) :
	"""Given a staple footprint, claims bases on the origami scaffold where ALL routes of the staple cross that base
	""" 
	for i53, count in enumerate(footprint) :
	
		if count == num_staple_routes :		# all routes of staple go through this scaffold base
			origami["footprints"][i53].add(staple_id)





def scaffold_segment_get_claimed_percent(origami, staple_id, i53_start, runlen) :
	"""On the strip (i53_start, runlen) on the origami scaffold, returns how many other staples (other than staple_id) have claimed these bases
	"""

	bases_claimed_by_other_staples = 0

	for i53 in range(i53_start, i53_start-runlen, -1) :	

		# note: i53 is staple frag 5'
		# it is decremented, to get to staple frag 3'

		# correct index in case of circular scaffold
		i53 = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53)

		occ = origami["footprints"][i53]

		if len(occ) > 0 :
			# some staples already occupying base i53

			if staple_id in occ :
				if len(occ) > 1 :
					# this staple, but also other staple(s) occupy i53
					bases_claimed_by_other_staples += 1
			else :
				# this staple is not at i53
				bases_claimed_by_other_staples += 1

	return (bases_claimed_by_other_staples / runlen)




def scaffold_segment_fully_nests_a_staple_section(origami, staple_id, i53_start, runlen) :
	"""Returns True if claiming scaffold segment (i53, runlen) would result in fully nesting AT LEAST one other staple section placed already

	Defn: Another staple section is "fully nested" in the scaffold segment i53, runlen
	when (1) there is a staple id (not equal to the one passed) in the segment
	that (2) does not go outside the boundaries of segment

	Note: a staple section is not expected to nest other smaller sections of the same staple,
	because srt.get_match_candidates() will truncate the larger section immediately
	when another section of the same staple is placed

	Note: it is possible to become nested by 2 staples, in the pathological limit. 
	However, this case will cause a fail at stage 3, and so is not considered here.
	"""	

	# step 1: find the ids of staples already in this segment
	staple_ids_in_segment = set()

	for i53 in range(i53_start, i53_start-runlen, -1) :	
		# correct index in case of circular scaffold
		i53 = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53)
		staple_ids_in_segment = staple_ids_in_segment.union( origami["footprints"][i53] )
	
	# step 2: these staple ids must move outside the segment, in order NOT to be fully nested by it
	i53_before_seg_start = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53_start+1)
	i53_after_seg_end = linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53_start-runlen)

	for staple_id2 in staple_ids_in_segment :
		if staple_id2 != staple_id :
			
			if (i53_before_seg_start != -1) and (staple_id2 in origami["footprints"][i53_before_seg_start]) :
				continue  # ok: staple_id2 goes outside the start boundary and so is not nested

			if (i53_after_seg_end != -1) and (staple_id2 in origami["footprints"][i53_after_seg_end]) :
				continue  # ok: staple_id2 goes outside the end boundary and so is not nested

			# other staple id does not go out of start or end, so is completely nested

			# note: when i53 at start or end is -1, then this is the end of a linear scaffold
			# and by definition in this case, the nested staple cannot pass over the physical scaffold end
			return True

	return False




def remove_staple_footprint_from_scaffold(origami, staple_id) :
	"""Note: does not work with detailed footprint
	"""
	for i53 in origami["footprints"] :

		claimants = origami["footprints"][i53]   # deep copy

		if staple_id in claimants :
			claimants.remove(staple_id)




def remaining_overlaps(origami) :
	"""The current number of staple bases overlapping on the scaffold 
	"""

	r = 0
	for i53 in origami["footprints"] :
		ov = len(origami["footprints"][i53])
		if ov > 0 : 	# base position occupied by staple or staples
			r += (ov-1)

	# count number of scaffold bases hybridised with staples

	return r



























#
#
# INTERMEDIATE STATS
#
#

def count_scaffold_bases_claimed(origami) :
	"""Returns how many scaffold bases are currently claimed by staples
	in the footprints dictionary
	"""
	count = 0

	for i53 in origami["footprints"] :
		if len(origami["footprints"][i53]) > 0 :
			count += 1

	return count




def count_staples_placed(origami) :
	"""The number of staples that have been routed
	"""

	count = 0

	for staple_id in origami["staples_placed"] :
		if len(origami["staples_placed"][staple_id]) > 0 :
			count += 1

	return count






















#
#
# FINAL STATS
#
#


def final_stats(contactmap) :

	st, sc, _ = contactmap

	# total scaffold bases hybridised
	sbh = 0
	for i53 in sorted(sc.keys()) :
		if sc[i53]["staple_id"] != -1 :
			sbh += 1

	return sbh, len(st)				# total scaffold bases hybridised, total physical staples placed





















