"""REVNANO : Stage 3 : Resolve staple-staple overlaps

Ben Shirt-Ediss, 2020, 2023
"""

import sys
sys.path.append("../")

from REVNANO import dna







#
#
# OVERLAP METHODS
#
#

def overlap3plus_exists(origami) :
	"""Returns true if there exists a base on the scaffold where 3 or more staples overlap
	"""

	for i53 in origami["footprints"] :
		claimants = origami["footprints"][i53]
		if len(claimants) > 2 :
			return True
	
	return False




def make_detailed_footprints_dictionary(origami) :
	"""A more detailed footprints datastructure is constructed just for this stage, using the St data structure
	origami["footprints"] is replaced by the new dictionary
		i53 : set( (staple_id, section_id), .......)

	Stage 3, needs a more detailed footprints dictionary since
	an overlap can cause a staple to intersect itself again, and this needs to be detected
	"""

	origami["footprints"] = {}
	for i53 in range(0, origami["scaffold_nt"]) :
		origami["footprints"][i53] = set()

	for staple_id in origami["staples_placed"] :
		all_sections = origami["staples_placed"][staple_id]
		for section_id, section in enumerate(all_sections) :
			# add all bases in this staple section to origami
			i53_section_5prime, runlen = section	
			r = runlen
			while r > 0 :
				i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - r + 1))
				origami["footprints"][i53].add( (staple_id, section_id) )
				r -= 1




def whois_overlapped_at_section_3prime(origami, staple_id, section_id) :
	"""Finds which other staple (and section of that staple) the 3 prime of the given staple section overlaps

	Returns:
		- (None, None, None) if this section does not overlap another staple at its 3 prime
		- (overlap_id, overlap_section_id, num_bases_overlapped), if this section does overlap another staple at its 3 prime

	Note that another staple may be overlapping this section at its 5 prime: we dont care about this
	We just care about what staple THIS section is overlapping at its 3 prime
	"""

	num_bases_overlapped = 0

	section = origami["staples_placed"][staple_id][section_id]
	
	i53_section_5prime, runlen = section
	i53_section_3prime = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - runlen + 1))

	# if no overlap at 3 prime, return
	claimants_3prime = set(origami["footprints"][i53_section_3prime])  # shallow copy
	if len(claimants_3prime) == 1 :  # i.e. only this section
		return (None, None, None)

	# if 3 or more overlaps at 3 prime, report error (pathological)
	if len(claimants_3prime) > 2 : 
		errormsg = "REVNANO Error #3 [Stage 3]: Three-staple overlap found. Stop."
		raise ValueError(errormsg)

	# so just overlap of 2 staples at 3 prime, this staple section and another staple section
	# find the identity of that other staple section (can be on the same staple!)
	claimants_3prime.remove( (staple_id, section_id) )
	overlap_id = list( claimants_3prime )[0][0]
	overlap_section_id = list( claimants_3prime )[0][1]

	# go from 3 prime to 5 prime of section
	# see how long overlap is
	r = runlen
	while r > 0 :
		i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - r + 1))

		if (overlap_id, overlap_section_id) in origami["footprints"][i53] :
			num_bases_overlapped += 1
		else :
			break  # this section no longer overlaps the other staple
		r -= 1

	return overlap_id, overlap_section_id, num_bases_overlapped




def fix_overlaps_at_3prime_overshoots_of_staple(origami, staple_id) :
	"""Fixes the overlaps that the passed staple makes at the 3 prime ends of its sections

	Returns:
		The number of bases where this staple overlapped other staples or itself -- which are now fixed (0 if no overlaps)
		-1 if an overlap could not be fixed (due to non-complementarity with scaffold)
		-2 if the 3 prime of the LAST section of the staple overlaps another staple. This is only possible if the last section is placed INcorrectly.
		
		The minus conditions signal that the staple should be omitted.
	"""

	overlapped_bases_fixed = 0

	all_sections = origami["staples_placed"][staple_id]

	for section_id, section in enumerate(all_sections) :

		overlap_id, overlap_section_id, num_bases_overlapped = whois_overlapped_at_section_3prime(origami, staple_id, section_id)

		if overlap_id == None :  # section does not overlap anything at its 3 prime, go to next staple section
			continue

		# this staple section DOES overlap another staple at its 3 prime
		
		# if its the last section, its an error: a correctly placed last section cannot overlap
		if section_id == (len(all_sections) - 1) :
			return -2

		# A) get the scaffold sequence (in direction scaffold 3 prime to 5 prime) 
		# under the overlap at the staple 3 prime
		i53_section_5prime, runlen = section

		scaffold35_A = ""
		i53_A = []
		r = runlen
		while r > (runlen - num_bases_overlapped) :
			i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_section_5prime - r + 1))
			scaffold35_A += origami["scaffold53"][i53]
			i53_A.append(i53)
			r -= 1
		scaffold35_A = scaffold35_A[::-1]

		# B) get the scaffold sequence (in direction scaffold 3 prime to 5 prime) 
		# before the 5 prime of the NEXT section on the staple
		i53_next_section_5prime, runlen_next = all_sections[section_id+1]

		scaffold35_B = ""
		i53_B = []
		for i53raw in range(i53_next_section_5prime + num_bases_overlapped, i53_next_section_5prime, -1) :
			i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53raw))
			scaffold35_B += origami["scaffold53"][i53]
			i53_B.append(i53)

		# scaffold sequences A) and B) must match to fix the overlap
		if scaffold35_A != scaffold35_B :
			return -1

		# scaffold bases in B) must not be claimed by more than 1 other staple, or else a 3-overlap will happen, when this staple is moved there
		for i53 in i53_B :
			if len(origami["footprints"][i53]) > 1 :
				# 3-overlap
				errormsg = "REVNANO Error #3 [Stage 3]: Three-staple overlap found. Stop."
				# Staple <staple_id> moving <num_bases_overlapped> bases from 3' of section <section_id> to 5' of next section causes 3-overlap
				raise ValueError(errormsg)

		# 1) update St dictionary
		I53_INDEX = 0
		RUNLEN_INDEX = 1

		# current section of the staple needs num_bases_overlapped taking off its 3 prime
		# i.e. needs just its runlen changing. Its i53 index remains the same.
		origami["staples_placed"][staple_id][section_id][RUNLEN_INDEX] = runlen - num_bases_overlapped

		# following section of the staple needs num_bases_overlapped added before the existing 5 prime
		# so its needs both its i53 and runlen changing
		origami["staples_placed"][staple_id][section_id+1][I53_INDEX] = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], (i53_next_section_5prime + num_bases_overlapped))
		origami["staples_placed"][staple_id][section_id+1][RUNLEN_INDEX] = runlen_next + num_bases_overlapped		

		# 2) update footprint of staple
		for i53 in i53_A :
			origami["footprints"][i53].remove( (staple_id, section_id) )
		for i53 in i53_B :
			origami["footprints"][i53].add( (staple_id, section_id) )

		overlapped_bases_fixed += num_bases_overlapped

	return overlapped_bases_fixed




def overlap_propagation(origami) :
	"""
	Fix staple overlaps in an iterative cycle, until no further overlaps exist
	"""

	iteration = 1
	overlapped_bases_fixed_last_iteration = 1

	while overlapped_bases_fixed_last_iteration > 0 :
		# loop while overlaps were still fixed on the last iteration
		# (because fixing these overlaps may have caused other ones to pop up)

		print("Iteration %d" % (iteration))
		overlapped_bases_fixed_last_iteration = 0

		# loop through all staples, and fix 3 prime overlaps of each
		for staple_id, _ in enumerate(origami["staples53"]) :

			overlapped_bases_fixed = fix_overlaps_at_3prime_overshoots_of_staple(origami, staple_id)

			if overlapped_bases_fixed == -1 :
				# bases from 3' of a section of this staple could not be added on to before 5 prime of next section
				# because scaffold was not complementary there
				# I take this as a sign this staple is mis-routed, and the staple is removed from the origami
				print("\tOmitting staple %d -- could not move its overlapped bases" % staple_id)
				
				origami["staples_placed"][staple_id] = []
				make_detailed_footprints_dictionary(origami)
				overlapped_bases_fixed = 0

				origami["staple_ids_with_overlap_problems"].append(staple_id)

			elif overlapped_bases_fixed == -2 :
				print("\tOmitting staple %d -- staple 3' overlaps another staple" % staple_id)
				
				origami["staples_placed"][staple_id] = []
				make_detailed_footprints_dictionary(origami)
				overlapped_bases_fixed = 0

				origami["staple_ids_with_overlap_problems"].append(staple_id)

			overlapped_bases_fixed_last_iteration += overlapped_bases_fixed

		print("\tOverlapped bases fixed: %d" % (overlapped_bases_fixed_last_iteration))

		iteration += 1




















#
#
# EXECUTE
#
#

def process(origami) :

	print("-----------------------------------------------------------\n")
	print("REVNANO STAGE 3 (Resolve staple-staple overlaps)           \n")
	print("-----------------------------------------------------------\n")

	# --------------------
	# Check that at least 2 staples have been placed (or else no staples can have overlaps checked)
	# Note: following error cannot happen if Stage 2 always precedes Stage 3, hence commented out
	# --------------------
	#if dna.count_staples_placed(origami) < 2 :
	#	errormsg = "REVNANO Error #: [Stage 3] No staples placed, no staple overlaps exist. Stop."
	#	raise ValueError(errormsg)

	# --------------------
	# Quietly check for pathological overlap cases
	# --------------------
	
	if overlap3plus_exists(origami) :
		errormsg = "REVNANO Error #3 [Stage 3]: Three-staple overlap found. Stop."
		raise ValueError(errormsg)

	# --------------------
	# Build detailed footprints dict, with staple section numbers also
	# --------------------

	make_detailed_footprints_dictionary(origami)

	# --------------------
	# Resolve overlaps
	# --------------------	
	
	print("Resolving staple-staple overlaps...")
	
	overlapped_bases_fixed_total = overlap_propagation(origami)

	print("[Done]\n", flush=True)

	r = dna.remaining_overlaps(origami)
	if r == 0 :
		print("\n ** All staple-staple overlaps resolved ** \n")
	else :
		errormsg = "REVNANO Error #4 [Stage 3]: %d staple-staple overlaps could not be resolved. Stop." % r
		raise ValueError(errormsg)		
	
	# --------------------
	# Stats
	# --------------------
	
	c1 = dna.count_scaffold_bases_claimed(origami)
	staples_placed = dna.count_staples_placed(origami)

	origami["stage3_result"] = (c1, staples_placed)

	

	



