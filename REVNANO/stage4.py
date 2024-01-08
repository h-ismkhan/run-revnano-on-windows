"""REVNANO : Stage 4 : Re-compose Staples and Make Contact Map 

Ben Shirt-Ediss, 2023
"""

import sys
sys.path.append("../")

from REVNANO import dna

from contactmap import contactutils






#
#
# CONVERT origami DICT TO CONTACT MAP FORMAT
# 
#

def convert_to_contactmap(origami) :

	st = {}   	# staples to scaffold mapping
	sc = {}		# scaffold to staples mapping

	St = origami["staples_placed"]

	# Part 0: Remove staples that REVNANO was unable to map
	# (in this case, staple ids are not continuous)
	all_staple_ids = list(St.keys())

	for staple_id in all_staple_ids :
		if len(St[staple_id]) == 0 :
			del St[staple_id]

	# Part 1: Staples to scaffold mapping
	for staple_id in St :

		base_id = 0 			# the base number on the staple
		section_lengths = []
		i53_list = []

		staple_seq = origami["staples53"][staple_id]

		# go through all sections of this staple
		all_sections = St[staple_id]  			# [[i53, runlen],[i53, runlen],[i53, runlen]]

		for section_id, section in enumerate(all_sections) :
			
			i53_section_5prime, runlen = section
			section_lengths.append(runlen)

			for i53 in range(i53_section_5prime, i53_section_5prime-runlen, -1) :

				i53 = dna.linear_circular(origami["scaffold_type"], origami["scaffold_nt"], i53)

				if i53 < 0 :
					# staple base NOT hyb with scaffold, this section is free
					i53_list.append(-1)
				else :
					# staple base IS hyb with scaffold	
					i53_list.append(i53)

					base = {}
					base["base_letter"] 		= origami["scaffold53"][i53] 	# index 0
					base["scaffold_domain_id"] 	= None							# index 1 # domain ids are decided below
					base["staple_id"] 			= staple_id 					# index 2
					base["staple_section_id"]	= section_id 					# index 3
					base["staple_base_id"]		= base_id						# index 4

					sc[i53] = base
																	
				base_id += 1

		staple = {}
		staple["sequence53"] 					= staple_seq
		staple["num_sections"] 					= len(all_sections)
		staple["section_lengths"] 				= tuple(section_lengths)
		staple["section_scaffold_domain_ids"] 	= None	# add the scaffold domain ids hybridised, when they are decided below
		staple["base_i53"] 						= tuple(i53_list)

		st[staple_id] = staple

	# Part 2: Scaffold to staples mapping
	# add scaffold bases which are not hybridised to staples
	for i53, _ in enumerate(origami["scaffold53"]) :
		if sc.get(i53, -1) == -1 :

			base = {}
			base["base_letter"] 		= origami["scaffold53"][i53] 	# index 0
			base["scaffold_domain_id"] 	= None		# index 1
			base["staple_id"] 			= -1 		# index 2
			base["staple_section_id"]	= -1  		# index 3
			base["staple_base_id"]		= -1		# index 4

			sc[i53] = base

	# add domain numbers to sc
	domain_id = 0

	staple_id0 = sc[0]["staple_id"]
	section_id0 = sc[0]["staple_section_id"]
	last = (staple_id0, section_id0)

	mapback = {}

	for i53 in sorted(sc.keys()) :

		staple_id = sc[i53]["staple_id"]
		section_id = sc[i53]["staple_section_id"]

		# if staple section changed, we are on new domain
		if (staple_id, section_id) != last :	
			domain_id += 1
			last = (staple_id, section_id)

		# domain id of first and last domain is the same, if 
		# scaffold is circular, and they are connected to same staple section		
		sc[i53]["scaffold_domain_id"] = domain_id

		mapback[(staple_id, section_id)] = domain_id


	# if the first and last domains on the scaffold hybridise the same staple (or are not connected to staples)
	# and the scaffold is circular, number them as domain 0
	if origami["scaffold_type"] == "CIRCULAR" :
		for i53 in range(len(sc)-1, 0, -1) :

			staple_id = sc[i53]["staple_id"]
			section_id = sc[i53]["staple_section_id"]		

			if (staple_id, section_id) == (staple_id0, section_id0) :
				sc[i53]["scaffold_domain_id"] = 0 
				mapback[(staple_id, section_id)] = 0
			else :
				break


	# Part 3: Now scaffold domain numbers are set, go through staples one last time and add them
	for staple_id in St :

		domains_list = []
		all_sections = St[staple_id]  			# [[i53, runlen],[i53, runlen],[i53, runlen]]
												# i53 is -1 when section does not bind scaffold
		
		for section_id, section in enumerate(all_sections) :
	
			i53_section_5prime, runlen = section

			if i53_section_5prime < 0 :
				# this staple section is not hybridised to any scaffold domain
				domains_list.append(-1)
			else :
				domains_list.append(mapback[(staple_id, section_id)])

		st[staple_id]["section_scaffold_domain_ids"] = tuple(domains_list)
	
	circular = (origami["scaffold_type"] == "CIRCULAR")

	return st, sc, circular



























#
# 
# RE-COMPOSE STAPLES WITH LOOPOUTS AND DANGLES
#
#

def recompose_staples(origami, contactmap) :
	# note: contactmap argument is updated by reference

	st, _, _ = contactmap

	# 1. re-compose staples with loopouts

	for staple_id in origami["extra_staples"] :

		extra_staples_list = origami["extra_staples"][staple_id]
		loopout_seq_list = origami["loopouts"][staple_id]

		# first check that all substaples of the loopout staple have been placed on the origami
		substaples = set( [staple_id] + extra_staples_list ) 
		allstaples = set( st.keys() )

		if not substaples.issubset(allstaples) :
			# some of the substaples have not been placed on the origami
			# this entire staple with loopouts hence CANNOT be placed, and all substaples in the complete loopout staple must be deleted from contact map
			contactutils.delete_staples(contactmap, list( substaples ))

			print("Staple %d with loopouts could not be recomposed, as some individual staple sections are not placed. Staple deleted" % staple_id)
		else :
			# OK: re-construct this staple with loopouts
			contactutils.recompose_staple_with_loopouts(contactmap, staple_id, extra_staples_list, loopout_seq_list)

	# 2. add 5' and 3' single stranded dangles back to staples (contactmap updated by reference)

	for staple_id, dangle5 in origami["dangles5"].items() :
		if staple_id in st :	# the staple with 5' dangle may not be placed, so check first
			contactutils.add_staple_end_dangles(contactmap, staple_id, dangle5, "")

	for staple_id, dangle3 in origami["dangles3"].items() :
		if staple_id in st :	# the staple with 3' dangle may not be placed, so check first
			contactutils.add_staple_end_dangles(contactmap, staple_id, "", dangle3)
























#
#
# EXECUTE
#
#


def process(origami) :

	print("-----------------------------------------------------------\n")
	print("REVNANO STAGE 4 (Re-compose staples and make contact map)  \n")
	print("-----------------------------------------------------------\n")

	# --------------------
	# Make origami contact map
	# --------------------

	print("Creating origami contact map...")

	contactmap = convert_to_contactmap(origami)

	print("[Done]\n", flush=True)

	# --------------------
	# Re-compose staples with dangles and loopouts
	# --------------------

	print("Re-composing staples with loopouts and dangles...")

	recompose_staples(origami, contactmap)

	print("[Done]\n", flush=True)		

	# --------------------
	# Final Reverse Engineering Stats
	# --------------------

	origami["stage4_result"] = dna.final_stats(contactmap)		# total scaffold bases hybridised, total physical staples placed

	return contactmap




