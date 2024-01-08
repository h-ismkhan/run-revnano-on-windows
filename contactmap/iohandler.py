"""
I/O operations

Ben Shirt-Ediss
"""

import json



#
#
# READ oxDNA TOPOLOGY AND CONFIGURATION FILES
#
#

def read_oxdna_topology_file(top_filename) :

	ln = 0
	nucleotide_id = 0

	top = {}	# nucleotide id maps to tuple (strand_id base 3'_neighbour_nucleotide_id 5'_neighbour_nucleotide_id)

	try:
		with open("%s" % (top_filename), "r") as read_file:

			for line in read_file:
				
				stripped_line = str(line).strip()
				ln += 1

				if ln == 1 :
					# first line is total nucleotides, total strands -- not significant here
					continue

				# the second line, is nucleotide id 0

				line_data = stripped_line.split(" ")
				line_data[0] = int(line_data[0])	# second item is a string
				line_data[2] = int(line_data[2])
				line_data[3] = int(line_data[3])

				top[nucleotide_id] = tuple( line_data )				

				nucleotide_id += 1
				
				#	Each line represents a single nucleotide on a strand, and says what its adjacent neighbours are
				# 	line number - 2 is the id of the nucleotide
				#
				#	Example:
				#
				#	1 T 6203 1
				#	strand_id 	base_letter 	3'_neighbour_nucleotide_id 		5'_neighbour_nucleotide_id
				#
				# 	3'_neighbour_nucleotide_id and 5'_neighbour_nucleotide_id are -1, when no neighbour exist
				#
				#   Note that the top file does not change, as backbones do not get broken in simulation
				#   (there are no nicking polymerases etc). However, the config file with 
				#	nucleotide positions, direction vectors and velocities does change.

	except(FileNotFoundError) :
		return "oxDNA topology file '%s' does not exist." % (top_filename)

	return top





def read_oxdna_config_file(config_filename) :
	
	ln = 0
	nucleotide_id = 0

	config = {}

	try:
		with open("%s" % (config_filename), "r") as read_file:

			for line in read_file:
				
				stripped_line = str(line).strip()	
				ln += 1

				if ln <= 3 :
					# ignore the first 3 lines of simulation settings
					continue

				line_data = stripped_line.split(" ")[:9]
				line_data = [float(d) for d in line_data]

				#  Each line represents position and velocity properties of a single nucleotide
				#  Here, we are only interested in the positional properties below
				#
				#  rx, ry, rz,    bx, by, bz,   nx, ny, nz
				#
				#  r is the vector representing the point of centre of mass of the whole nucleotide monomer (not the centre of mass of the nucleobase at the end)
				#  b is the vector pointing FROM the backbone of the DNA strand to the base, and appears to already be a unit vector
				#  n is the vector normal to the face of the base, pointing in helix direction 5' to 3'

				config[nucleotide_id] = tuple( line_data )

				nucleotide_id += 1

	except(FileNotFoundError) :
		return "oxDNA config file '%s' does not exist." % (config_filename)

	return config

		





















#
#
# READ SCADNANO JSON
#
#

def read_scadnano_json(json_filename) :

	try:
		with open("%s" % (json_filename), "r") as read_file:
			data = json.load(read_file)
	except(json.decoder.JSONDecodeError) :
		return "scadnano JSON file is not valid JSON."
	except(FileNotFoundError) :
		return "scadnano JSON file does not exist."

	return data























#
#
# READ AND WRITE ORIGAMI CONTACT MAP
#
#


def to_int_list(mylist) :

	mylist = [int(i) for i in mylist] 

	return mylist




def read_origami_contactmap(csv_filename) :

	st = {}   	# staples to scaffold mapping
	sc = {}		# scaffold to staples mapping

	scaffold_circular = False

	try:
		with open("%s" % (csv_filename), "r") as read_file:

			data_region = 0
			waiting_for_data_region = True

			for line in read_file:
				stripped_line = line.strip()		

				if len(stripped_line) == 0 :		# skip blank lines
					continue
				
				if stripped_line.startswith(">") :

					if waiting_for_data_region :
						data_region += 1
						waiting_for_data_region = False
				else :
					waiting_for_data_region = True

					line_data = stripped_line.split(",")

					if data_region == 1 :

						staple_id = int(line_data[0])
						seq = line_data[1] 										# 0 
						num_sec = int(line_data[2])								# 1
						seclen_list = to_int_list(line_data[3].split(";"))		# 2
						dom_list = to_int_list(line_data[4].split(";"))			# 3
						base_list = to_int_list(line_data[5].split(";"))		# 4

						st[staple_id] = {"sequence53" 					: seq,
										 "num_sections"		 			: num_sec,
										 "section_lengths" 				: seclen_list,
										 "section_scaffold_domain_ids" 	: dom_list,
										 "base_i53"						: base_list}

					elif data_region == 2 :
				
						if line_data[0].upper() in ["LINEAR", "CIRCULAR"] :
							scaffold_circular = (line_data[0].upper() == "CIRCULAR")
							continue

						i53 = int(line_data[0]) 
						base_letter = line_data[1]					# 0
						dom_id = int(line_data[2])					# 1
						staple_id = int(line_data[3])				# 2
						staple_section_id = int(line_data[4])		# 3
						staple_base_id = int(line_data[5])			# 4

						sc[i53] = 		{"base_letter"			: base_letter,
										 "scaffold_domain_id"	: dom_id,
										 "staple_id"			: staple_id,
										 "staple_section_id"	: staple_section_id,
										 "staple_base_id"		: staple_base_id}

	# a string is returned on failure
	except(FileNotFoundError) :
		return "Contact map file '%s' does not exist." % (csv_filename)
	except :
		return "Contact map file '%s' contains errors and could not be parsed." % (csv_filename)

	# contact map is a 3-element list. A list not a tuple, so that its mutable, and can have e.g. scaffold_circular attribute changed
	return [st, sc, scaffold_circular]




def tuple_to_delimstr(tup) :

	if tup == None :
		return ""

	output = [str(x) for x in tup]
	return ";".join(output)




def write_origami_contactmap(contactmap, directory, filename, madeby="") :	

	st, sc, scaffold_circular = contactmap

	if not directory.endswith("/") :
		directory += "/"

	f = open("%s%s" % (directory, filename), "w")

	# Part 1: Write how staples are connected to scaffold
	f.write("> %s origami contact map [made by %s]\n" % (filename, madeby))
	f.write("> ----------------------------------------------------------- \n")
	f.write("> STAPLES TO SCAFFOLD MAPPING. Each staple listed 5' --> 3' \n")
	f.write("> column 1: staple id\n")
	f.write("> column 2: staple sequence\n")
	f.write("> column 3: number of staple sections\n")
	f.write("> column 4: staple section lengths (; delimited list)\n")
	f.write("> column 5: scaffold domains respective staple sections hybridise to (; delimited list. Domain 0 at scaffold 5')\n")
	f.write("> column 6: scaffold bases respective staple bases hybridise to (; delimited list. Base 0 at scaffold 5')\n")
	f.write("> ----------------------------------------------------------- \n")
	f.write("> %s,%s,%s,%s,%s,%s\n" % ("staple_id", "sequence", "num_sections", "section_lengths", "section_dom_id", "scaffold_base_ids"))

	for staple_id in sorted(st.keys()) :

		staple = st[staple_id]

		f.write("%d,%s,%d,%s,%s,%s\n" % (	staple_id,
											staple["sequence53"],
											staple["num_sections"],
											tuple_to_delimstr(staple["section_lengths"]),
											tuple_to_delimstr(staple["section_scaffold_domain_ids"]),
											tuple_to_delimstr(staple["base_i53"])
										))

	# Part 2: Write how scaffold is connected to staples
	f.write("> ----------------------------------------------------------- \n")
	f.write("> SCAFFOLD TO STAPLES MAPPING. Scaffold listed 5' --> 3'\n")
	f.write("> column 1: scaffold base id (Base 0 at scaffold 5') \n")
	f.write("> column 2: scaffold base letter\n")
	f.write("> column 3: scaffold domain id (Domain 0 at scaffold 5')\n")
	f.write("> column 4: staple id hybridised here (-1 for no staple)\n")
	f.write("> column 5: staple section id (of above staple) hybridised here\n")
	f.write("> column 6: staple base id (of above staple) hybridised here\n")
	f.write("> ----------------------------------------------------------- \n")
	f.write("> %s,%s,%s,%s,%s,%s\n" % ("scaffold_base_id", "base", "dom_id", "staple_id", "staple_section_id", "staple_base_id"))	
	
	if scaffold_circular :
		f.write("%s\n" % "CIRCULAR")
	else :
		f.write("%s\n" % "LINEAR")

	for i53 in sorted(sc.keys()) :
		base = sc[i53]

		f.write("%d,%s,%d,%d,%d,%d\n" % (	i53,
											base["base_letter"],
											base["scaffold_domain_id"],
											base["staple_id"],
											base["staple_section_id"],
											base["staple_base_id"]
										))

	f.close()




























#
#
# CONTACT MAP TO FASTA SEQUENCE FILE, or REVNANO FILE
#
#

def write_FASTA(contactmap, directory, filename) :

	if not directory.endswith("/") :
		directory += "/"

	st, sc, scaffold_circular = contactmap

	scaffold_type = "LINEAR"
	if scaffold_circular : scaffold_type = "CIRCULAR"

	f = open("%s%s" % (directory, filename), "w")

	scaffold53 = "".join([sc[i53]["base_letter"] for i53 in sorted(sc.keys())])
	f.write(">scaffold_%dnt_%s_5to3prime\n" % (len(scaffold53), scaffold_type))
	f.write("%s\n" % (scaffold53))

	for staple_id in st :
		seq53 = st[staple_id]["sequence53"]
		f.write(">staple%d_%dnt_5to3prime\n" % (staple_id, len(seq53)))
		f.write("%s\n" % (seq53))

	f.close()



def write_REVNANO(contactmap, directory, filename) :

	if not directory.endswith("/") :
		directory += "/"

	st, sc, scaffold_circular = contactmap

	scaffold_type = "LINEAR"
	if scaffold_circular : scaffold_type = "CIRCULAR"

	f = open("%s%s" % (directory, filename), "w")

	scaffold53 = "".join([sc[i53]["base_letter"] for i53 in sorted(sc.keys())])
	f.write("# Scaffold, %d nt, 5' to 3'\n\n" % len(scaffold53))
	f.write("%s" % scaffold53)

	f.write("\n\n%s\n\n" % scaffold_type)

	f.write("# %d Staples, all listed 5' to 3'\n" % len(st))
	f.write("# Staple dangles and interior loopouts enclosed between * pairs\n\n")
	for staple_id in sorted(st.keys()) :

		# important: star * delimiters need marking around staple dangles, and staple interior loopouts, in the staple sequence
		# so REVNANO has a chance at reversing the design

		seq53 = ""

		inhyb = True
		for i, i53 in enumerate(st[staple_id]["base_i53"]) :
			if i53 == -1 and inhyb == True : # just entered into a loopout or dangle
				seq53 += "*"
				inhyb = False

			if i53 != -1 and inhyb == False : # just exited a loop or dangle
				seq53 += "*"
				inhyb = True
			
			seq53 += st[staple_id]["sequence53"][i]

		if inhyb == False :	# if ends unhybrised, add final star
			seq53 += "*"

		f.write("%s\n" % (seq53))

	f.close()

































#
#
# READ ONE-LINE SCAFFOLD SEQUENCE FILE
#
#

def read_sequence(filename) :

	seq = None

	try:
		with open(filename, "r") as read_file:

			for line in read_file:
				seq = line.strip().replace(" ", "")	# sequence is set to string on last line of the file, with internal spaces removed

	except(FileNotFoundError) :
		return None

	return seq.upper()



























#
#
# READ STAPLES CSV FILE
#
#

def read_staples_csv(staples_filename) :
	
	# Note: staples CSV file is one sequence per line
	# Note: * characters are removed from sequences

	staple_list = []

	try:
		with open("%s" % (staples_filename), "r") as read_file:

			for line in read_file:
					
				stripped_line = str(line).strip()
				stripped_line = stripped_line.replace(" ","")
				stripped_line = stripped_line.replace(",","")
				stripped_line = stripped_line.replace(";","")
				stripped_line = stripped_line.replace("*","")

				if len(stripped_line) > 0 :
					staple_list.append(stripped_line.upper())

	except(FileNotFoundError) :
		return "Staples CSV file '%s' does not exist." % (staples_filename)

	return staple_list





























#
#
# WRITE CONTACT MAP DIFF FILE
#
#

def write_contactmap_diff(origami_name1, diff1, origami_name2, diff2, diff_filename) :

	f = open(diff_filename, "w")

	f.write("%s\n" % origami_name1)
	f.write("staple_id, sequence53, num_sections, section_lengths, section_scaffold_domain_ids, base_i53\n")

	for staple_id in sorted(diff1.keys()) :
		sequence53, num_sections, section_lengths, section_scaffold_domain_ids, base_i53 = diff1[staple_id]
		f.write("%d,%s,%d,%r,%r,%r\n" % (staple_id,sequence53, num_sections, tuple_to_delimstr(section_lengths), tuple_to_delimstr(section_scaffold_domain_ids), tuple_to_delimstr(base_i53)))
	
	f.write("%s\n" % origami_name2)
	f.write("staple_id, sequence53, num_sections, section_lengths, section_scaffold_domain_ids, base_i53\n")

	for staple_id in sorted(diff2.keys()) :
		sequence53, num_sections, section_lengths, section_scaffold_domain_ids, base_i53 = diff2[staple_id]
		f.write("%d,%s,%d,%r,%r,%r\n" % (staple_id,sequence53, num_sections, tuple_to_delimstr(section_lengths), tuple_to_delimstr(section_scaffold_domain_ids), tuple_to_delimstr(base_i53)))

	f.close()






















#
#
# WRITE NUCLEOTIDES PAIRED LIST
#
#


def write_bases_paired(bases_paired, oxview_filename, bp_filename) :
	"""This javascript is for pasting into the Browser Console, while oxView is running and is displaying the oxDNA design
	It highlights all nucleotides which are base paired. 
	Tested on Google Chrome. Use Cmd-Option-J to enter Browser Console on Mac.
	"""

	f = open(oxview_filename, "w")
	g = open(bp_filename, "w")

	g.write("scaffold nucleotide id, staple nucleotide id paired with\n")

	nuc_ids = ""

	for scaffold_nucleotide_id in sorted(bases_paired["scaffold_to_staple"].keys()) :
		staple_nucleotide_id = bases_paired["scaffold_to_staple"][scaffold_nucleotide_id]
		g.write("%d,%d\n" % (scaffold_nucleotide_id, staple_nucleotide_id))
		nuc_ids += "%d,%d," % (scaffold_nucleotide_id, staple_nucleotide_id)		

	nuc_ids	= nuc_ids[:-1]

	f.write("clearSelection();\n")
	f.write("api.selectElementIDs([%s]);\n" % (nuc_ids))
	f.write("render();\n")
	f.close()

	g.write("staple nucleotide id, scaffold nucleotide id paired with\n")

	for staple_nucleotide_id in sorted(bases_paired["staple_to_scaffold"].keys()) :
		scaffold_nucleotide_id = bases_paired["staple_to_scaffold"][staple_nucleotide_id]
		g.write("%d,%d\n" % (staple_nucleotide_id, scaffold_nucleotide_id))

	g.close()
































#
#
# READ STAPLE REPLACEMENT FILE
#
#

def read_staple_dangles_file(filename) :

	from_to = {} 		# old_sequence: new_sequence

	try:
		with open("%s" % (filename), "r") as read_file:

			for line in read_file:

				parts = line.split(" ")

				if len(parts) != 2 :
					return "Error: all lines in file must be of format: old_sequence new_sequence."
					
				old_sequence = parts[0].strip().upper()
				new_sequence = parts[1].strip().upper()

				from_to[old_sequence] = new_sequence

	except(FileNotFoundError) :
		return "Staple replacement file '%s' does not exist." % (filename)


	# check new_sequence is always (1) longer than old_sequence and (2) contains it as a subset

	dangles = {}		# old_sequence: (5' dangle to add, 3' dangle to add)

	for old_sequence in from_to :
		new_sequence = from_to[old_sequence]

		if (len(new_sequence) > len(old_sequence)) and (old_sequence in new_sequence) :
			# extract 5' end, and 3' end dangles to add

			start_idx = new_sequence.find(old_sequence)
			end_idx = start_idx + len(old_sequence)

			dangle5 = new_sequence[0:start_idx]
			dangle3 = new_sequence[end_idx:]

			dangles[old_sequence] = (dangle5, dangle3)
		else :
			return "Error: for each line, new_sequence should contain old_sequence, and also be longer than it."

	return dangles



































#
#
# REVNANO FILE AUTO-STARRING
#
#

def mark_staple_sequence_polyT_islands(seq) :

	# 1. find the distinct polyT islands in the sequence and store in record
	# (islands must be TTTT or longer)

	record = list("." * len(seq))

	seq_list = list(seq)

	push = 0

	while True :

		try :
			i = "".join(seq_list).index("TTTT")
		except ValueError :
			break

		# TTT is still in the string
		record[i+push] = "T"
		record[i+1+push] = "T"
		record[i+2+push] = "T"
		record[i+3+push] = "T"

		seq_list = seq_list[1:]

		push += 1

	# 2. mark the original sequence with * surrounding polyT islands

	seq_marked = ""
	phase = "looking for T region"

	for i, base in enumerate(list(seq)) :

		if (record[i] == "T") and (phase == "looking for T region") :
			seq_marked += "*"
			phase = "in T region"
			
		if (record[i] == ".") and (phase == "in T region") :
			seq_marked += "*"
			phase = "looking for T region"

		seq_marked += base

	if record[-1] == "T" :
		seq_marked += "*"

	return seq_marked





def REVNANO_auto_star(filename) :
	"""Takes an existing (valid) REVNANO sequence input file, and marks all polyT islands on staples with * delimiters
	"""

	scaffold_type = "LINEAR"
	scaffold53 = ""
	staples53 = []

	scaffold_found = False	

	with open(filename) as csv_file:	

		for line_num, line in enumerate(csv_file):

			line = line.strip().replace(" ", "").upper()	# sequences with internal spaces have these spaces removed

			if len(line) == 0 or line[0] == "#" :	# skip blank lines and comments
				continue

			if line in ["LINEAR", "CIRCULAR"] :		# scaffold type can be anywhere (and appear multiple times)
				scaffold_type = line
				continue

			if scaffold_found == False :			# catch first sequence as scaffold sequence
				scaffold53 = line
				scaffold_found = True
				continue

			staples53.append(line)					# otherwise append to staples


	# now re-output the file (devoid of comments)

	f = open(filename, "w")
	f.write("%s\n" % scaffold53)
	f.write("%s\n" % scaffold_type)

	for seq in staples53 :

		seq = seq.replace("*","")		# get rid of any existing * delimiters
		
		f.write("%s\n" % mark_staple_sequence_polyT_islands(seq))

	f.close()







