""" REVNANO file I/O operations

Ben Shirt-Ediss, 2020, 2023
"""

import sys
sys.path.append("../")

from REVNANO import dna





#
#
# LOAD .rev CSV FILE
#
#

def load_input_csv(filename) :
	"""Read a CSV file with origami sequences, and parse scaffold and staples sequences and scaffold type
	- # signal comment lines
	- First sequence encountered is the scaffold
	- Following sequences are the staples
	- "linear" or "circular" specifies type of scaffold. Scaffold assumed linear if it is not stated
	- Non-sequences which are not comments are not allowed (for general tidiness)
	- * pairs can appear in staple sequences marking dangles or interior loopouts
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

			# treat line as a sequence
			sequence_type = dna.sequence_type(line)

			if sequence_type == "* DELIMITER MISTAKE" :
				errormsg = "Initialisation Error: Sequence at line %d of %s contains a mistake with * delimiters. Stop." % (line_num+1, filename)
				raise ValueError(errormsg)				

			elif sequence_type == "UT MISTAKE" :
				errormsg = "Initialisation Error: Sequence at line %d of %s contains both U and T bases. Individual sequences must have either DNA or RNA bases, not mix both. Stop." % (line_num+1, filename)
				raise ValueError(errormsg)

			elif sequence_type == "NUCLEIC ACID" :

				if scaffold_found == False :			# catch first sequence as scaffold sequence
					scaffold53 = line
					scaffold_found = True
					continue

				staples53.append(line)					# otherwise append to staples

			else :
				errormsg = "Initialisation Error: Illegal text at line %d of %s. Stop." % (line_num+1, filename)
				raise ValueError(errormsg)
	
	return (scaffold_type, scaffold53, staples53)




















