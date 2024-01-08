"""DNA and RNA sequence operations

Ben Shirt-Ediss, April 2021
"""

import random


wc_complements = 	[{"A","T"}, {"C","G"}, {"A","U"}]

bases = {}
bases["DNA"] = ["A","C","G","T"]
bases["RNA"] = ["A","C","G","U"]

complements = {}	# watson-crick complements only
complements["DNA"] = ["T","G","C","A"]
complements["RNA"] = ["U","G","C","A"]

wobble = {}			# extra wobble pairs
wobble["G"] = "U"
wobble["U"] = "G"



# given staple material _ scaffold material, and a staple base, which is a Watson-Crick complement to a scaffold base
# returns the complete set of staple bases that could hybridise with the scaffold base
# Note: assumes only watson-crick complements are used in the DNA or RNA origami design

staple_base_equivalent = {}	

staple_base_equivalent["DNA_DNA"] = {"A" : ["A"],
									 "C" : ["C"],
									 "G" : ["G"],
									 "T" : ["T"]}

staple_base_equivalent["RNA_DNA"] = {"A" : ["A"],		# A on RNA staple section, maps to T on DNA scaffold (WC complement), which only matches with A on RNA staple
									 "C" : ["C", "U"],	# C on RNA staple section, maps to G on DNA scaffold (WC complement), which matches to C or U on RNA staple (latter case is G-U wobble)
									 "G" : ["G"],
									 "U" : ["U"]}

staple_base_equivalent["DNA_RNA"] = {"A" : ["A", "G"],	# A on DNA staple section, maps to U on RNA scaffold (WC complement), which matches to A or G on DNA staple (latter case is G-U wobble)
									 "C" : ["C"],
									 "G" : ["G"],
									 "T" : ["T"]}

staple_base_equivalent["RNA_RNA"] = {"A" : ["A", "G"],
									 "C" : ["C", "U"],
									 "G" : ["G"],
									 "U" : ["U"]}





def bases_are_WC_compliments(base1, base2) :
	"""For DNA, or RNA
	"""

	return {base1,base2} in wc_complements



def is_valid_base(b, material) :

	return b.upper() in bases[material.upper()]



def is_valid_base_set(bset, material) :

	return bset.issubset( set(bases[material.upper()]) )



def is_valid_sequence(s, material) :

	unique_bases = set(s)
	return is_valid_base_set(unique_bases, material.upper())




def complementary_staple_base(scaffold_base, scaffold_material, staples_material, watson_crick = True) :
	"""Given a scaffold base, returns possible complementary staple bases, in the staples material

	If watson_crick = True, Watson-Crick complements are enforced and each base has exactly ONE complement
	"""

	scaffold_base = scaffold_base.upper()

	idx = bases[scaffold_material].index(scaffold_base)

	if watson_crick or (scaffold_base not in ["G","U"]):
		return ( complements[staples_material][idx], )
	else :
		# watson crick = False, and scaffold base is G or U, and so can form wobble pair
		# (GU wobble pairs are present in RNA-RNA, DNA-RNA, RNA-DNA scaffold-staples)
		return ( complements[staples_material][idx], wobble[scaffold_base] )




def random_sequence(material, slen) :
	# note: assumes 4 letter alphabet
	
	seq53 = [ random.choice(bases[material]) for n in range(0, slen) ] 

	return "".join(seq53)



def to_material(sequence, into_material) :

	if into_material == "DNA" :
		return sequence.replace("U", "T")	# if sequence is RNA, will get converted to DNA
	else :
		return sequence.replace("T", "U")	# if sequence is DNA, will get converted to RNA




def linear_circular(scaffold_type, scaffold_len, i53) :
	"""Allows a scaffold sequence -- specified as a linear string in the 5' to 3' direction -- to be used in a linear or a circular way

	Given i53, returns the effective i53 index to use to access the linear string of scaffold
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


		







if __name__ == "__main__" :

	"""
	print( bases_are_WC_compliments("A", "U") )
	print( bases_are_WC_compliments("A", "T") )
	print( bases_are_WC_compliments("A", "A") )
	print( bases_are_WC_compliments("A", "C") )
	print( bases_are_WC_compliments("A", "G") )
	print( bases_are_WC_compliments("G", "C") )
	print( bases_are_WC_compliments("U", "A") )
	quit()
	"""

	# testing complements

	# works OK -- DNA scaffold, DNA staple, G on scaffold will bind to U on staple.
	# but U on staple cannot occur when staple is DNA

	# DNA scaffold
	print("--- DNA SCAFFOLD ---")

	print("DNA scaffold, DNA staple :: watson-crick only")
	print( "A", complementary_staple_base("A", "DNA", "DNA", watson_crick = True) ) 
	print( "C", complementary_staple_base("C", "DNA", "DNA", watson_crick = True) ) 
	print( "G", complementary_staple_base("G", "DNA", "DNA", watson_crick = True) ) 
	print( "T", complementary_staple_base("T", "DNA", "DNA", watson_crick = True) ) 

	print("DNA scaffold, DNA staple :: all")
	print( "A", complementary_staple_base("A", "DNA", "DNA", watson_crick = False) ) 
	print( "C", complementary_staple_base("C", "DNA", "DNA", watson_crick = False) ) 
	print( "G", complementary_staple_base("G", "DNA", "DNA", watson_crick = False) ) 
	print( "T", complementary_staple_base("T", "DNA", "DNA", watson_crick = False) ) 

	print("DNA scaffold, RNA staple :: watson-crick only")
	print( "A", complementary_staple_base("A", "DNA", "RNA", watson_crick = True) ) 
	print( "C", complementary_staple_base("C", "DNA", "RNA", watson_crick = True) ) 
	print( "G", complementary_staple_base("G", "DNA", "RNA", watson_crick = True) ) 
	print( "T", complementary_staple_base("T", "DNA", "RNA", watson_crick = True) ) 

	print("DNA scaffold, RNA staple :: all")
	print( "A", complementary_staple_base("A", "DNA", "RNA", watson_crick = False) ) 
	print( "C", complementary_staple_base("C", "DNA", "RNA", watson_crick = False) ) 
	print( "G", complementary_staple_base("G", "DNA", "RNA", watson_crick = False) ) 
	print( "T", complementary_staple_base("T", "DNA", "RNA", watson_crick = False) ) 

	# RNA 
	print("--- RNA SCAFFOLD ---")

	print("RNA scaffold, RNA staple :: watson-crick only")
	print( "A", complementary_staple_base("A", "RNA", "RNA", watson_crick = True) ) 
	print( "C", complementary_staple_base("C", "RNA", "RNA", watson_crick = True) ) 
	print( "G", complementary_staple_base("G", "RNA", "RNA", watson_crick = True) ) 
	print( "U", complementary_staple_base("U", "RNA", "RNA", watson_crick = True) ) 

	print("RNA scaffold, RNA staple :: all")
	print( "A", complementary_staple_base("A", "RNA", "RNA", watson_crick = False) ) 
	print( "C", complementary_staple_base("C", "RNA", "RNA", watson_crick = False) ) 
	print( "G", complementary_staple_base("G", "RNA", "RNA", watson_crick = False) ) 
	print( "U", complementary_staple_base("U", "RNA", "RNA", watson_crick = False) ) 

	print("RNA scaffold, DNA staple :: watson-crick only")
	print( "A", complementary_staple_base("A", "RNA", "DNA", watson_crick = True) ) 
	print( "C", complementary_staple_base("C", "RNA", "DNA", watson_crick = True) ) 
	print( "G", complementary_staple_base("G", "RNA", "DNA", watson_crick = True) ) 
	print( "U", complementary_staple_base("U", "RNA", "DNA", watson_crick = True) ) 

	print("RNA scaffold, DNA staple :: all")
	print( "A", complementary_staple_base("A", "RNA", "DNA", watson_crick = False) ) 
	print( "C", complementary_staple_base("C", "RNA", "DNA", watson_crick = False) ) 
	print( "G", complementary_staple_base("G", "RNA", "DNA", watson_crick = False) ) 
	print( "U", complementary_staple_base("U", "RNA", "DNA", watson_crick = False) ) 


