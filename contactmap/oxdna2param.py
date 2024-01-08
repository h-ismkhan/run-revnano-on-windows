"""Calculates the optimal range for the distance and angle parameters, for detecting base pairs in a particular oxDNA config file.
Displays a heatmap of results.

This procedure can sometimes be necessary, as oxDNA can be quite a "wild" format.
Different tools have difference idiosyncracies in how they export to oxDNA format: 
	- some (like tacoxDNA) compress base pair distances when insertions are present, and lengthen bp distances when deletions are present
		others (like scadnano) don't
	- tacoxDNA outputs helices as much closer together than does scadnano. There may also be differences in nucleotide spacings.

The detection parameters in oxdna2contact were derived by running this script on many
origamis produced by both scadnano and tacoxDNA tools, and seeing what parameter region
gives reliable base pair detection in ALL cases.

Ben Shirt-Ediss, November 2021
"""

import sys
sys.path.append("../")

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt 
import numpy
import os
import pickle

from contactmap import iohandler
from contactmap import oxdna2contact


# ----- Heat Map Parameters ----------

# Note: these standard parameters work with most oxDNA files

MAX_DETECTION_SPHERE_RADIUS = 3.0 		# radius of detection sphere used around each scaffold nucleotide
										# to get local staple nucleotides which could be pase pair candidates
										# this is the upper limit of the heatmap x-axis

THETA_MAX = 35.0

GRID_SIZE = 21
# ------------------------------------










#
#
# CALCULATE OPTIMAL REGION OF BASE PAIR DETECTION PARAMETERS, FOR oxDNA CONFIG FILE
#
#

def stats_base_pairings(N, config, euclidean_cutoff, angle_max) :
	"""
	For the parameter pair given, returned is:
		(1) the number of monogamous scaffold nucleotides that are paired to a "monogamous" staple nucleotide
				(the scaffold base, pairs only to that staple base, and vice versa)
		(2) the number of non-monogamous scaffold nucleotides (scaffold base has multiple staple base partners)
		(3) the number of monogamous scaffold nucleotides that are paired to a "non-monogamous" staple nucleotide (staple base has multiple scaffold base partners)

	The idea is to find the parameter pairs maximising (1) while leaving (2) zero and (3) zero

	(1) and (2) and (3) are returned, rather than returning just (1) alone because
	if (2) or (3) is non-zero, then the base pair detector in oxdna2contact will not work properly.
	"""

	scaffold_to = {}
	staple_to = {}
	
	cheat_staple = set()	
	
	cheat = 0

	# each time a scaffold base pairs uniquely to a staple base, record the staple base
	for scaffold_nucleotide_id in N :

		candidates = 0
		candidate_id = None

		sc_r = config[scaffold_nucleotide_id][0:3]
		sc_b = config[scaffold_nucleotide_id][3:6]

		for staple_nucleotide_id in N[scaffold_nucleotide_id] :

			st_r = config[staple_nucleotide_id][0:3]
			st_b = config[staple_nucleotide_id][3:6]

			if oxdna2contact.is_candidate_base_pair(sc_r, sc_b, st_r, st_b, euclidean_cutoff, angle_max) :
				candidates += 1
				candidate_id = staple_nucleotide_id

		if candidates == 1 :
			# scaffold nucleotide is monogamous ... but is the staple nucleotide its paired to also monogamous?
			scaffold_to[scaffold_nucleotide_id] = candidate_id

			if candidate_id in staple_to :
				# staple not monogamous, its already paired with someone
				staple_to[candidate_id].append(scaffold_nucleotide_id)
				cheat_staple.add(candidate_id)
			else :
				# staple is monogamous (for the moment)
				staple_to[candidate_id] = [scaffold_nucleotide_id]
		elif candidates > 1 :
			# scaffold nucleotide is not monogamous
			cheat += 1
			
		# otherwise, scaffold base is not paired with anything

	# remove scaffold nucleotides which have "cheating" non-monogamous staple nucleotide partners
	monog_cheat = 0

	for staple_nucleotide_id in cheat_staple :
		for scaffold_nucleotide_id in staple_to[staple_nucleotide_id] :
			del scaffold_to[scaffold_nucleotide_id]
			monog_cheat += 1

	monog_monog = len(scaffold_to)

	return monog_monog, cheat, monog_cheat




def optimal_parameter_region(N, config, origami_name) :

	# 1. set parameter grid range
	theta_list = list(numpy.linspace(0,THETA_MAX,GRID_SIZE))		
	d_list = list(numpy.linspace(0,MAX_DETECTION_SPHERE_RADIUS,GRID_SIZE))					
	
	# 2. data matrices 1, 2a and 2b -- details of base pairing
	datamatrix1 = [] 		# monog-monog
	datamatrix2a = [] 		# cheat
	datamatrix2b = [] 		# monog-cheat
	optimal = 0 			# (the maximal value of monog-monog encountered, when cheat and monog-cheat are both zero)
							# not known until the end

	for theta_max in theta_list :	# rows y
		row1 = []; row2a = []; row2b = []

		print("  Theta %2.2f" % (theta_max))

		for d_max in d_list :		# columns x

			monog_monog, cheat, monog_cheat = stats_base_pairings(N, config, d_max, theta_max)
			
			row1.append( monog_monog )
			row2a.append( cheat )
			row2b.append( monog_cheat )

			if (cheat == 0) and (monog_cheat == 0) and monog_monog > optimal :
				optimal = monog_monog

		datamatrix1.append(row1); datamatrix2a.append(row2a); datamatrix2b.append(row2b)

	# 3. data matrix 3 -- optimal parameter region
	datamatrix3 = []
	optimal_region_exists = False

	for i, row1 in enumerate(datamatrix1) :
		
		row2a = datamatrix2a[i]; row2b = datamatrix2b[i]; row3 = []

		for j, r1 in enumerate(row1) :
			if (r1 == optimal) and (row2a[j] == 0) and (row2b[j] == 0) :
				row3.append(1)
				optimal_region_exists = True
			else :
				row3.append(0)

		datamatrix3.append(row3)

	if optimal_region_exists :
		print("%d\tMonogamous scaffold-staple base pairs in optimal parameter region" % (optimal))	
		
		optimal_region = {}
		optimal_region["datamatrix"] = datamatrix3
		optimal_region["theta_list"] = theta_list
		optimal_region["d_list"] = d_list

		pickle.dump( optimal_region, open( "_assets/%s.p" % (origami_name), "wb" ) )	
	else :
		print("(No optimal parameter region exists)")

	# 5. heatmap ticks
	x_ticks = list(range(0, GRID_SIZE))
	x_tick_labels = ["%2.2f" % item for item in d_list]
	y_ticks = list(range(0, GRID_SIZE))
	y_tick_labels = ["%2.2f" % item for item in theta_list]
	
	# heatmap 1: monogamous-monogamous scaffold-staple base pairs
	plt.figure(figsize=(8,8))
	ax = plt.subplot(2,2,1)
	plt.imshow(datamatrix1, interpolation='nearest', cmap="viridis")
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_tick_labels, rotation=45)
	ax.set_yticks(y_ticks)
	ax.set_yticklabels(y_tick_labels)
	ax.tick_params(labelsize = 2.0) # font size of tick labels
	ax.set_xlabel("d_max", fontdict = {'fontsize' : 6})
	ax.set_ylabel("theta_max", fontdict = {'fontsize' : 6})
	ax.set_title("Scaffold (Monog) - Staple (Monog) ", fontdict = {'fontsize' : 6})
	cbar = plt.colorbar(orientation="vertical", shrink=0.5)
	cbar.ax.tick_params(labelsize=6) 

	# heatmap 2: monogamous-cheat scaffold-staple base pairs
	ax = plt.subplot(2,2,2)
	plt.imshow(datamatrix2b, interpolation='nearest', cmap="viridis")
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_tick_labels, rotation=45)
	ax.set_yticks(y_ticks)
	ax.set_yticklabels(y_tick_labels)
	ax.tick_params(labelsize = 2.0) # font size of tick labels
	ax.set_xlabel("d_max", fontdict = {'fontsize' : 6})
	ax.set_ylabel("theta_max", fontdict = {'fontsize' : 6})
	ax.set_title("Scaffold (Monog) - Staple (Cheat)", fontdict = {'fontsize' : 6})
	cbar = plt.colorbar(orientation="vertical", shrink=0.5)
	cbar.ax.tick_params(labelsize=6) 

	# heatmap 3: cheat scaffold base pairs
	ax = plt.subplot(2,2,3)
	plt.imshow(datamatrix2a, interpolation='nearest', cmap="viridis")
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_tick_labels, rotation=45)
	ax.set_yticks(y_ticks)
	ax.set_yticklabels(y_tick_labels)
	ax.tick_params(labelsize = 2.0) # font size of tick labels
	ax.set_xlabel("d_max", fontdict = {'fontsize' : 6})
	ax.set_ylabel("theta_max", fontdict = {'fontsize' : 6})
	ax.set_title("Scaffold (Cheat)", fontdict = {'fontsize' : 6})
	cbar = plt.colorbar(orientation="vertical", shrink=0.5)
	cbar.ax.tick_params(labelsize=6) 	

	# heatmap 4: optimal region where only monog-monog base pairs exist, and there is the maximum number of them
	ax = plt.subplot(2,2,4)
	plt.imshow(datamatrix3, interpolation='nearest', cmap="Greens")
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_tick_labels, rotation=45)
	ax.set_yticks(y_ticks)
	ax.set_yticklabels(y_tick_labels)
	ax.tick_params(labelsize = 2.0) # font size of tick labels
	ax.set_xlabel("d_max", fontdict = {'fontsize' : 6})
	ax.set_ylabel("theta_max", fontdict = {'fontsize' : 6})
	ax.set_title("Optimal Parameter Region (%d bp)" % (optimal), fontdict = {'fontsize' : 6})

	#plt.show()
	plt.savefig("_assets/%s.png" % (origami_name), dpi = 600, bbox_inches='tight')




















#
#
# CALCULATE CONSENSUS OPTIMAL REGION
#
#

def consensus_optimal_parameter_region() :
	"""Note: a consensus region won't exist for some sets of origamis. 
	"""

	datamatrix = None
	d_list = None
	theta_list = None

	num_origamis = 0

	# open all files ending .p in _assets/ dir, loads the data matrix,
	# and makes a data matrix which is the COMMON optimal parameter region of all data matrices

	for file in os.listdir("_assets/"):
		if file.endswith(".p"):

			optimal_region = pickle.load( open( "_assets/%s" % (file), "rb" ) )
			num_origamis += 1

			if datamatrix == None :
				datamatrix = optimal_region["datamatrix"]
				d_list = optimal_region["d_list"]
				theta_list = optimal_region["theta_list"]
				continue

			# check the axes of this data matrix are the same as the first data matrix added
			if optimal_region["d_list"] == d_list and optimal_region["theta_list"] == theta_list :

				for j, row in enumerate(optimal_region["datamatrix"]) :
					row_prev = datamatrix[j]
					# only keep elements "1" in row of datamatrix, that
					# are 1 both in datamatrix, and in the new datamatrix
					datamatrix[j] = list((a and b) for a,b in zip(row_prev,row))

	if datamatrix != None :

		x_ticks = list(range(0, len(d_list) ))
		x_tick_labels = ["%2.2f" % item for item in d_list]
		y_ticks = list(range(0, len(theta_list) ))
		y_tick_labels = ["%2.2f" % item for item in theta_list]
	
		plt.figure(1)
		ax = plt.subplot(1,1,1)
		plt.imshow(datamatrix, interpolation='nearest', cmap="Greens")
		ax.set_xticks(x_ticks)
		ax.set_xticklabels(x_tick_labels, rotation=45)
		ax.set_yticks(y_ticks)
		ax.set_yticklabels(y_tick_labels)
		ax.tick_params(labelsize = 3.0) # font size of tick labels
		ax.set_xlabel("d_max", fontdict = {'fontsize' : 8})
		ax.set_ylabel("theta_max", fontdict = {'fontsize' : 8})
		ax.set_title("Consensus Optimal Parameter Region (%d Origamis)" % (num_origamis), fontdict = {'fontsize' : 8})

		plt.show()




















#
#
# VISUAL DEMO TO HELP UNDERSTAND BASE PAIR DETECTION
#
#

def base_pair_detection_visual_debug(N, scaffold_nucleotide_id, euclidean_cutoff, angle_max) :
	"""Draw base pair detection, for a single scaffold nucleotide
	"""

	print("Scaffold nucleotide: %d\n" % (scaffold_nucleotide_id))

	print("Angle 1 = angle between b vector of scaffold nucleotide, and connecting vector draw from scaffold nucleotide to staple nucleotide")
	print("Angle 2 = angle between b vector of staple nucleotide, and connecting vector draw from staple nucleotide to scaffold nucleotide\n")

	plt.figure()
	ax = plt.axes(projection='3d')

	sc_r = config[scaffold_nucleotide_id][0:3]
	sc_b = config[scaffold_nucleotide_id][3:6]

	# 1. draw scaffold base
	ax.scatter3D(sc_r[0],sc_r[1],sc_r[2],color="xkcd:blue")
	# draw b vector of scaffold base, for convenience starting at centre of scaffold base
	x = [sc_r[0], sc_r[0]+sc_b[0]*0.3]
	y = [sc_r[1], sc_r[1]+sc_b[1]*0.3]
	z = [sc_r[2], sc_r[2]+sc_b[2]*0.3]
	ax.plot3D(x,y,z,"k")
	# label scaffold base
	ax.text(sc_r[0],sc_r[1],sc_r[2],  '%s' % (str(scaffold_nucleotide_id)), size=8, zorder=1, color='k') 

	# 2. draw local staple base neighbours
	for staple_nucleotide_id in N[scaffold_nucleotide_id] :
			
		st_r = config[staple_nucleotide_id][0:3]
		st_b = config[staple_nucleotide_id][3:6]

		# draw staple base
		ax.scatter3D(st_r[0],st_r[1],st_r[2],color="xkcd:red")
		# draw b vector of staple base, for convenience starting at centre of staple base
		x = [st_r[0], st_r[0]+st_b[0]*0.3]
		y = [st_r[1], st_r[1]+st_b[1]*0.3]
		z = [st_r[2], st_r[2]+st_b[2]*0.3]
		ax.plot3D(x,y,z,"k")
		# label staple base
		ax.text(st_r[0],st_r[1],st_r[2],  '%s' % (str(staple_nucleotide_id)), size=8, zorder=1, color='k') 

		# draw the vector, from the scaffold base to the staple base
		sc_st_vector = oxdna2contact.vector_between_points_3d(sc_r, st_r)
		x = [sc_r[0], sc_r[0]+sc_st_vector[0]]
		y = [sc_r[1], sc_r[1]+sc_st_vector[1]]
		z = [sc_r[2], sc_r[2]+sc_st_vector[2]]		
		ax.plot3D(x,y,z,"xkcd:orange")

		# (vector from staple base to scaffold base is just the reverse)
		st_sc_vector = oxdna2contact.vector_between_points_3d(st_r, sc_r)

		# output stats
		print("Staple nucleotide: %d" % (staple_nucleotide_id))
		print("-- Distance ", oxdna2contact.euclidean_distance_3d(sc_r, st_r) )
		print("-- Angle1 ", oxdna2contact.vector_angle_degrees(sc_b, sc_st_vector) )
		print("-- Angle2 ", oxdna2contact.vector_angle_degrees(st_b, st_sc_vector) )
		print("-- Paired? ", oxdna2contact.is_candidate_base_pair(sc_r, sc_b, st_r, st_b, euclidean_cutoff, angle_max) )


	plt.title("Scaffold Nucleotide %d: Staple Neighbours" % (scaffold_nucleotide_id))
	ax.view_init(azim=67,elev=34) # set azimuth and elevation, to best show the 3D plot
	plt.show()


	













#
#
# EXECUTE
#
#

if __name__ == "__main__" :

	print("\n------------------------------------------------")
	print("oxDNA --> Optimal Base Pair Detection Parameters")
	print("------------------------------------------------\n")

	if len(sys.argv) != 2 :
		print("Usage: python3 oxdna2param.py <origami_name>\n")
		quit()

	origami_name = sys.argv[1]

	if origami_name == "consensus" :

		print("Finding consensus optimal parameter region...")
		consensus_optimal_parameter_region()
		print("[Done]\n")

	else :
		# perform parameter analysis for the oxDNA config file of the origami specified

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

		print("Converting oxDNA topology file...")
		topology = oxdna2contact.oxdna_topology_file_to_useful_format(top)
		print("[Done]")

		print("Building scaffold nucleotide neighbour index...")
		N = oxdna2contact.build_scaffold_nucleotide_neighbour_index(topology, config, MAX_DETECTION_SPHERE_RADIUS)
		print("[Done]")

		# get a certain point of parameter grid
		#monog_monog, cheat, monog_cheat = stats_base_pairings(N, config, 2.0, 10.0)
		#print(monog_monog)
		#print(monog_cheat)
		#print(cheat)
		#quit()

		# for a single chosen scaffold nucleotide, see visual demo of base pair detection,
		base_pair_detection_visual_debug(N, 0, 2.0, 10.0)	
		quit()

		print("Building parameter heatmaps and saving to _assets/ directory...")
		optimal_parameter_region(N, config, origami_name)
		print("[Done]\n")

	
	