"""REVNANO parameters

Ben Shirt-Ediss, 2020-2023
"""

# Three main affecting staple placement

MU_MIN 				= 6 		# [Set to either 5 or 6]
								# 6 for raster origami (recommended default)
								# 5 for wireframe origami (recommended default)

								# The minimum number of basepairs a hybridised scaffold domain (staple section) is allowed to be
								# The absolute minimum number of nucleotides is 1 less than this number (used for short staples)

SIGMA 				= 4 		# [Set in range 0 <= SIGMA <= 10]
								# 4 for raster origami (recommended default)
								# 2 for wireframe origami (recommended default)

								# Influences number of competitor routes considered during staple routing tree construction
								# higher value = more competitor routes considered
								# sigma too low: shallow routing tree and some correct staple routes will be missing. Over tendency to wrongly classify staples as single route.
								# sigma too high: deep tree, many alternative staple routes, either stage 0 becomes intractable, OR stage 1 cannot narrow down staple routes effectively enough to get into Stage 2

BETA 				= 0.3 		# [Set in range 0 < BETA <= 1]
								# 0.3 for raster origami (recommended default)
								# 0.2 for wireframe origami (recommended default)

								# Controls how much overlap between staple domains is tolerated during staple placement
								# A domain is counted as UNAVAILABLE if this fraction (or greater) of its bases are already claimed
								# note: If BETA = 0, then no staples can be placed.

# If REVNANO is deterministic or not

IS_DETERMINISTIC	= True 		# [Default: True]
								# True = staples polled by ascending staple id (for placement of definite 1-route staples at Stage 0, and placement of >1 route staples at Stage 1)
								# False = staples polled in RANDOM order (note: for some origamis, this still may not lead to non-determinism)

# Fixed parameters

SIGMA_1ROUTE 		= 10 		# [Default: 10]
								# Value of sigma used when building extra deep routing trees to determine if a staple is a definite 1-route staple

SHORT_STAPLE_LEN 	= 14 		# [Default: 14]
								# "Short" staples are this number of nt, or shorter

GAMMA 				= 0.65 		# [Default: 0.65] -- but can be lowered
								# Fraction of all staples required to be placed in Stage 1 for Stage 2 to be able to execute
								# this ensures that the origami mesh is well connected and thus a shortest path through the mesh gives meaningful information
								# (No reverse engineered origami will have less than this fraction of staples placed)

MAX_SUBTREES 		= 10000	 	# [Default: 10000]
								# Maximum number of subtrees allowed in a single staple routing tree (trees growing over this number are killed at Stage 0)

VERSION 			= 1.1   	# REVNANO version

PAPER_DOI 			= "https://doi.org/10.1016/j.csbj.2023.07.011"

