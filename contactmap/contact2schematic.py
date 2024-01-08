"""Attempts to re-construct an approximate origami "guiding schematic" (2D or 3D) 
from the topological origami contact map. This "guiding schematic" helps the manual
entry of the reverse engineered origami design into a design tool of choice.

The NetworkX implementation of the 2D or 3D Kamada Kawai spring embedder at
networkx.drawing.layout.kamada_kawai_layout() is used to
get an initial geometric layout of graph nodes and edges. This particular implementation
of the Kamada-Kawai algorithm seems to be a good one.

Then, html/javascript engine https://github.com/vasturiano/3d-force-graph
is used to apply further forces to the graph nodes in a physical simulaton
and display the graph as an interactive 3D web page (allows the graph to be zoomed, inspected, nodes can be re-positioned etc)

Ben Shirt-Ediss, September-October 2021, Feb 2023
"""

import sys
sys.path.append("../")


import networkx
import os
import pickle
import math
import numpy

from contactmap import iohandler
from contactmap import contact2dlgraph
from contactmap import contact2ambig
from contactmap import contactutils









#
#
# KAMADA-KAWAI SPRING FORCE GRAPH LAYOUT (in 2 or 3 dimensions)
#
#

def kk_spring_embedder(G, dim) :
	"""NetworkX seems to contain many implementations of KK, but this one is by far the most full-featured:
		
	https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.kamada_kawai_layout.html
	
	Initial node positions are not specified, so KK will use a circular layout of nodes and move nodes from that
	however, this seems to lead to an acceptable optimum layout in most cases (especially in 3D)

	For origamis, which are essentially lattice graphs, KK is much better than Fruchtermann Reingold layout
	For an origami lattice graph, the optimal distance between two nodes is correlated to the shortest path between those nodes
	"""

	# 1. do raw KK layout

	node_coords = networkx.drawing.layout.kamada_kawai_layout(G, weight="schematiclen", dim=dim)
		# note: weight is name of edge attribute containing the ideal edge length

		# returns: node_id (its i53 index) : (tuple of its coordinates; x,y or x,y,z)

	""" Notes:
	- kk layout uses schematiclen as link weighting, and gives a geometric layout in arbitrary coordinates (I call this "KK space")
	- edge distances are not schematiclen, but the ratio of actual sizes of two edges in KK space should approximate the ratio of their schematiclen's
	- below we scale up arbitrary kk coordinates, so that length of edges approaches their schematiclen (then axes are, approximately, in bases)
	- this makes layout of the points by Force Graph / 3D Force Graph much more predictable
	- (some of the coordinates are likely negative as some of the original kk coordinates are likely negative, but this is no problem)
	"""

	# 2. scale KK coordinates returned
	scaled_node_coords = {}

	# get the ends (i53 indexes) of the first graph edge
	# (the edge may be a domain or crossover, does not matter)
	u, v = list(G.edges())[0]
	schematiclen = G.edges[u, v]["schematiclen"]
	
	# get the arbitrary KK coordinates of each end of the edge, 
	# and work out the length of the edge vector in KK space (vectors can be 2D or 3D)
	# efficient euclidean distance can be computed with numpy (not scipy):
	# https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy

	point1 = numpy.array(node_coords[u])
	point2 = numpy.array(node_coords[v])

	kklen = numpy.linalg.norm(point1 - point2)

	# scale all points in KK space
	scale_factor = schematiclen / kklen

	for node_id in node_coords :
		scaled_node_coords[node_id] = tuple([c * scale_factor for c in node_coords[node_id]])

	return scaled_node_coords		


	


















#
#
# LABEL DOMAIN-LEVEL GRAPH READY FOR DISPLAY AS HTML PAGE
#
#


# these settings work well with most origamis

THIN_LINE_WIDTH = 1
THICK_LINE_WIDTH = 5
NODE_SIZE = 5
ENGINE_COOLDOWN_MS = 10000		# JS force-layout engine time before no more updates = 10 seconds

xkcd__teal_blue = "#01889f" 		# scaffold nodes
xkcd__bright_yellow = "#fffd01" 	# scaffold ssdna
xkcd__steel = "#738595" 			# scaffold dsdna
xkcd__tangerine = "#ff9408" 		# scaffold nicks
xkcd__hot_pink = "#ff028d" 			# 3 prime and 5 prime nodes on linear scaffold
xkcd__beige = "#e6daa6"				# de-selected staple/scaffold crossovers
xkcd__green = "#15b01a" 			# non-ambiguous
xkcd__red = "#e50000" 				# ambiguous


def calc_charge_force(G) :
	"""Small origamis (10 graph nodes) are given -100 force, so they spring apart and you can see all links
	Larger origamis (from 130 graph nodes) are given -5 force. (An M13 origami typically has 1200 nodes on the domain-level graph)
	Force is linearly interpolated between these bounds and constant outside them
	"""

	# linear zone
	x1 = 10; y1 = -100 		# x in nodes, y is force value
	x2 = 130; y2 = -5

	m = (y2-y1) / (x2-x1)

	# c = y-mx for one point
	c = y2 - (m * x2)

	x = G.number_of_nodes()
	
	force = (m * x) + c

	# outside linear zone
	if force < -100 :
		force = -100
	elif force > -5 :
		force = -5

	return force




def node_is_ambiguous(node_id, A) :
	# node_id is assumed to be within A
	return A[node_id] != (0,0)




def node_ambiguous_colour(node_id, A) :

	if node_is_ambiguous(node_id, A):
		return xkcd__red
	else :
		return xkcd__green




def edge_ambiguous_colour(u, v, A) :

	if node_is_ambiguous(u, A) or node_is_ambiguous(v, A) :
		return xkcd__red		# ambiguous edge: one or both ends of the edge are ambiguous
	else :
		return xkcd__green




def truncate_long_seq(seq) :
	# sequences over about 40nt don't play nice in tool tip info boxes: they are truncated to fit

	if len(seq) > 43 :
		return seq[0:20] + "..." + seq[-20:]	# clip to 43nt, listing only first 20 and last 20 bases
	
	return seq



def add_interactive_display_attributes_to_G(G, A, colours) :
	"""Adds extra node and edge attributes to allow interactive display of the graph G
	in three different view modes: scaffold view, staples view and ambig view

	A is the staple split ambiguity dictionary for the ambiguous view
	"""

	# 1. EDGES

	# edge tool tips
	# (same in all SCAFFOLD, STAPLE and AMBIGUOUS CROSSOVER views)
	for u,v,a in G.edges(data=True) :

		if a["type"] == "ss_domain" :
			G.edges[u, v]["tooltip"] = "<p class='tt'>%d nt single stranded scaffold</p><p class='acgt_sc'>3'-%s-5'</p><p class='tt'>[scaffold domain %d]</p>" % (a["runlen"], truncate_long_seq(a["scaffold_sequence53"][::-1]), a["scaffold_domain_id"])
		elif a["type"] == "ds_domain" :
			G.edges[u, v]["tooltip"] = "<p class='tt'>%d bp helix<br />staple %d, section %d</p><p class='acgt'>5'-%s-3'</p><p class='acgt_sc'>3'-%s-5'</p><p class='tt'>[scaffold domain %d]</p>" % (a["runlen"], a["staple_id"], a["section_id"], truncate_long_seq(a["staple_sequence53"]), truncate_long_seq(a["scaffold_sequence53"][::-1]), a["scaffold_domain_id"])
		elif a["type"] == "ghost_nick" :
			G.edges[u, v]["tooltip"] = "<p class='tt'>scaffold nick</p>"
		elif a["type"] == "crossover" :
			G.edges[u, v]["tooltip"] = "<p class='tt'>staple %d crossover</p>" % (a["staple_id"])
		elif a["type"] == "loopout" :
			G.edges[u, v]["tooltip"] = "<p class='tt'>%d nt loopout<br />staple %d, section %d</p><p class='acgt'>5'-%s-3'</p>" % (a["runlen"], a["staple_id"], a["section_id"], truncate_long_seq(a["staple_sequence53"]))
		elif a["type"] == "dangle" :
			G.edges[u, v]["tooltip"] = "<p class='tt'>%d nt dangling end<br />staple %d, section %d</p><p class='acgt'>5'-%s-3'</p>" % (a["runlen"], a["staple_id"], a["section_id"], truncate_long_seq(a["staple_sequence53"]))

	# edge widths
	for u,v,a in G.edges(data=True) :
		if a["type"] == "ss_domain" :
			G.edges[u, v]["width_sc"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_st"] = THIN_LINE_WIDTH
			G.edges[u, v]["width_ac"] = THIN_LINE_WIDTH
		elif a["type"] == "ds_domain" :
			G.edges[u, v]["width_sc"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_st"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_ac"] = THICK_LINE_WIDTH
		elif a["type"] == "ghost_nick" :
			G.edges[u, v]["width_sc"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_st"] = THIN_LINE_WIDTH
			G.edges[u, v]["width_ac"] = THIN_LINE_WIDTH
		elif a["type"] == "crossover" :
			G.edges[u, v]["width_sc"] = THIN_LINE_WIDTH
			G.edges[u, v]["width_st"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_ac"] = THICK_LINE_WIDTH
		elif a["type"] == "loopout" :			
			G.edges[u, v]["width_sc"] = THIN_LINE_WIDTH
			G.edges[u, v]["width_st"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_ac"] = THICK_LINE_WIDTH
		elif a["type"] == "dangle" :
			G.edges[u, v]["width_sc"] = THIN_LINE_WIDTH
			G.edges[u, v]["width_st"] = THICK_LINE_WIDTH
			G.edges[u, v]["width_ac"] = THICK_LINE_WIDTH
		
		G.edges[u, v]["width"] = G.edges[u, v]["width_sc"] 	# scaffold view is default

	# edge colours
	for u,v,a in G.edges(data=True) :

		colour_id = G.edges[u, v].get("staple_id", 0) % len(colours)

		if a["type"] == "ss_domain" :
			G.edges[u, v]["colour_sc"] = xkcd__bright_yellow
			G.edges[u, v]["colour_st"] = xkcd__beige
			G.edges[u, v]["colour_ac"] = edge_ambiguous_colour(u, v, A)
		elif a["type"] == "ds_domain" :
			G.edges[u, v]["colour_sc"] = xkcd__steel
			G.edges[u, v]["colour_st"] = colours[colour_id]
			G.edges[u, v]["colour_ac"] = edge_ambiguous_colour(u, v, A)
		elif a["type"] == "ghost_nick" :
			G.edges[u, v]["colour_sc"] = xkcd__tangerine
			G.edges[u, v]["colour_st"] = xkcd__beige
			G.edges[u, v]["colour_ac"] = edge_ambiguous_colour(u, v, A)
		elif a["type"] == "crossover" :
			G.edges[u, v]["colour_sc"] = xkcd__beige
			G.edges[u, v]["colour_st"] = colours[colour_id]
			G.edges[u, v]["colour_ac"] = edge_ambiguous_colour(u, v, A)
		elif a["type"] == "loopout" :			
			G.edges[u, v]["colour_sc"] = xkcd__beige
			G.edges[u, v]["colour_st"] = colours[colour_id]
			G.edges[u, v]["colour_ac"] = edge_ambiguous_colour(u, v, A)
		elif a["type"] == "dangle" :						
			G.edges[u, v]["colour_sc"] = xkcd__beige
			G.edges[u, v]["colour_st"] = colours[colour_id]
			G.edges[u, v]["colour_ac"] = edge_ambiguous_colour(u, v, A)

		G.edges[u, v]["colour"] = G.edges[u, v]["colour_sc"] 	# scaffold view is default

	# 2. NODES

	# node tool tips
	# (same in all SCAFFOLD, STAPLE and AMBIG views)
	for i53 in G.nodes() :

		# ambiguous tooltip info to add IF node is ambiguous
		amb_text = ""
		if node_is_ambiguous(i53, A) :
			changeminus, changeplus = A[i53]
			amb_text = "<br/>{ambig %d to +%d}" % (changeminus, changeplus)	

		if G.nodes[i53].get("type", -1) != -1 :
			if G.nodes[i53]["type"] == "5prime" :
				G.nodes[i53]["tooltip"] = "<p class='tt'>5' SCAFFOLD (base 0)</p>%s" % (amb_text)
			elif G.nodes[i53]["type"] == "3prime" :
				G.nodes[i53]["tooltip"] = "<p class='tt'>3' SCAFFOLD (base %d)</p>%s" % (i53, amb_text)
			elif G.nodes[i53]["type"] == "dangle":
				G.nodes[i53]["tooltip"] = "<p class='tt'>%s</p>" % (amb_text)	# no label on nodes at end of staple dangles
		else :
			G.nodes[i53]["tooltip"] = "<p class='tt'>scaffold base %d%s</p>" % (i53, amb_text)				
						
	# node colours
	for i53 in G.nodes() :
		if G.nodes[i53].get("type", -1) != -1 and G.nodes[i53]["type"] in ["5prime", "3prime"] :
			# node is 5' or 3'
			G.nodes[i53]["colour_sc"] = xkcd__hot_pink
			G.nodes[i53]["colour_st"] = xkcd__hot_pink
			G.nodes[i53]["colour_ac"] = node_ambiguous_colour(i53, A)
		else :
			# node is any other node
			colour = xkcd__beige
			if G.nodes[i53].get("staple_id", -1) != -1 :
				# node is on a staple, and should be coloured as that staple
				colour_id = G.nodes[i53]["staple_id"] % len(colours)
				colour = colours[colour_id]

			G.nodes[i53]["colour_sc"] = xkcd__teal_blue
			G.nodes[i53]["colour_st"] = colour
			G.nodes[i53]["colour_ac"] = node_ambiguous_colour(i53, A)

		G.nodes[i53]["colour"] = G.nodes[i53]["colour_sc"] 		# scaffold view is default




def get_attributes_for_force_graph_HTML(G, node_coords, kk_dim, space_dim) :
	"""Converts the graph G into a set of attributes that
	force-graph/3d-force-graph uses to display the graph

	force-graph: 	https://github.com/vasturiano/force-graph
	3d-force-graph: https://github.com/vasturiano/3d-force-graph
	"""

	attr = {}

	# 1. nodes and links of graph, as javascript object
	# graph nodes are set at initial positions calculated by KK spring embedder
	nodes = "nodes: [ "
	for i53 in G.nodes() :

		tt =	G.nodes[i53]["tooltip"]
		
		c_sc =	G.nodes[i53]["colour_sc"]
		c_st =	G.nodes[i53]["colour_st"]
		c_ac =	G.nodes[i53]["colour_ac"]
		c = 	G.nodes[i53]["colour"]
		
		x = 0; y = 0; z = 0;
		if kk_dim == 2 :
			x, y = node_coords[i53]
			nodes += "{id: %d, x: %2.4f, y: %2.4f, c_sc: \"%s\", c_st: \"%s\", c_ac: \"%s\", c: \"%s\", tt: \"%s\"}," % \
			 	(i53, x, y, c_sc, c_st, c_ac, c, tt)
		elif kk_dim == 3 :
			x, y, z = node_coords[i53]
			nodes += "{id: %d, x: %2.4f, y: %2.4f, z: %2.4f, c_sc: \"%s\", c_st: \"%s\", c_ac: \"%s\", c: \"%s\", tt: \"%s\"}," % \
				(i53, x, y, z, c_sc, c_st, c_ac, c, tt)

	nodes = nodes[:-1] + " ]"

	links = "links: [ "
	for u,v,a in G.edges(data=True):		# u and v are node ids, a is a dict of link attributes
		
		tt = 	a["tooltip"]

		w_sc = 	a["width_sc"]
		w_st =	a["width_st"]
		w_ac =	a["width_ac"]
		w =		a["width"]

		c_sc =	a["colour_sc"]
		c_st =	a["colour_st"]
		c_ac =	a["colour_ac"]
		c = 	a["colour"]

		links += "{source: %d, target: %d, type: \"%s\", schematiclen: %d, w_sc: %d, w_st: %d, w_ac: %d, w: %d, c_sc: \"%s\", c_st: \"%s\", c_ac: \"%s\", c: \"%s\", tt: \"%s\"}," % \
					(a["source"], a["target"], a["type"], a["schematiclen"], w_sc, w_st, w_ac, w, c_sc, c_st, c_ac, c, tt)
	
	links = links[:-1] + " ]"

	attr["nodes"] = nodes
	attr["links"] = links

	# 2. differences in javascript for 2D and 3D force graph display
	attr["js_include"] = "<script src=\"https://unpkg.com/force-graph\"></script>"
	attr["js_force_graph"] = ""
	attr["js_num_dimensions"] = ""
	attr["js_linkopacity"] = ""

	if space_dim == 3 :	# redefine, if graph is 3D
		attr["js_include"] = "<script src=\"https://unpkg.com/3d-force-graph\"></script>"
		attr["js_force_graph"] = "3D"
		attr["js_num_dimensions"] = "\n.numDimensions(%d)" % (kk_dim)
		attr["js_linkopacity"] = "\n.linkOpacity(0.8)"

	# 3. force and cooldown attributes
	attr["node_rel_size"] = NODE_SIZE
	attr["cooldown_ms"] = ENGINE_COOLDOWN_MS
	attr["charge_force_strength"] = calc_charge_force(G)

	return attr






























#
#
# WRITE GUIDE SCHEMATIC HTML/JS PAGE
# 		(based on engine: https://github.com/vasturiano/3d-force-graph)
#
#

def write_force_graph_HTML(attr, origami_name, assets_dir) :

	f = open("%s%s.html" % (assets_dir, origami_name), "w")		

	html = """<!DOCTYPE html>
<html lang="en">
	<head>
		<meta charset="utf-8">
		<title>%s</title>
		<style>
			body { 
				margin: 0;
				background-color: #fff;
				font-size: 10pt;
			} 

			.info {
				padding: 5px 0 0 0;
				color: #000;
				font-family: Helvetica, sans-serif;
			}

			.tt {
				padding: 2px;
				margin: 0;
				font-family: Helvetica, sans-serif;
				background: #222;
			}

			.acgt {
				padding: 2px;
				margin: 0;
				font-family: Courier New, Courier, Lucida Sans Typewriter, Lucida Typewriter, monospace;
				background: #000;
			}

			.acgt_sc {
				padding: 2px;
				margin: 0;
				font-family: Courier New, Courier, Lucida Sans Typewriter, Lucida Typewriter, monospace;
				background: #000;
				color: #95d0fc;
			}	

			#currentview {
				position: absolute;
				top: 0;
				left: 10px;
			}

			#currentforces {
				position: absolute;
				top: 0;
				right: 10px;
			}			

			#controls {
				width: 100%%;
				text-align: center;
			}

			.inline {
				display: inline-block;
			}
		</style>
		%s
	</head>
	<body>
		<div id="currentview" class="info">
			Scaffold View
		</div>

		<div id="currentforces" class="info">
			Forces ON
		</div>

		<div id="controls">
			<div class="inline"><button id="buttonview" onclick="toggle_view()">Toggle View</button></div>
			<div class="inline"><button onclick="recentre()">Re-centre</button></div>
			<div class="inline"><button id="buttonforces" onclick="toggle_forces()">Turn Forces OFF</button></div>
		</div>

		<div id="origami-domain-level-graph"></div>
		
		<script>

			// definition of origami graph
			const gData = {
				%s,
				%s
			};

			const Graph = ForceGraph%s()
				(document.getElementById('origami-domain-level-graph'))%s
					.graphData(gData)
					.backgroundColor("#ffffff")
					.nodeRelSize(%d)
					.nodeColor(node => node.c)
					.nodeLabel(node => node.tt)
					.linkWidth(link => link.w)
					.linkColor(link => link.c)%s
					.linkLabel(link => link.tt)
					.linkDirectionalArrowLength(10)
					.linkDirectionalArrowRelPos(0.5)
					.cooldownTime(%d);
					
			Graph.d3Force('link').distance(link => link.schematiclen);
				// note: the default function for link strength (rigidity in range 0 to 1) below works well
				// https://github.com/vasturiano/d3-force-3d#link_distance

			Graph.d3Force('charge').strength(%d);

			function recentre() {
				Graph.zoomToFit(1000);   
			}

			var view = 0;	

			function toggle_view() {
				if (view == 0) {
					// currently scaffold view, go to staples view

					gData.nodes.forEach(node => {
						node.c = node.c_st;
					});
					gData.links.forEach(link => {
						link.w = link.w_st;
						link.c = link.c_st;

						if (link.type == "ds_domain") {
							// reverse arrows on scaffold-hybridised staple sections 
							// link.source returns a javascript node object
							// and thus its easiest just to swap source and target node objects to reverse the link
							node_source = link.source;
							link.source = link.target;
							link.target = node_source;
						}
					});
					view = 1;
					document.querySelector('#currentview').innerText = "Staples View";
					//document.querySelector('#buttonview').innerText = "Next: Ambig View";
				}
				else if (view == 1) {
					// currently staples view, go to ambig view
					
					gData.nodes.forEach(node => {
						node.c = node.c_ac;
					});
					gData.links.forEach(link => {
						link.w = link.w_ac;
						link.c = link.c_ac;		
					});
					view = 2;
					document.querySelector('#currentview').innerText = "Sequence-Ambiguous Junction View";
					//document.querySelector('#buttonview').innerText = "Next: Scaffold View";
				}
				else if (view == 2) {
					// currently ambig view, go to scaffold view

					gData.nodes.forEach(node => {
						node.c = node.c_sc;
					});
					gData.links.forEach(link => {
						link.w = link.w_sc;
						link.c = link.c_sc;

						if (link.type == "ds_domain") {
							// reverse back the arrows on scaffold-hybridised staple sections 
							node_source = link.source;
							link.source = link.target;
							link.target = node_source;
						}								
					});					
					view = 0;
					document.querySelector('#currentview').innerText = "Scaffold Routing View";					
					//document.querySelector('#buttonview').innerText = "Next: Staples View";
				}

				Graph.nodeColor(Graph.nodeColor());
				Graph.linkColor(Graph.linkColor());
			}

			var forces_on = 1;

			function toggle_forces() {
				if (forces_on) {
					// switch forces off
					Graph.cooldownTicks(0);
					forces_on = 0;
					document.querySelector('#currentforces').innerText = "Forces OFF";
					document.querySelector('#buttonforces').innerText = "Turn Forces ON";
				}
				else {
					// turn forces back on (..after the next node is dragged, see issue: https://github.com/vasturiano/3d-force-graph/issues/475)
					Graph.cooldownTicks(Infinity);
					forces_on = 1;
					document.querySelector('#currentforces').innerText = "Forces ON";
					document.querySelector('#buttonforces').innerText = "Turn Forces OFF";
				}
			}

		</script>
	</body>
</html>""" % (origami_name, \
							attr["js_include"], \
							attr["nodes"], attr["links"], \
							attr["js_force_graph"], attr["js_num_dimensions"], \
							attr["node_rel_size"], \
							attr["js_linkopacity"], \
							attr["cooldown_ms"], \
							attr["charge_force_strength"])

	f.write(html)
	f.close()



































#
#
# contact2schematic API
#
#

def build_guide_schematic(origami_name, assets_dir, kk_dim, space_dim, verbose = False) :

	if verbose == False :
		sys.stdout = open(os.devnull, 'w')

	if not assets_dir.endswith("/") :
		assets_dir += "/"

	# check dimensions are correct
	if (kk_dim not in [2, 3]) :
		print("Error: Origami spring layout must be in 2 or 3 dimensions.")
		if verbose == False : sys.stdout = sys.__stdout__
		return False

	if (space_dim not in [2, 3]) :
		print("Error: Origami must be displayed in 2 or 3 dimensional space.")
		if verbose == False : sys.stdout = sys.__stdout__
		return False

	if space_dim < kk_dim :
		print("Error: Cannot display %dD origami in %dD space.\n" % (kk_dim, space_dim))
		if verbose == False : sys.stdout = sys.__stdout__
		return False

	# load contact map
	contactmap = iohandler.read_origami_contactmap("%s%s.csv" % (assets_dir, origami_name))

	if isinstance(contactmap, str) :
		print("Error: %s\n" % (contactmap))
		if verbose == False : sys.stdout = sys.__stdout__
		return False

	# check no problematic 1nt scaffold domains are present in the contact map
	I = contactutils.scaffold_1nt_domain_locations(contactmap)
	if len(I) > 0 :
		print("Error. Could not build the domain-level origami graph.\n1nt scaffold domains found at scaffold bases %r.\n1nt domains cannot be represented in the origami domain-level graph.\n" % I)
		print("Try to remove the 1nt domains by e.g:")
		print("python3 contactutils.py del1nt origami_name")
		print("python3 contactutils.py rs1nt origami_name\n")
		
		if verbose == False : sys.stdout = sys.__stdout__
		return False

	print("Generating domain-level graph from contact map...")
	G = contact2dlgraph.origami_dl_graph(contactmap)
	print("[Done]")
	
	if contact2dlgraph.schematiclen_corrected(G) == True :
		print("Warning: Some edges in the domain-level graph have display length of 0 (indicating a single base domain). These have been adjusted to have a small positive display length, to stop the KK spring embedder from crashing.\n")

	print("Building staple split ambiguity information...")
	A = contact2ambig.build_ambig_dict(G)
	print("[Done]")
	
	print("Updating domain-level graph for display...")

	colours = [	"#e96721", 
				"#f39c2e",
				"#f9b92b",
				"#6cb556",
				"#19a586",
				"#5ab9a0",
				"#1db2ca",
				"#3873b7",
				"#d05c9b",
				"#b04490",
				"#e8575d",
				"#e64835" ]
 
	add_interactive_display_attributes_to_G(G, A, colours)
	print("[Done]")

	print("Running Kamada-Kawai spring embedder in %d dimensions..." % (kk_dim))
	node_coords = kk_spring_embedder(G, kk_dim)
	print("[Done]")

	# if 2D, save node-coords (they may need to be tweaked later)
	if kk_dim == 2 :
		pickle.dump( (G, node_coords), open( "%s%s.layout" % (assets_dir, origami_name), "wb" ) )

	print("Generating force-graph HTML page to display origami schematic in %d dimensions..." % (space_dim))
	attr = get_attributes_for_force_graph_HTML(G, node_coords, kk_dim, space_dim) 
	write_force_graph_HTML(attr, origami_name, assets_dir)
	print("[Done]\n")

	if verbose == False : sys.stdout = sys.__stdout__
	return True






def rotate(origin, point, rotate_degrees_clockwise):
    """ Rotate point around origin. Modified from:
    https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
    """

    angle = math.radians(rotate_degrees_clockwise)

    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    
    return qx, qy






def tweak_existing_2D_guide_schematic(origami_name, assets_dir, space_dim, xreflect = False, yreflect = False, rotate_degrees_clockwise = 0) :
	""" When a 2D guide schematic exists, this method allows reflecting and rotating the coordinates
	to make the schematic orientation more suitable

	Reflections are applied before rotations

	yreflect = True
		reflect the shape about the y-axis
	xreflect = True
		reflect the shape about the x-axis

	rotate_degrees_clockwise = 90
		rotate the shape 90 degrees clockwise, about the centre of the shape

	Returns False on failure

	This function produces no verbose output
	"""

	# 1. check layout pickle file exists for this origami (if it does, it is a 2D origami)
	layout_file = "%s%s.layout" % (assets_dir, origami_name)

	if not os.path.exists(layout_file) :
		return False

	# 2. get origami graph and node coordinates
	G, node_coords = pickle.load( open( layout_file, "rb" ) )

		# node_coords   node_id : (x,y)

	# 3. apply reflections
	minx = 1e12
	miny = 1e12
	node_coords2 = {}

	for node_id in node_coords :

		x,y = node_coords[node_id]

		if yreflect :
			x = x * -1
		if xreflect :
			y = y * -1

		# find bottom left coord of shape bounding box (this will be the rotation origin point, later)
		if x < minx :
			minx = x
		if y < miny :
			miny = y

		node_coords2[node_id] = (x, y)

	# 4. apply rotation
	if rotate_degrees_clockwise in [0, 360] :
		node_coords3 = node_coords2
	else :
		node_coords3 = {}
		origin = (minx, miny)

		for node_id in node_coords2 :
		
			point = node_coords2[node_id]
			node_coords3[node_id] = rotate(origin, point, rotate_degrees_clockwise)

	# 5. re-display the 2D origami in 2D or 3D space
	attr = get_attributes_for_force_graph_HTML(G, node_coords3, 2, space_dim) 
	write_force_graph_HTML(attr, origami_name, assets_dir)






























#
#
# contact2schematic CLI
#
#


if __name__ == "__main__" :

	# note: the CLI version cannot reflect and rotate 2D structures, like the API version can

	print("")
	print("----------------------------------------------------")
	print("Reconstruct Origami Guide Schematic from Contact Map")
	print("----------------------------------------------------")
	print("")

	# check args
	if len(sys.argv) < 4 :
		print("Example usage: python3 contact2schematic.py rectangle <layout dimension> <space dimension>\n")
		quit()

	try :
		origami_name = sys.argv[1]
		kk_dim = int(sys.argv[2]) 		# kk embed dimension
		space_dim = int(sys.argv[3])	# space_dim dimension to display kk solution
	except ValueError :
		print("Example usage: python3 contact2schematic.py rectangle 2 2\n")
		quit()

	if not os.path.exists("_assets/%s.csv" % (origami_name)) :
		print("Error: contact map file '%s.csv' does not exist in the _assets/ directory.\n" % (origami_name))
		quit()		

	print("Contact map:   %s.csv\n" % (origami_name))

	success = build_guide_schematic(origami_name, "_assets/", kk_dim, space_dim, verbose = True)

	if success :
		print("Open %s in a web browser.\n" % ("_assets/%s.html" % origami_name))

	


