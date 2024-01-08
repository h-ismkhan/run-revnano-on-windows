"""Uses the scadnano python API, to convert cadnano JSON files into scadnano SC files
(and vice versa)

Ben Shirt-Ediss
"""

import scadnano as sc
import sys
import os


if __name__ == "__main__" :

	print("\n--------------------------------")
	print("Cadnano <--> Scadnano Conversion")
	print("--------------------------------\n")
	
	if len(sys.argv) != 2 :
		print("Usage: python3 cadnano2scadnano.py <origami_name.sc> or <origami_name.json>\n")
		quit()

	filename = sys.argv[1]

	if not os.path.exists("_assets/%s" % (filename)) :
		print("Error: file '%s' does not exist in the _assets/ directory.\n" % (filename))
		quit()


	if filename.endswith(".sc") :
		# convert scadnano to cadnano
		
		print("Converting %s to cadnano..." % (filename))
		
		design = sc.Design.from_scadnano_file("_assets/%s" % filename)
		
		try :
			design.export_cadnano_v2(directory="_assets", filename=filename.replace(".sc",""))
		except ValueError as ve:
			print("scadnano Error: %s\n" % str(ve))
			quit()


	elif filename.endswith(".json") :
		# convert cadnano to scadnano

		print("Converting %s to s-cadnano..." % (filename))

		design = sc.Design.from_cadnano_v2(directory="_assets", filename=filename)
		design.write_scadnano_file(directory="_assets", filename=filename.replace(".json",""))

	print("[Done]\n")



