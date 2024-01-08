import sys
# "D:/" is parent folder of REVNANO folder and contactmap, 
# but if you saved each of REVNANO and contactmap in different folder,
# then should use sys.path.append twice.
# if this .py file, REVNANO and contactmap are in same parent directory,
# then ignore this instruction
sys.path.append("D:/")
import os

from REVNANO import naive
from REVNANO import revnano
from contactmap import iohandler as cio
from contactmap import contact2schematic
from contactmap import contactutils

input_folder = "input"
output_folder = "output"

origami_dimension = 2
    
guide_schematic_dimension = 2

MU_MIN = 6

SIGMA  = 4

BETA   = 0.3

folder_path = input_folder  # Specify the folder path here
#txt_files = [file for file in os.listdir(folder_path) if file.endswith(".rev")]
txt_files = [os.path.splitext(file)[0] for file in os.listdir(folder_path) if file.endswith(".rev")]

for origami_name in txt_files:
    
    sequence_filename = input_folder + "/%s.rev" % origami_name
    
    origami = revnano.read_rev_file(origami_name, input_folder + "/", verbose = True)

    if origami == None :
        #raise SystemExit("Stop.")
        continue

    #contactmap, errormsg = revnano.reverse_engineer(origami, MU_MIN, SIGMA, BETA, revnano_deterministic = True, verbose = True)
    contactmap, errormsg = naive.reverse_engineer(origami, MU_MIN)

    if contactmap == None :
        print(errormsg)
        #raise SystemExit("Stop.")
        continue


    cio.write_origami_contactmap(contactmap, output_folder + "/", "%s.csv" % origami_name, madeby="REVNANO")


    contactmap_status = revnano.make_modified_contact_maps_for_guide_schematic_display(origami_name, output_folder + "/", verbose = True)

    if contactmap_status == 0 :
        print("There are no unhybridised 1nt domains on the scaffold.")
    elif contactmap_status == 1 :
        origami_name += "D"    
    elif contactmap_status == 2 :
        raise SystemExit("Hybridised 1nt domains exist on the origami scaffold. Stop.")


    print("Making guide schematic for origami contact map: %s" % origami_name)


    success = contact2schematic.build_guide_schematic(origami_name, output_folder + "/", origami_dimension, guide_schematic_dimension, verbose = True)

    if not success :
        raise SystemExit("Stop.")

    if (origami_dimension == 2) and (guide_schematic_dimension == 2) :
    
        reflect_horizontal = False
        reflect_vertical = False

        rotate_degrees_clockwise = 0 

        contact2schematic.tweak_existing_2D_guide_schematic(origami_name, output_folder + "/", guide_schematic_dimension, xreflect = reflect_horizontal, yreflect = reflect_vertical, rotate_degrees_clockwise = rotate_degrees_clockwise)
