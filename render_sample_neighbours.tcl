mol load gro [lindex $argv 0]

# define colour dictionaries for dengue and flu residue type of interest

#turn off the axes:
axes location Off 

#switch to Orthographic view:
display projection Orthographic

#switch to white background:
color Display Background white

# set all atoms to the approximate vdw radius in MARTINI FF
set all_atoms [atomselect top "all"]
$all_atoms set radius 2.3

# as a background, colour all particles in a non-imposing silver / gray colour
mol modselect 0 0 all
mol modstyle 0 0 QuickSurf 1.0 0.5 1.0 1.0
mol modmaterial 0 0 Ghost
mol modcolor 0 0 ColorID 6 

set num_resid_values [lindex $argv 2]
set virus_type [lindex $argv [expr {$num_resid_values + 6}]]

if {$virus_type == "dengue"} {
	set colour_dict [dict create DUPC "11" PPCE "0" PPCS "9" DPPE "7"]
} else {
	set colour_dict [dict create CHOL "7" FORS "5" PPCH "1" DOPX "11"]
}

set argv_position 2
for {set repnum 1} {$repnum < [expr {$num_resid_values * 3}]} {incr repnum 3} {
	set argv_position [expr {$argv_position + 1}]
	set resid_value [lindex $argv $argv_position]
	puts "resid_value: $resid_value"
	set sel [atomselect top "resid $resid_value"]
	set current_resname [lindex [$sel get resname] 0]
	puts "current_resname: $current_resname"
	set current_colour_ID [dict get $colour_dict $current_resname]
	mol addrep 0

	if {$virus_type == "dengue"} {
		mol modselect $repnum 0 (resid $resid_value) and name PO4
} else {
		mol modselect $repnum 0 (resid $resid_value) and ((resname PPCH and name PO4) or (resname CHOL and name ROH) or (resname FORS and name AM2) or (resname DOPX and name PO4))
}
	mol modstyle $repnum 0 VDW 1.0 20.0
	mol modmaterial $repnum 0 AOShiny
	mol modcolor $repnum 0 ColorID $current_colour_ID 

	# display the rest of the (non-headgroup) particles in the highlight lipids in gray
	set second_rep [expr {$repnum + 1}]
	mol addrep 0
	if {$virus_type == "dengue"} {
		mol modselect $second_rep 0 (resid $resid_value) and not name PO4 NC3 NH3
} else {
		mol modselect $second_rep 0 (resid $resid_value) and not name PO4 NC3 NH3 BC B2 B3 INV ROH AM2 B1 B4
}
	mol modstyle $second_rep 0 VDW 1.0 20.0
	#mol modstyle $second_rep 0 DynamicBonds 6.0 0.3 6.0
	mol modmaterial $second_rep 0 AOShiny
	mol modcolor $second_rep 0 ColorID 6

	# connect the CG particles
	set third_rep [expr {$repnum + 2}]
	mol addrep 0
	mol modselect $third_rep 0 (resid $resid_value) and not name NC3 NH3 BC B2 B3 INV AM2 B1 B4
	mol modstyle $third_rep 0 DynamicBonds 6.0 0.9 9.0
	mol modmaterial $third_rep 0 AOShiny
	mol modcolor $third_rep 0 ColorID 16

}

#color lipids orange:
#mol addrep 0
#mol modselect 1 0 (resname POPC or resname PPCE or resname DPPE or resname CER or resname DUPC or resname DOPS or resname PPCS) and z < 300
#mol modstyle 1 0 VDW 1.0 20.0
#mol modmaterial 1 0 AOShiny
#mol modcolor 1 0 ColorID 3

#load another copy of the coord file and use for XC view
#mol load gro dengue_final_snapshot.gro
#set sel [atomselect top "all"]
#$sel moveby {500 0 0}
#
##colour proteins blue:
#mol modselect 0 1 (index 0 to 221039) and z < 330
#mol modstyle 0 1 VDW 1.0 20.0
#mol modmaterial 0 1 AOShiny
#mol modcolor 0 1 ColorID 10 
#
##color lipids orange:
#mol addrep 1
#mol modselect 1 1 (resname POPC or resname PPCE or resname DPPE or resname CER or resname DUPC or resname DOPS or resname PPCS) and z < 300
#mol modstyle 1 1 VDW 1.0 20.0
#mol modmaterial 1 1 AOShiny
#mol modcolor 1 1 ColorID 3

#color FORS green:
#mol addrep 0
#mol modselect 2 0 not (index 0 to 344389) and (resname FORS)
#mol modstyle 2 0 VDW 1.0 20.0
#mol modmaterial 2 0 AOShiny
#mol modcolor 2 0 ColorID 7

#draw some scale bars
#draw color black
#draw cylinder {260 110 0} {660 110 0} radius 20.0

#zoom:
#scale by 2.5
#rotate the view:
# read in the rotation values from argv

rotate x by [lindex $argv [expr {$argv_position + 1}]]
rotate y by [lindex $argv [expr {$argv_position + 2}]]
rotate z by [lindex $argv [expr {$argv_position + 3}]]
#tranlate the view a bit (in screen units):
#translate by -0.2 0.0 0


#activate ambient occlusion lighting and shadows:
display shadows on
display ambientocclusion on

#use Tachyon to generate a .dat file, which can then be rendered externally:
render TachyonInternal [lindex $argv 1]
#render Tachyon dengue_final_snapshot.dat
