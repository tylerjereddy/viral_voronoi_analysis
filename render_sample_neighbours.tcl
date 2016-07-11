mol load gro [lindex $argv 0]

#turn off the axes:
axes location Off 

#switch to Orthographic view:
display projection Orthographic

#switch to white background:
color Display Background white

#colour proteins blue:
mol modselect 0 0 all
mol modstyle 0 0 VDW 1.0 20.0
mol modmaterial 0 0 AOShiny
mol modcolor 0 0 ColorID 10 

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
#rotate x by 50
#rotate y by -60
#rotate z by 194
#tranlate the view a bit (in screen units):
#translate by -0.2 0.0 0


#activate ambient occlusion lighting and shadows:
display shadows on
display ambientocclusion on

#use Tachyon to generate a .dat file, which can then be rendered externally:
render TachyonInternal test.tga
#render Tachyon dengue_final_snapshot.dat
