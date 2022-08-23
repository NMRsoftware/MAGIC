# MAGIC Act MAGIC View
Prepare files for automated methyl assignment 

# Instalation 
Download scripts into the Poky modleues folder:
  /Applications/poky_mac/poky.app/Contents/Resources/modules/poky

To create the two letter command codes (MA, MV) to open the modules in sparky edit the poky_site.py file and add the following lines under the assignment_menu:

    ('MA', 'MAGIC-Act',                 ('MAGIC_Act','show_dialog')),
    ('MV', 'MAGIC-View',                ('MAGIC_View','show_dialog')),
