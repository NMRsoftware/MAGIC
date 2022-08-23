# MAGIC Act MAGIC View
Prepare files for automated methyl assignment 

# Instalation 
Download scripts into the NMRFAM-Sparky modleues folder:
  /Applications/nmrfam-sparky-mac/NMRFAM-SPARKY.app/Contents/Resources/python/sparky/

To create the two letter command codes (MA, MV) to open the modules in sparky edit the sparky_site.py file and add the following lines under the assignment_menu:

    ('MA', 'MAGIC-Act',                 ('MAGIC_Act','show_dialog')),
    ('MV', 'MAGIC-View',                ('MAGIC_View','show_dialog')),



