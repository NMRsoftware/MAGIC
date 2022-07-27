import Tkinter
from Tkinter import *
import types
import re
import math as m
import pyutil
import sparky
import sputil
import tkutil
import tkMessageBox
import myseq
import os
import shutil
import tkFileDialog
import string


# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
#
# Last updates: July 27, 2022
#
#
# ------------------------------------------------------------------------------
#
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
 "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L",
 "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
 "TRP": "W", "TYR": "Y", "VAL": 'V' }
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
MEA = {'I':['CD1','HD1'],'L':['CD','HD'],'V':['CG','HG'],'M':['CE','HE'], 'A':['CB','HB'], 'T':['CG2','HG2']}
#### Me_dic = Me:[#peaks in 2D, Catom-Hatom]
Me_dic = {"A":["CB-HB"], "I":["CD1-HD1"], "L":["CD1-HD1", "CD2-HD2"], "M":["CE-HE"], "T":["CG2-HG2"],"V":["CG1-HG1", "CG2-HG2"]}
#### Dictionary average (mu) methyl carbon and proton chemical shifts and standard deviation (sd) Me : [muC, sdC, muH, sdH]
MeCSdict = {"I":[13.40, 1.67, 0.68, 0.29], "L":[24.37, 1.70, 0.74, 0.28], "V":[21.395, 1.54, 0.815, 0.28], 
"M":[17.11, 1.70, 1.89, 0.40], "A":[18.97, 1.76, 1.36, 0.24], "T":[21.55, 1.11, 1.14, 0.21]}
MeTypedict = {"I":0.2, "L":0.1, "V":0.1, "M":0.2, "A":0.025, "T":0.01}

class HMQC_spectrum_menu(tkutil.option_menu):

	def __init__(self, session, parent, label, allow_no_choice = 0):

		self.session = session
		self.allow_no_choice = allow_no_choice
		self.names = self.spectrum_names()
		tkutil.option_menu.__init__(self, parent, label,
																self.names, self.default_spectrum_name())
		self.menu['postcommand'] = self.update_menu_cb

	# --------------------------------------------------------------------------
	#
	def spectrum_names(self):

		Names = []
		spectra = self.session.project.spectrum_list()
		for spectrum in spectra:
			if spectrum.dimension == 2:
				Names.append(spectrum)
		names = pyutil.attribute_values(Names, 'name')
		if self.allow_no_choice:
			names = ('',) + names
		return names
	# --------------------------------------------------------------------------
	def default_spectrum_name(self):
		return pyutil.attribute(self.session.selected_spectrum(), 'name', '')
	# --------------------------------------------------------------------------
	def update_menu_cb(self):

		current_names = self.spectrum_names()
		if current_names != self.names:
			self.names = current_names
			self.remove_all_entries()
			for name in self.names:
				self.add_entry(name)
			if not self.get() in self.names:
				self.set(self.default_spectrum_name())
	# --------------------------------------------------------------------------
	def spectrum(self):

		return sputil.name_to_spectrum(self.get(), self.session)

# --------------------------------------------------------------------------
class spectrum_menu_3D(tkutil.option_menu):

	def __init__(self, session, parent, label, allow_no_choice = 0):

		self.session = session
		self.allow_no_choice = allow_no_choice
		self.names = self.spectrum_names()
		tkutil.option_menu.__init__(self, parent, label,
																self.names, self.default_spectrum_name())
		self.menu['postcommand'] = self.update_menu_cb

	# --------------------------------------------------------------------------
	#
	def spectrum_names(self):
		Names = []
		spectra = self.session.project.spectrum_list()
		for spectrum in spectra:
			if spectrum.dimension == 3 and 'sim' not in spectrum.name:
				Names.append(spectrum)
		names = pyutil.attribute_values(Names, 'name')
		if self.allow_no_choice:
			names = ('',) + names
		return names

	# --------------------------------------------------------------------------
	#
	def default_spectrum_name(self):

		return pyutil.attribute(self.session.selected_spectrum(), 'name', '')

	# --------------------------------------------------------------------------
	#
	def update_menu_cb(self):

		current_names = self.spectrum_names()
		if current_names != self.names:
			self.names = current_names
			self.remove_all_entries()
			for name in self.names:
				self.add_entry(name)
			if not self.get() in self.names:
				self.set(self.default_spectrum_name())
	# --------------------------------------------------------------------------
	#
	def spectrum(self):

		return sputil.name_to_spectrum(self.get(), self.session)

class obj:
	pass

class methyl_dialog(tkutil.Dialog, tkutil.Stoppable):

	# ------------------------------------------------------------------------------
	#
	def __init__(self, session):

		self.session = session

		## Get the sequence for the protein, if one does not exist prompt user to load one 
		try:
			persist_path = self.session.project.save_path + '.seq'
		except:
			persist_path = ''
		sequence = []
		## If the session.project.save_path + '.seq' file exist then produce sequence = [(idx, szA)]
		if os.path.exists(persist_path):
			pszLines =[line.rstrip() for line in open(persist_path).readlines() if line.rstrip() and ">" not in line and "#" not in line]
			for line in pszLines:
				if line.split()[0] in AAA_dict.keys():
					sequence.append((int(line.split()[1]), AAA_dict[line.split()[0]]))

		if not os.path.exists(persist_path):
			tkMessageBox.showinfo('Input Error', "No Sequence file was found\n Please load a sequence using 'sq' command\nSave the project and relaunch MAGIC-Act")

			return	self.close_cb

		self.sequence = sequence
		self.selection_notice = None
		tkutil.Dialog.__init__(self, session.tk, 'MAGIC-Act')

		explain = ('Preparation for Automated Methyl Assignment:')
		w = Tkinter.Label(self.top, text = explain, justify = 'left')
		w.pack(side = 'top', anchor = 'w',padx=2)
		sl = tkutil.scrolling_list(self.top, 'Status Report', 15)
		sl.frame.pack(fill = 'both', expand = 1)
		self.summary_list = sl


	# Primary Methyl HMQC spectrum 
		self.lf2 = Tkinter.LabelFrame(self.top, text='Primary HMQC (defines w2, w3 dimensions of NOESY)')
		self.lf2.pack(fill=X, pady=5)
		self.pmframe= Tkinter.Frame(self.lf2)
		self.pmframe.pack(pady=5)
		lb = tkutil.entry_field(self.pmframe, 'Labeled Methyls: ', 'ILVMAT', 8)
		lb.frame.pack(side =LEFT, anchor='w')
		self.labeling = lb.variable
		self.sc_HMQC_2D = HMQC_spectrum_menu(session, self.pmframe, 'Spectrum: ')
		self.sc_HMQC_2D.frame.pack(side =LEFT,padx=3, expand = 1)

	# Type Labeled Methyl HMQC spectrum 
		self.lf3 = Tkinter.LabelFrame(self.top, text='Type Labeled HMQC')
		self.lf3.pack(fill=X, pady=5)
		self.tmframe= Tkinter.Frame(self.lf3)
		self.tmframe.pack(pady=5)
		tlb = tkutil.entry_field(self.tmframe, 'Selected Methyls: ', '', 8)
		tlb.frame.pack(side =LEFT, anchor='w')
		self.typelabeling = tlb.variable
		self.Type_HMQC_2D = HMQC_spectrum_menu(session, self.tmframe, 'Spectrum: ')
		self.Type_HMQC_2D.frame.pack(side =LEFT,padx=3, expand = 1)

	# Select Cm-CmHm NOESY Spectrum 
		self.sc_CCH_NOESY = spectrum_menu_3D(session, self.top, 'Cm-CmHm NOESY Spectrum: ')
		self.sc_CCH_NOESY.frame.pack(side = 'top', anchor = 'w',padx=2, pady=5)

	# Select HMBC-HMQC 3D Spectrum
		self.sc_HMBC = spectrum_menu_3D(session, self.top, 'HMBC-HMQC Spectrum: ')
		self.sc_HMBC.frame.pack(side = 'top', anchor = 'w',padx=2, pady=5)

	# PDB Information 
		self.lf1 = Tkinter.LabelFrame(self.top, text='PDB')
		self.lf1.pack(fill=X, pady=5)
		self.pdbframe= Tkinter.Frame(self.lf1)
		self.pdbframe.pack(pady=5)
		ep = tkutil.file_field2(self.pdbframe, 'PDB file', 'Browse...', file_type=[('Protein Data Bank File', '.pdb')], default_ext='.pdb')
		self.pdb_path = ep.variable
		ep.frame.pack(side=LEFT, anchor='w')
		e = tkutil.entry_field(self.pdbframe, 'distance cutoff: ', '6', 5)
		self.max_dist = e.variable
		e.frame.pack(side=LEFT,fill=X, expand=1)
		e = tkutil.entry_field(self.pdbframe, 'chain: ', 'A', 5)
		self.chians = e.variable
		e.frame.pack(side=LEFT,fill=X, expand=1)

	# Type of Methyl Labeling 
		explain = ('LV Labeling')
		w = Tkinter.Label(self.top, text = explain, justify = 'left')
		w.pack(side = 'top', anchor = 'w',padx=2)

		di = tkutil.checkbutton(self.top, 'Dimethyl', 1)
		di.button.pack(side = 'top', anchor = 'w',padx=2)
		self.dimethyl = di

		im = tkutil.checkbutton(self.top, 'Mono-methyl', 0)
		im.button.pack(side = 'top', anchor = 'w',padx=2)
		self.mono = im

		ps = tkutil.checkbutton(self.top, 'proS', 0)
		ps.button.pack(side = 'top', anchor = 'w',padx=2)
		self.proS = ps

		pr = tkutil.checkbutton(self.top, 'proR', 0)
		pr.button.pack(side = 'top', anchor = 'w',padx=2)
		self.proR = pr
	# Rename Resonances option 
		rename = tkutil.checkbutton(self.top, 'Keep Existing Assignments', 1)
		rename.button.pack(side = 'top', anchor = 'w',padx=2)
		self.rename = rename

	# PPM Tolerances for filtering noesy and finding geminal pairs
		er = tkutil.entry_row(self.top, 'PPM tolerance: ',
																		('1H', '0.01', 5),
																		('13C', '0.10', 5))
		er.frame.pack(side = 'top', anchor = 'w',padx=2)
		self.ppm_range = er

		progress_label = Tkinter.Label(self.top, anchor = 'nw')
		progress_label.pack(side = 'top', anchor = 'w',padx=2)

		br = tkutil.button_row(self.top,
													 ('Type Peaks', self.Type_peaks),
													 ('Find Geminal Pairs', self.Geminal_cb),
													  ('Check NOESY', self.NOESY_cb),
													 ('Update', self.update_cb))

		br.frame.pack(side = 'top', anchor = 'w')

		br2 = tkutil.button_row(self.top, ('Generate MAGIC Input', self.Generate_MAGIC_cb),
											('Generate FLYA Input', self.Generate_FLYA_cb),
											('Read Result', self.Read_MAGIC_results_cb),
											('Close', self.close_cb),
											('Stop', self.stop_cb))
		br2.frame.pack(side = 'top', anchor = 'w')

		self.settings = self.get_settings()
		tkutil.Stoppable.__init__(self, progress_label, br2.buttons[4])

	# ------------------------------------------------------------------------------
	#

	def get_settings(self):

		Seq ,pdb_atoms, Labels, flya_atoms =[],[], [], []
		pairs = 0
		
		for (resi,resn) in self.sequence:
			if resn in self.labeling.get():
				Seq.append(resn + str(resi))
				if self.dimethyl.state() == True or self.mono.state() == True:
					for me in Me_dic[resn]:
						Labels.append(resn + str(resi) + me)
						pdb_atoms.append(resn + str(resi) +'-' + me.split('-')[0])
						flya_atoms.append(me.split('-')[0] + '.' + resn + str(resi))
						if resn == "V" or resn == "L":
							pairs = pairs +0.5
				if self.proS.state() == True:
					Labels.append(resn + str(resi) + Me_dic[resn][0])
					pdb_atoms.append(resn + str(resi) +'-' + Me_dic[resn][0].split('-')[0])
				if self.proR.state() == True:
					Labels.append(resn + str(resi) + Me_dic[resn][-1])
					pdb_atoms.append(resn + str(resi) +'-' + Me_dic[resn][-1].split('-')[0])
		magiclabeldict = {'I':'I,CD1;', 'L':'L,CD1,CD2;', 'V':'V,CG1,CG2;', 'M':'M;', 'A':'A;', 'T':'T;'}
		magiclabeling = ''
		for me in self.labeling.get():
			magiclabeling = magiclabeling + magiclabeldict[me]
		magiclabeling = magiclabeling[:-1]
		if self.proS.state() == True:
			magiclabeling = magiclabeling.replace('L,CD1,CD2', 'L,CD1').replace('V,CG1,CG2', 'V,CG1')
		if self.proR.state() == True:
			magiclabeling = magiclabeling.replace('L,CD1,CD2', 'L,CD2').replace('V,CG1,CG2', 'V,CG2')
			
		settings = pyutil.generic_class()
		settings.magiclabeling = magiclabeling
		settings.sequence = self.sequence
		settings.hmqc_spectrum = self.sc_HMQC_2D.spectrum()
		settings.noesy_spectrum = self.sc_CCH_NOESY.spectrum()
		settings.HMBC_spec = self.sc_HMBC.spectrum()
		settings.pdb_path = self.pdb_path.get()
		settings.labeling = self.labeling.get()
		settings.typelabeling = self.typelabeling.get()
		settings.Typing_HMQC = self.Type_HMQC_2D.spectrum()
		settings.dimethyl = self.dimethyl.state()
		settings.mono = self.mono.state()
		settings.proS = self.proS.state()
		settings.proR = self.proR.state()
		settings.rename = self.rename.state()
		settings.Htol = float(self.ppm_range.variables[0].get())
		settings.Ctol = float(self.ppm_range.variables[1].get())
		settings.Labels = Labels
		settings.Seq = Seq
		settings.pdb_atoms = pdb_atoms
		settings.flya_atoms = flya_atoms
		settings.possible_pairs = int(pairs)
		return settings

	# ---------------------------------------------------------------------------
	#
	def update_cb(self):
		s = self.get_settings()
		self.stoppable_call(self.show_summary)


	# ---------------------------------------------------------------------------
	#
	def show_summary(self):
		s = self.get_settings()
		self.summary_list.clear()
		self.stoppable_loop('shifts', 100)
		sI, sL, sV, sM, sA, sT = 0, 0, 0, 0, 0, 0
		uI, uL, uV, uM, uA, uT = 0, 0, 0, 0, 0, 0
		aI, aL, aV, aM, aA, aT = 0, 0, 0, 0, 0, 0
		cI, cL, cV, cM, cA, cT = 0, 0, 0, 0, 0, 0
		for assignment in s.Labels:
			if assignment[0] == 'I':sI+=1
			if assignment[0] == 'L':sL+=1
			if assignment[0] == 'V':sV+=1
			if assignment[0] == 'M':sM+=1
			if assignment[0] == 'A':sA+=1
			if assignment[0] == 'T':sT+=1
		aTotal = 0
		pairedV = 0
		pairedL = 0
		pairedLV = 0
		Overlapped ,Used = [], []
		Isolated = []
		for peak in s.hmqc_spectrum.peak_list():
			if peak.color == 'dark red':
				Isolated.append(peak)
			if len(peak.note) != 0:
				metype = peak.note.split()[0]
				if len(metype) == 1:
					if metype == '-':metype = peak.assignment[0]
					if metype == 'I': uI= uI + 1
					if metype == 'L': 
						uL+= 1
						if 'C' in peak.note: pairedL += 0.5
					if metype == 'V': 
						uV+=1
						if 'C' in peak.note: pairedV += 0.5
					if metype == 'M': uM+=  1
					if metype == 'A': uA+=  1
					if metype == 'T': uT+=  1
				if len(metype) >= 2:
					aTotal = aTotal +1
					if 'I'in metype: aI+=1
					if 'L'in metype: 
						aL+= 1
						if 'C' in peak.note: pairedLV += 0.25
					if 'V'in metype: 
						aV+= 1
						if 'C' in peak.note: pairedLV += 0.25
					if 'M'in metype: aM+= 1
					if 'A'in metype: aA+= 1
					if 'T'in metype: aT+= 1
				if peak.assignment in s.Labels:
					if 'I'in peak.assignment[0]: cI+=1
					if 'L'in peak.assignment[0]: cL+=1
					if 'V'in peak.assignment[0]: cV+=1
					if 'M'in peak.assignment[0]: cM+=1
					if 'A'in peak.assignment[0]: cA+=1
					if 'T'in peak.assignment[0]: cT+=1
			for peak2 in s.hmqc_spectrum.peak_list():
				if peak2 != peak and peak2 not in Used:
					if ((abs(peak.frequency[0] - peak2.frequency[0]) < s.Ctol) and (abs(peak.frequency[1] - peak2.frequency[1]) < s.Htol)):
						Overlapped.extend([peak,peak2])
						Used.extend([peak,peak2])
		Reciprocated, Diagonal, Orphan, Bad= [], [], [], []
		for peak in s.noesy_spectrum.peak_list():
			if 'green' in peak.color: Reciprocated.append(peak)
			if peak.color == 'cyan':Diagonal.append(peak)
			if peak.color == 'gold':Orphan.append(peak)
			if peak.color == 'red':Bad.append(peak)

		uTotal = uI+uL+uV+uM+uA+uT
		cTotal = cI+cL+cV+cM+cA+cT
		sTotal = sI+sL+sV+sM+sA+sT

		self.summary_list.append('Number of Methyl Peaks:')
		self.summary_list.append(('Expected: %3d' %sTotal))
		self.summary_list.append(('Observed: %3d' %len(s.hmqc_spectrum.peak_list())))
		self.summary_list.append(('Overlapped: %3d' %len(Overlapped)))
		self.summary_list.append(('Assigned: %3d' %cTotal))
		self.summary_list.append(('Ambiguously Typed: %3d' %aTotal))
		self.summary_list.append('   Expected  Unambiguous Ambiguous Assigned ')
		if 'I' in s.labeling:self.summary_list.append('I    %3d         %3d       %3d       %3d' %(sI, uI, aI,cI))
		if 'L' in s.labeling:self.summary_list.append('L    %3d         %3d       %3d       %3d' %(sL, uL, aL,cL))
		if 'V' in s.labeling:self.summary_list.append('V    %3d         %3d       %3d       %3d' %(sV, uV, aV,cV))
		if 'M' in s.labeling:self.summary_list.append('M    %3d         %3d       %3d       %3d' %(sM, uM, aM,cM))
		if 'A' in s.labeling:self.summary_list.append('A    %3d         %3d       %3d       %3d' %(sA, uA, aA,cA))
		if 'T' in s.labeling:self.summary_list.append('T    %3d         %3d       %3d       %3d' %(sT, uT, aT,cT))
		if s.dimethyl == True or s.mono == True:
			self.summary_list.append('Found %d of %d possible LV pairs' %((pairedV+pairedL+pairedLV),s.possible_pairs))
			self.summary_list.append('Found %d of %d possible L pairs' %(pairedL,sL/2))
			self.summary_list.append('Found %d of %d possible V pairs' %(pairedV,sV/2))
			if pairedLV !=0: self.summary_list.append('Found %d LV pairs' %(pairedLV))
		if s.noesy_spectrum.dimension >= 3:
			self.summary_list.append('Found %d Isolated methyls' %(len(Isolated)))
			self.summary_list.append('Experimental NOEs in Cm-CmHm 3D %d' %(len(s.noesy_spectrum.peak_list())))
			self.summary_list.append('   %d Reciprocated peaks' %(len(Reciprocated)))
			if len(Diagonal) > 0:
				self.summary_list.append('   %d Diagonal peaks' %(len(Diagonal)))
			if len(Orphan) > 0:
				self.summary_list.append('   %d Orphan peaks' %(len(Orphan)))
			if len(Bad) > 0:
				self.summary_list.append('   %d Bad peaks need removal' %(len(Bad)))
		if s.pdb_path:
			NOEcount, ULN, LN, Isolated2 = self.Cacl_theo_NOE_networks()
			self.summary_list.append('')
			self.summary_list.append('PDB Statistics using a %sA cutoff' %(self.max_dist.get()))
			self.summary_list.append('   Expected # of NOEs in Cm-CmHm 3D %d' %NOEcount)
			self.summary_list.append('   %d of %d have local networks' %(LN, len(s.Labels)))
			self.summary_list.append('   %d of %d of local networks are unique' %(ULN, LN))
			self.summary_list.append('   Minimum assignment %2.1f%%' %((float(ULN)/len(s.Labels))*100.0))
			self.summary_list.append('   Maximum assignment %2.1f%%' %((float(LN)/len(s.Labels))*100.0))
			self.summary_list.append('   Found %d isolated methyls' %(len(Isolated2)))
			sx = range(0,len(Isolated2),4)[-1]
			PDBiso = '   '
			for x in range(0, len(Isolated2),4)[:-1]:
				self.summary_list.append('   %s %s %s %s' %(Isolated2[x],Isolated2[x+1],Isolated2[x+2],Isolated2[x+3]))
			lx = len(Isolated2) - sx
			for lx in range(lx):
				PDBiso = PDBiso + Isolated2[sx+lx] + ' '
			self.summary_list.append(PDBiso)

	# ---------------------------------------------------------------------------
	#
	def Cacl_theo_NOE_networks(self):
		s = self.get_settings()
		PDB_entries ={}
		Netowrks = {}
		NOEatoms = []
		tseq = ''
		if s.pdb_path:
			for line in open(s.pdb_path).readlines():
				if line[0:4] == "ATOM" or line[0:4] == 'HETA': 
					if line[21] == self.chians.get():
						if line[17:20].strip() in AAA_dict.keys():
							resid = AAA_dict[line[17:20].strip()] + line[22:26].strip() + "-" + line[12:16].strip()
							if resid in s.pdb_atoms:
								tseq = tseq + resid[0]
								PDB_entries[resid] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
								NOEatoms.append(resid)
		for atom1 in PDB_entries.keys():
			temp = []
			for atom2 in PDB_entries.keys():
				cor1 = PDB_entries[atom1]
				cor2 = PDB_entries[atom2]
				d = round(m.sqrt((cor1[0]-cor2[0])**2 + (cor1[1]-cor2[1])**2 + (cor1[2]-cor2[2])**2),1)
				if d != 0.0 and d <= float(self.max_dist.get()):
					if atom1.split('-')[0] != atom2.split('-')[0]:
						if atom2 not in temp:temp.append(atom2)
					if s.dimethyl == True and atom1.split('-')[0] == atom2.split('-')[0]:
						if atom2 not in temp:temp.append(atom2)
			if atom1 not in temp:temp.append(atom1)
			if len(temp) >= 2: Netowrks[atom1] = temp

		if s.dimethyl == True:
			for acceptor in NOEatoms:
				if acceptor in Netowrks.keys():
					donors = Netowrks[acceptor]
					for donor in donors:
						if donor[0] == 'L':
							if donor.split('-')[0] + '-CD1' not in donors and donor.split('-')[0] + '-CD1' != acceptor:
								Netowrks[acceptor].append(donor.split('-')[0] + '-CD1')
								Netowrks[donor.split('-')[0] + '-CD1'].append(acceptor)
							if donor.split('-')[0] + '-CD2' not in donors and donor.split('-')[0] + '-CD2' != acceptor:
								Netowrks[acceptor].append(donor.split('-')[0] + '-CD2')
								Netowrks[donor.split('-')[0] + '-CD2'].append(acceptor)
						if donor[0] == 'V':
							if donor.split('-')[0] + '-CG1' not in donors and donor.split('-')[0] + '-CG1' != acceptor:
								Netowrks[acceptor].append(donor.split('-')[0] + '-CG1')
								Netowrks[donor.split('-')[0] + '-CG1'].append(acceptor)
							if donor.split('-')[0] + '-CG2' not in donors and donor.split('-')[0] + '-CG2' != acceptor:
								Netowrks[acceptor].append(donor.split('-')[0] + '-CG2')
								Netowrks[donor.split('-')[0] + '-CG2'].append(acceptor)

		## Determine the number of unique local networks based on specified distance cutoff and PDB.
		## A Local network is defined by donor type count (dtc) which is the number of each methyl 
		## type present in the local network
		Isolated, allULN, ULN, fully_Unique, Used = [], [], [], [], []
		NOEcount = 0 
		LNID = {}
		for acceptor in NOEatoms:
			if acceptor not in Netowrks.keys(): Isolated.append(acceptor)
			if acceptor in Netowrks.keys():
				if s.dimethyl == True:
					if acceptor[0] == 'V' and len(Netowrks[acceptor]) -2 == 0:
						Isolated.append(acceptor)
					if acceptor[0] == 'L' and len(Netowrks[acceptor]) -2 == 0:
						Isolated.append(acceptor)
				NOEcount = NOEcount + len(Netowrks[acceptor]) -1 
				types = ''
				donors = Netowrks[acceptor]
				used = []
				for donor in donors:
					types = types+ donor[0]
				dtc = ''
				for me in s.labeling:
					if types.count(me) != 0:
						dtc = dtc + me+str(types.count(me))
				LNID[acceptor] = dtc
				if acceptor.split('-')[0] not in Used and acceptor not in Isolated:
					allULN.append(dtc)
					if s.dimethyl == True:
						Used.append(acceptor.split('-')[0])
					if dtc not in ULN: ULN.append(dtc)

		for dtc in ULN:
			nDTC = allULN.count(dtc)
			if nDTC == 1:
				fully_Unique.append(dtc)
		ULNcount = 0
		for acceptor in LNID.keys():
			if LNID[acceptor] in fully_Unique: 
				ULNcount = ULNcount +1 

		total = len(s.Labels) - tseq.count('L')/2 - tseq.count('V')/2
		return NOEcount, ULNcount , len(Netowrks), Isolated

	# ------------------------------------------------------------------------------
	#
	'''                                   Type_Peaks                                  
	
	   This section assigns amino acid type to each peak in the specified Methyl HMQC
	   and labeled methyls, storing the result in the peak.note.
	   Typing is based on probability of peaks is a given methyl PMe:
	 
	   PMe = m.exp(-((m.pow((Wc-muMEc),2)/(2*m.pow(sdMEc,2)))+(m.pow((Wh-muMEh),2)/(2*m.pow(sdMEh,2))))),5)
	
	       muMEc = average 13C chemical shift of a given methyl group 
	       sdMEc = standard deviation of the 13C chemical shift of a given methyl group 
	       myMEh = average 1H chemical shift of a given methyl group 
	       sdMEh = standard deviation of the 1H chemical shift of a given methyl group
	       Wc = observed 13C chemical shift 
	       Wh = observed 1H chemical shift
	
	   MeCSdict is a dictionary storing these values with the format Me:[muMEc. sdMec, muMeh, sdMeh]
	   If PMe/ total probability is greater than threshold specified in MeTypedict then the peaks is 
	   assigned that methyl type
	'''
	# ------------------------------------------------------------------------------
	# 
	def Assign_HMQC(self, peak, group, a1, a2):
		peak.assign(0,group, a1)
		peak.assign(1,group, a2)
		peak.show_assignment_label()
		return
	# ------------------------------------------------------------------------------
	# 
	def get_Metype(self, peak, labeling):
		total = 0
		mtype = ''
		TypeProb = {}
		for me in labeling:
			Pme = round(m.exp(-((m.pow((peak.frequency[0] -MeCSdict[me][0]),2)/(2*m.pow(MeCSdict[me][1],2)))+(m.pow((peak.frequency[1]-MeCSdict[me][2]),2)/(2*m.pow(MeCSdict[me][3],2))))),6)
			total = total + Pme
			TypeProb['P'+ me] = Pme
		total = total + 1E-07
		for me in labeling:
			if (TypeProb['P'+ me] / total) >= MeTypedict[me]:
				mtype = mtype + me
		if len(mtype) == 0:
			mtype = labeling
		return mtype

	# ------------------------------------------------------------------------------
	# 
	def Type_peaks(self):
		s = self.get_settings()
		HMQCpeaks = sorted(s.hmqc_spectrum.peak_list(), key = lambda x: (x.frequency[0]), reverse=True)
		#### Assign potential methyl type for each peak in HMQC based on labeling 
		if len(s.typelabeling) == 0:
			if len(self.labeling.get()) == 0:
				tkinter.messagebox.showinfo('Input Error', "No Labeling Specified \n Specify Labeled methyls and try again")
				return
			x = 0
			comment = ''
			for peak in HMQCpeaks:
				x = x+1
				if len(peak.note) == 0:
					if peak.is_assigned == 1 and peak.assignment in s.Labels:
						mtype = peak.assignment[0]
					else:
						mtype = self.get_Metype(peak,s.labeling)
					peak.note = mtype
				## If the peak has no assignment give it assingment with group.symbol = metype and group.number = list index, with atom1 = C and atom2 = H
				if peak.is_assigned == 0:
					group = peak.note.split()[0].lower() + str(x)
					self.Assign_HMQC(peak, peak.note.split()[0].lower() + str(x), 'C','H')
				## If the methyl type in the note does not match the group.symbol ubdate it 
				if len(peak.note) >= 1 and peak.assignment not in s.Labels:
					if peak.note.split()[0].lower() != peak.resonances()[0].group.symbol:
						self.Assign_HMQC(peak, peak.note.split()[0].lower() + str(x), 'C','H')
						peak.label.color = 'white'
				if not re.search('[L,V]', peak.note.split()[0]):
					peak.label.color = 'white'

		if len(s.typelabeling) != 0:
			if s.Typing_HMQC.name == s.hmqc_spectrum.name:
				tkinter.messagebox.showinfo('Input Error', "Please select a type labeled Me HMQC")
				return
			if s.Typing_HMQC.name != s.hmqc_spectrum.name:
				retyped = []
				used = []
				for peak in s.Typing_HMQC.peak_list():
					if peak.note.split():
						if len(peak.note.split()[0]) ==1 and peak.note.split()[0] in s.typelabeling:
							mtype = peak.note.split()[0]
					if len(peak.note) == 1 and peak.note.split()[0] in s.typelabeling:
						mtype = peak.note
					if peak.is_assigned == 1 and peak.assignment in s.Labels:
						mtype = peak.assignment[0]
					if len(peak.note) == 0:
						mtype = self.get_Metype(peak,s.typelabeling)
						peak.note = mtype
					Matches = [hmqc for hmqc in HMQCpeaks if (abs(hmqc.frequency[0] -peak.frequency[0]) < s.Ctol) and (abs(hmqc.frequency[1] -peak.frequency[1]) < s.Htol) and hmqc not in retyped]
					diff = []
					if len(Matches) >= 1:
						for hmqc in Matches: 
							diff.append(abs(peak.frequency[0] - hmqc.frequency[0]) + abs(peak.frequency[1] - hmqc.frequency[1]))
						peak2 = Matches[diff.index(min(diff))]
						if peak2.assignment in s.Labels:
							self.Assign_HMQC(peak, peak2.resonances()[0].group.name, peak2.resonances()[0].atom.name, peak2.resonances()[1].atom.name)
						if peak2.assignment not in s.Labels:
							self.Assign_HMQC(peak, mtype.lower() + str(idx), 'C','H')
						retyped.append(peak2)
						if peak2.note.split(): peak2type = peak2.note.split()[0]
						else: peak2type = peak2.note
						peak2.note = peak2.note.replace(peak2type, mtype)
				for peak3 in s.hmqc_spectrum.peak_list():
					if peak3 not in retyped:
						group = peak3.resonances()[0].group.name
						for me in s.typelabeling:
							peak3.note = peak3.note.replace(me,'')
							group = group.replace(me.lower(),'')
						self.Assign_HMQC(peak3, group, peak3.resonances()[0].atom.name,peak3.resonances()[1].atom.name)

		self.stoppable_call(self.show_summary)
	# ------------------------------------------------------------------------------
	# 
	def Geminal_cb(self):
		s = self.get_settings()
		self.stoppable_call(self.Geminal_Pairs_cb)
		message = ('Checking possible combination %d\nFound %d of possible geminal pairs %d' % (self.count, self.LVcount, s.possible_pairs))
		self.progress_report(message)

	# ------------------------------------------------------------------------------
	#
	#                               Geminal_Pairs_cb                                
	#
	#   This section uses the HMBC-HMQC 3D and the typed methyl HMQC 2D to find geminal 
	#   pairs
	#
	#   Only HMQC peaks with L or V in their type are considered 
	#   
	#   To find the geminal pair (Me1, Me2) search the HMBC3D to find the acceptor cross peak HMBC_A. 
	#   A cross peaks with (abs(w2 - Me1c) < Ctol) and (abs(w3 -Me1h) < Htol) 
	#   In the event that more than one cross peak matches this condition, common in overlapped spectra,
	#   the cross peaks with the smallest deviation from Me1 is selected. 
	#   The w1 of HMBC_A is then used to find the donor methyl (Me2) and cross peak in the HMBC (HMBC_D). 
	#   A short list of possible candidates for Me2 is found by filtering the HMQC peak list according to 
	#   abs(HMBC_A_w1 - Me2c) < Ctol. 
	#   Then the correct methyl is found by further searching the HMBC 3D to find peaks with 
	#       (abs(w1 - Me1c) < Ctol) & (abs(w2 -Me2c) < Ctol) & (abs(w3-Me2h) < Htol)
	#   generating HMBC_D
	#   In the event that more than one HMBC cross peaks satisfies this condition the one with the smallest deviation 
	#   from  w1 - Me1c and w2 - Me2c is selected. 
	#   Once a single Me1, Me2, HMBC_A, are HBC_D are found they are assigned and paired in the notes section, 
	#   and the peaks are added to list to prevent attempted reuse.
	#   If the peak is not already assigned it is give then group.symbol u and group.number = LVcount


	def Geminal_Pairs_cb(self):
		s = self.get_settings()
		if s.hmqc_spectrum.dimension != 2:
			tkMessageBox.showinfo('Input Error', "Please Select 2D Methyl Spectrum")
			return
		if s.HMBC_spec.dimension != 3:
			tkMessageBox.showinfo('Input Error', "Please Select 3D NOESY Spectrum")
			return
		self.Type_peaks()	# run type peaks so that all the peaks will have a type and temporary assignment
		HMQCpeaks = sorted(s.hmqc_spectrum.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
		LVpeaks = []
		for me in range(len(HMQCpeaks)):
			if 'L' in HMQCpeaks[me].note or 'V' in HMQCpeaks[me].note :
				if "C" in HMQCpeaks[me].note:
					HMQCpeaks[me].note = HMQCpeaks[me].note.split()[0]
				LVpeaks.append(me)
		Used = []
		Paired = []
		idxlist = []
		self.LVcount = 0
		HMBCpeaks = s.HMBC_spec.peak_list()
		for peak in HMQCpeaks:
			if (abs(peak.frequency[0] - peak.frequency[1]) < s.Ctol):
				peak.note = 'Diagonal'
				peak.color = 'cyan'
		self.count = 0
		for x in LVpeaks:
			if x not in Paired: 
				me1 =  HMQCpeaks[x]
				HMBC_A = [peak for peak in range(len(HMBCpeaks)) if (abs(HMBCpeaks[peak].frequency[1] - me1.frequency[0]) < s.Ctol) and (abs(HMBCpeaks[peak].frequency[2] - me1.frequency[1]) < s.Htol) and peak not in Used]
				if len(HMBC_A) >=1:
					diff = []
					for idx in HMBC_A:
						diff.append(abs(HMBCpeaks[idx].frequency[1] - me1.frequency[0]) + abs(HMBCpeaks[idx].frequency[2] - me1.frequency[1]))
					HMBC_A = HMBC_A[diff.index(min(diff))]
					# now find the possible peaks in the HMQC that have the same w1 frequency as this single HMBC_A peak
					Me2list = [me for me in LVpeaks if (abs(HMBCpeaks[HMBC_A].frequency[0] - HMQCpeaks[me].frequency[0]) < s.Ctol) and me not in Paired]
					HMBC_D,diff2,Me2list_2 = [], [], []
					for y in Me2list:
						me2 = HMQCpeaks[y]
						self.count = self.count + 1
						message = ('Checking possible combination %d \nFound %d of possible geminal pairs %d' % (self.count, self.LVcount, s.possible_pairs))
						self.progress_report(message)
						if x not in Paired and y != x:
							HMBC_D_short = [peak2 for peak2 in range(len(HMBCpeaks)) if (abs(HMBCpeaks[peak2].frequency[0] - me1.frequency[0]) < s.Ctol) and (abs(HMBCpeaks[peak2].frequency[1] - me2.frequency[0]) < s.Ctol) and (abs(HMBCpeaks[peak2].frequency[2] - me2.frequency[1]) < s.Htol) and peak2 not in Used]
							if len(HMBC_D_short) >= 1:
								for idy in HMBC_D_short:
									d = abs(HMBCpeaks[idy].frequency[0] - me1.frequency[0]) + abs(HMBCpeaks[idy].frequency[1] - me2.frequency[0])
									diff2.append(d)
									Me2list_2.append(y)
									HMBC_D.append(idy)
					if len(Me2list_2) >= 1:
						me2 = HMQCpeaks[Me2list_2[diff2.index(min(diff2))]]
						HMBC_D = HMBC_D[diff2.index(min(diff2))]
						Paired.append(x)
						Paired.append(Me2list_2[diff2.index(min(diff2))])
						Used.append(HMBC_A)
						Used.append(HMBC_D)
						self.LVcount = self.LVcount + 1 
						if me1.resonances()[0].group.number not in idxlist:
							idxlist.append(me1.resonances()[0].group.number)
							index = me1.resonances()[0].group.number
						if me1.resonances()[0].group.number in idxlist:
							index = me2.resonances()[0].group.number 
							idxlist.append(me2.resonances()[0].group.number)
						self.add_note_tempassignment(HMBCpeaks[HMBC_A], HMBCpeaks[HMBC_D], me1, me2, index)
						message = ('Checking possible combination %d \nFound %d of possible geminal pairs %d' % (self.count, self.LVcount, s.possible_pairs))
						self.progress_report(message)

		## Clean up LV assignment, only LV peaks which are paired should have C1-H1/C2-H2 as atom names 
		HMQCpeaks = sorted(s.hmqc_spectrum.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
		x = 0
		for peak in HMQCpeaks:
			x= x+1 
			if 'L' in peak.note or 'V' in peak.note:
				if "C" not in peak.note:
					peak.label.color = 'purple'
					if peak.assignment not in s.Labels:
						group = peak.note.split()[0].lower() + str(x)
						peak.assign(0,group,'C')
						peak.assign(1,group,'H')
						peak.show_assignment_label()

		self.stoppable_call(self.show_summary)

	# ------------------------------------------------------------------------------
	# 

	def add_note_tempassignment(self, hmbc1, hmbc2, me1, me2, count):
		MEA2 = {'L':['CD','HD'],'V':['CG','HG'],'LV':['C','H']}
		s = self.get_settings()
		group1 = me1.resonances()[0].group.name
		group2 = me2.resonances()[0].group.name
		## Option 1 either me1 or me2 has a residue specific assignment
		if me1.assignment in s.Labels or me2.assignment in s.Labels:
			if group1 == group2: # if me1 or me2 have the same sequential assignment, then set the methyl type to the assigned methyl type
				group = group1
				self.Assign_HMBC(hmbc1,group, me2.resonances()[0].atom.name, me1.resonances()[0].atom.name, me1.resonances()[1].atom.name)
				self.Assign_HMBC(hmbc2,group, me1.resonances()[0].atom.name, me2.resonances()[0].atom.name, me2.resonances()[1].atom.name)
				me1.note = me1.resonances()[0].group.symbol
				me2.note = me1.resonances()[0].group.symbol
				metype = group1[0]
			if group1 != group2: # if me1 or me2 do not have a residue specific assignment or different residue specific assignments
				if me1.resonances()[0].group.symbol.upper() != me2.resonances()[0].group.symbol.upper():
					metype = 'LV'
					group = 'LV' + str(count)
				if me1.resonances()[0].group.symbol.upper() == me2.resonances()[0].group.symbol.upper():
					metype = group1[0]
					group =metype.lower()+str(count)
				self.Assign_HMBC(hmbc1,group, 'C2', 'C1', 'H1')
				self.Assign_HMBC(hmbc2,group, 'C1', 'C2', 'H2')
				self.Assign_HMQC(me1, group, 'C1','H1')
				self.Assign_HMQC(me2, group, 'C2','H2')
				if me1.assignment in s.Labels:
					me1.note = metype + ':' + me1.assignment
				if me1.assignment not in s.Labels:
					me1.note = metype
				if me2.assignment in s.Labels:
					me2.note = metype + ':' + me2.assignment
				if me2.assignment not in s.Labels:
					me2.note = metype
		if me1.assignment not in s.Labels and me2.assignment not in s.Labels:
			me1_type = me1.note.replace('T','').replace('A','').split(' ')[0]
			me2_type = me2.note.replace('T','').replace('A','').split(' ')[0]
			if len(me1_type) < len(me2_type):
				metype = me1_type
			if len(me2_type) < len(me1_type):
				metype = me2_type
			if len(me1_type) == len(me2_type) and me1_type == me2_type:
				metype = me1_type
			if len(me1_type) == len(me2_type) and me1_type != me2_type:
				metype = 'LV'
			me1.note = metype
			me2.note = metype
			group =metype.lower()+str(count)
			self.Assign_HMBC(hmbc1,group, 'C2', 'C1', 'H1')
			self.Assign_HMBC(hmbc2,group, 'C1', 'C2', 'H2')
			self.Assign_HMQC(me1, group, 'C1','H1')
			self.Assign_HMQC(me2, group, 'C2','H2')
		me1.note = me1.note + ' ;'+ me2.assignment
		me1.label.color = 'white'
		me2.note = me2.note + ' ;'+ me1.assignment
		me2.label.color = 'white'
		return

	# ------------------------------------------------------------------------------
	# 
	def Assign_HMBC(self, peak, group, a1, a2, a3):
		peak.assign(0,group, a1)
		peak.assign(1,group, a2)
		peak.assign(2,group, a3)
		peak.show_assignment_label()
		return
	# ------------------------------------------------------------------------------
	# 
	def Assign_3D_NOESY(self, peak, group1, group2, a1, a2, a3):
		peak.assign(0,group1, a1)
		peak.assign(1,group2, a2)
		peak.assign(2,group2, a3)
		peak.show_assignment_label()
		return
	# ------------------------------------------------------------------------------
	# 
	# Find all possible donors for a given methyl in HMQC with NOESY w2 = me1C, and NOEYS w3 = me1H
	#
	def Find_3D_Donors(self, PeaksList, me1C, me1H, tList):

		s = self.get_settings()
		Donors = []
		for peak in range(len(PeaksList)):
			if (abs(PeaksList[peak].frequency[1] - me1C) < s.Ctol) and (abs(PeaksList[peak].frequency[2] - me1H) < s.Htol) and peak not in tList:
				Donors.append(peak)
		return Donors
	# ------------------------------------------------------------------------------
	# 
	# Find all possible reciprocated acceptor cross peak(s) serving as donor in original query strip. NOESY w1 = me1C, w2 = me2C, and NOEYS w3 = me2H
	#
	def Find_3D_Acceptor(self, PeaksList,me1C, me2C, me2H, tList):
		s = self.get_settings()
		Acceptors = []
		for peak in range(len(PeaksList)):
			if (abs(PeaksList[peak].frequency[0] - me1C) < s.Ctol) and (abs(PeaksList[peak].frequency[1] - me2C) < s.Ctol) and (abs(PeaksList[peak].frequency[2] - me2H) < s.Htol) and peak not in tList:
				Acceptors.append(peak)
		return Acceptors
	# ------------------------------------------------------------------------------
	# 
	#                               Check_NOESY_cb                                
	#
	#   This section uses the Cm-CmHm 3D NOESY and the types methyl HMQC 2D to: 
	#   	1) If dimethyl labeled sample check for symmetry between geminal pairs
	#		2) Check for reciprocity of donor acceptor pairs
	#
	#   1) Symmetry between geminal pairs
	#		Only paired L/V peaks in the 2D HMQC are considered, the symmetry of the 
	#		two NOESY strips is determined by the ratio of the maximum number of peaks observed
	#		in both stips max(len(Me1Donors), len(Me2Doonrs)) and the weighted sum of the peaks 
	#		observed in each individual strip (MeXscore). Peaks in both NOESY strips  that satisfy
	#		abs(Me1.NOESY.w1 - Me2.NOESY.w1) < Ctol, are considered matched and assigned a value of 1.
	#		The impact on unmatched peaks (p) is determined by p = meXmin/peak.data_height, where meXmin
	#		is the minimum data height observed in meX NOESY strip.
	#		If MeXscore is > 0.85 the geminal pairing is correct and the HMQC peak label is colored green. 
	#		If MeXscore is <0.85 but > 0.6 the geminal pairing is most likely symmetric but require further
	#		curation and the HMQC peak label color is set to gold. 
	#		If MeXscor is < 0.6 the strips are not symmetric and the geminal pairing may need reevaluation, the 
	#		HMQC peak label color is set to magenta. 
	#	2) Reciprocity of donor acceptor pairs 
	#		


	def color_geminal(self, peak, score):

		if score >= 0.85: 
			peak.label.color = 'green'
			peak.note = peak.note.split(':')[0] + ' : %0.2f' %score
		elif score < 0.85 and score >= 0.6: 
			peak.label.color = 'gold'
			peak.note = peak.note.split(':')[0] + ' : %0.2f' %score
		elif score < 0.60:
			peak.label.color = 'magenta'
			peak.note = peak.note.split(':')[0] + ' : %0.2f' %score
		return


	def NOESY_cb(self):
		s = self.get_settings()
		self.count = 0
		if s.noesy_spectrum.dimension != 3:
			tkMessageBox.showinfo('Input Error', "Please Select 3D NOESY Spectrum")
			return
		NOESYpeaks = s.noesy_spectrum.peak_list()
		Used = []
		Pairs, Pairs2, label2idx = {},{}, {}
		GPchecked = []
		HMQCpeaks = sorted(s.hmqc_spectrum.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
		## Add Geminal Pair Cross Peaks first 
		if s.dimethyl == True:
			LVpeaks, Overlapped,checked, paired = [], [], [], []
			for me in HMQCpeaks:
				if 'L' in me.note or 'V' in me.note:
					if 'C' in me.note:
						LVpeaks.append(me)
			paired = []
			## Make a dictionary to relate the geminal pairs
			for x in range(0,len(LVpeaks)-1,1):
				peakOverlap = []
				label2idx[LVpeaks[x].assignment] = x
				for y in range(1,len(LVpeaks),1):
					label2idx[LVpeaks[y].assignment] = y
					if x !=y:
						if x not in paired and y not in paired: 
							if LVpeaks[x].resonances()[0].group.name == LVpeaks[y].resonances()[0].group.name:
								Pairs[x] = y
								Pairs2[x] = y
								Pairs2[y] = x
								paired.extend([x,y])
						if ((abs(LVpeaks[x].frequency[0] - LVpeaks[y].frequency[0]) < s.Ctol) and (abs(LVpeaks[x].frequency[1] - LVpeaks[y].frequency[1]) < s.Htol)):
							if LVpeaks[x].assignment not in peakOverlap and LVpeaks[x].assignment not in checked: 
								peakOverlap.append(LVpeaks[x].assignment)
								checked.append(LVpeaks[x].assignment)
							if LVpeaks[y].assignment not in checked: 
								peakOverlap.append(LVpeaks[y].assignment)
								checked.append(LVpeaks[y].assignment)
				if len(peakOverlap) > 0:
					Overlapped.append(peakOverlap)
### Check the overlapped peaks.
			for olpeaks in Overlapped:
				olDonors, olMatched, resolved = [], [], []
				for peak in olpeaks:
					me2 = LVpeaks[label2idx[peak]] ## overlapped peak 
					for peak in self.Find_3D_Donors(NOESYpeaks,me2.frequency[0],me2.frequency[1],[]):
						if peak not in olDonors:
							olDonors.append(peak)
				# print(olDonors)
				olmin = min([NOESYpeaks[peak].data_height for peak in olDonors])
				multimatch = 'matches '
				for peak in olpeaks:
					me1Matched, me2Matched, me1score = [],[],[]
					me1 = LVpeaks[Pairs2[label2idx[peak]]] ## resolved peak 
					me2 = LVpeaks[label2idx[peak]] ## overlapped peak 
					me1Donors = self.Find_3D_Donors(NOESYpeaks,me1.frequency[0],me1.frequency[1],[])
					if len(me1Donors) == 0 or len(olDonors) == 0:
						me1.label.color = 'red'
						me2.label.color = 'red'
					if len(me1Donors) >=1 and len(olDonors) >=1:
						me1min = min([NOESYpeaks[peak].data_height for peak in me1Donors])
						for peak in me1Donors: 
							if (abs(NOESYpeaks[peak].frequency[0] - me2.frequency[0]) < s.Ctol) and peak not in me1Matched: ## Find geminal cross peak in me1 trajectory 
								me1Matched.append(peak)
								me1score.append(1.1)
								NOESYpeaks[peak].color = 'dark green'
								Used.append(peak)
								self.Assign_3D_NOESY(NOESYpeaks[peak],me2.resonances()[0].group.name, me1.resonances()[0].group.name,me2.resonances()[0].atom.name,me1.resonances()[0].atom.name,me1.resonances()[1].atom.name)
							for peak2 in olDonors:
								if (abs(NOESYpeaks[peak2].frequency[0] - me1.frequency[0]) < s.Ctol) and peak2 not in me2Matched: ## Find geminal cross peak in me2 trajectory
									NOESYpeaks[peak2].color = 'dark green'
									me2Matched.append(peak2)
									olMatched.append(peak2)
									Used.append(peak2)
									self.Assign_3D_NOESY(NOESYpeaks[peak2],me1.resonances()[0].group.name, me2.resonances()[0].group.name,me1.resonances()[0].atom.name,me2.resonances()[0].atom.name,me2.resonances()[1].atom.name)
								if peak not in me1Matched:
									if peak2 not in me2Matched: 
										if (abs(NOESYpeaks[peak].frequency[0] - NOESYpeaks[peak2].frequency[0]) < s.Ctol):  ## Find shared cross peaks in me1 and me2 trajectories 
											me1score.append(1.1)
											me1Matched.append(peak)
											me2Matched.append(peak2)
											olMatched.append(peak2)
						me1Unmatched = [peak for peak in me1Donors if peak not in me1Matched]
						for peak in me1Unmatched:
							me1score.append(me1min/NOESYpeaks[peak].data_height)
						resolved.append(me1score)
				olUnmatched = [peak for peak in olDonors if peak not in olMatched]
				olpenalty = sum([olmin/NOESYpeaks[peak].data_height for peak in olUnmatched])
				# print(olpenalty)
				for peak, score in zip(olpeaks,resolved):
					me1 = LVpeaks[Pairs2[label2idx[peak]]] ## resolved peak 
					me2 = LVpeaks[label2idx[peak]] ## overlapped peak
					nDonors = len(score)
					me1score = (sum(score) - score.count(1.1)*0.1)/nDonors
					me2score = (score.count(1.1) + olpenalty)/nDonors
					if me2score > 1.0: me2score = 1.00
					self.color_geminal(me1, me1score)
					self.color_geminal(me2, me2score)
					# print('checking symmetry of %s and %s noesy strips ' %(me1.assignment, me2.assignment))
					# print('%s Matched %d of %d possible matches with score of %0.2f' %(me1.assignment, score.count(1.1), len(me1Donors),me1score))
					# print('%s Matched %d of %d possible matches wiht score of %0.2f' %(me2.assignment, score.count(1.1), len(olDonors), me2score))
					# print()
			Overlapped2 = [peak for subset in Overlapped for peak in subset]
			i = 0
			for x in Pairs.keys():
				me1 = LVpeaks[x]
				me2 = LVpeaks[Pairs[x]]
				if me1.assignment not in Overlapped2 and me2.assignment not in Overlapped2:
					i= i+1
					#print('checking LV pair %d of %d ' %(i, len(LVpeaks)/2))
					#print('checking symmetry of %s and %s noesy strips ' %(me1.assignment, me2.assignment))
					me1Donors = self.Find_3D_Donors(NOESYpeaks,me1.frequency[0],me1.frequency[1],[])
					me2Donors = self.Find_3D_Donors(NOESYpeaks,me2.frequency[0],me2.frequency[1],[])
					me1Matched, me2Matched = [],[]
					if len(me1Donors) == 0 or len(me2Donors) == 0:
						me1.label.color = 'red'
						me2.label.color = 'red'
					if len(me1Donors) >=1 and len(me2Donors) >=1:
						me1min = min([NOESYpeaks[peak].data_height for peak in me1Donors])
						me2min = min([NOESYpeaks[peak].data_height for peak in me2Donors])
						for peak in me1Donors: 
							if (abs(NOESYpeaks[peak].frequency[0] - me2.frequency[0]) < s.Ctol) and peak not in me1Matched: ## Find geminal cross peak in me1 trajectory 
								me1Matched.append(peak)
								NOESYpeaks[peak].color = 'dark green'
								Used.append(peak)
								self.Assign_3D_NOESY(NOESYpeaks[peak],me2.resonances()[0].group.name, me1.resonances()[0].group.name,me2.resonances()[0].atom.name,me1.resonances()[0].atom.name,me1.resonances()[1].atom.name)
							for peak2 in me2Donors:
								if (abs(NOESYpeaks[peak2].frequency[0] - me1.frequency[0]) < s.Ctol) and peak2 not in me2Matched:  ## Find geminal cross peak in me2 trajectory
									me2Matched.append(peak2)
									NOESYpeaks[peak2].color = 'dark green'
									Used.append(peak2)
									self.Assign_3D_NOESY(NOESYpeaks[peak2],me1.resonances()[0].group.name, me2.resonances()[0].group.name,me1.resonances()[0].atom.name,me2.resonances()[0].atom.name,me2.resonances()[1].atom.name)
								if peak not in me1Matched:
									if peak2 not in me2Matched: 
										if (abs(NOESYpeaks[peak].frequency[0] - NOESYpeaks[peak2].frequency[0]) < s.Ctol):  ## Find shared cross peaks in me1 and me2 trajectories 
											me1Matched.append(peak)
											me2Matched.append(peak2)
						me1Unmatched = [peak for peak in me1Donors if peak not in me1Matched]
						me2Unmatched = [peak for peak in me2Donors if peak not in me2Matched]
						me1score = len(me1Matched)
						me2score = len(me2Matched)
						nDonors = max([len(me1Donors), len(me2Donors)])
						for peak in me1Unmatched:
							p = me1min/NOESYpeaks[peak].data_height
							#print('penalty = %3.2f for peak with SNR of %4.0f' %(p, NOESYpeaks[peak].data_height/s.noesy_spectrum.noise))
							me1score = me1score + p
						for peak in me2Unmatched:
							p = me2min/NOESYpeaks[peak].data_height
							#print('penalty = %3.2f for peak with SNR of %4.0f' %(p, NOESYpeaks[peak].data_height/s.noesy_spectrum.noise))
							me2score = me2score + p
						Me1score = me1score/float(nDonors)
						Me2score = me2score/float(nDonors)
						# print('%d of %d donors observed for %s are matched; score %0.2f ' %(len(me1Matched), nDonors, me1.assignment, Me1score))
						# print('%d of %d donors observed for %s are matched; score %0.2f \n ' %(len(me2Matched), nDonors, me2.assignment, Me2score))
						self.color_geminal(me1, Me1score)
						self.color_geminal(me2, Me2score)

		## NOESY strip terminology: the w2,w3 frequencies in the NOEYS represent the acceptor group, while w1 represents the donor group
		## For each methyl in the HMQC 2D extract the NOESY strip, each peak in the strip represents a possible donor
		for x in range(len(HMQCpeaks)):
			me1 =  HMQCpeaks[x]
			Donor_check =self.Find_3D_Donors(NOESYpeaks,me1.frequency[0],me1.frequency[1],[])
			Donors = self.Find_3D_Donors(NOESYpeaks,me1.frequency[0],me1.frequency[1],Used)
			if 'C' not in me1.note and len(Donor_check) == 0:
				me1.color = 'dark red'
			if 'C' in me1.note and len(Donor_check) <= 1:
				me1.color = 'dark red'
			if len(Donors) >=1:
				me1.color = 'white'
				for donor in Donors:
					## Find the possible peaks in the HMQC that have the same w1 frequency as the donor peak in the NOESY
					Me2list = [me for me in range(len(HMQCpeaks)) if (abs(NOESYpeaks[donor].frequency[0] - HMQCpeaks[me].frequency[0]) < s.Ctol)]
					NOESY_A, diff2, Me2list_2= [],[],[]
					### For each possible w1 match check the NOESY strip for a cross peak with w1 = acceptor 13C freq 
					for y in Me2list:
						self.count = self.count + 1
						message = ('Checking possible combination %d in NOESY' % (self.count))
						self.progress_report(message)
						if y != x:
							## Find all the cross peaks in the NOESY that have w1 = me1C and w2-w3 =  me2C-H values
							me2 = HMQCpeaks[y]
							NOESY_A_short = self.Find_3D_Acceptor(NOESYpeaks,me1.frequency[0], me2.frequency[0], me2.frequency[1], Used)
							if len(NOESY_A_short) >= 1:
								for idy in NOESY_A_short:
									d = abs(NOESYpeaks[idy].frequency[0] - me1.frequency[0]) + abs(NOESYpeaks[donor].frequency[0] - NOESYpeaks[idy].frequency[1]) + abs(HMQCpeaks[y].frequency[0]-NOESYpeaks[donor].frequency[0])
									diff2.append(d)
									Me2list_2.append(y)
									NOESY_A.append(idy)
					if len(Me2list_2) >= 1:
						me2 = HMQCpeaks[Me2list_2[diff2.index(min(diff2))]]
						NOESY_A = NOESY_A[diff2.index(min(diff2))]
						## Assign the 
						self.Assign_3D_NOESY(NOESYpeaks[donor],me2.resonances()[0].group.name, me1.resonances()[0].group.name,me2.resonances()[0].atom.name,me1.resonances()[0].atom.name,me1.resonances()[1].atom.name)
						self.Assign_3D_NOESY(NOESYpeaks[NOESY_A],me1.resonances()[0].group.name, me2.resonances()[0].group.name,me1.resonances()[0].atom.name,me2.resonances()[0].atom.name,me2.resonances()[1].atom.name)
						NOESYpeaks[NOESY_A].color = 'dark green'
						NOESYpeaks[NOESY_A].note = ''
						NOESYpeaks[donor].color = 'dark green'
						Used.append(NOESY_A)
						Used.append(donor)
					if len(Me2list_2) >= 2:
						dnote = ''
						for p,d  in zip(Me2list_2, diff2):
							dnote = dnote + HMQCpeaks[p].assignment.split('-')[0] +' %4.3f; ' %(d)
						NOESYpeaks[donor].note = dnote[:-1]
						NOESYpeaks[donor].color = 'light green'

		message = ('Found %d Reciprocated peaks in NOESY' % (len(Used)))
		self.progress_report(message)
		Unused = []
		for peak in range(len(NOESYpeaks)):
			if peak not in Used: Unused.append(NOESYpeaks[peak])
		Diagonal = []
		Orphans = []
		for peak in Unused:
			if (abs(peak.frequency[0]-peak.frequency[1]) < s.Ctol):
				Diagonal.append(peak)
				peak.color = 'cyan'
				peak.note = 'Diagonal'
			for me in HMQCpeaks:
				if (abs(peak.frequency[0]-me.frequency[0]) < s.Ctol):
					for me2 in HMQCpeaks:
						if me2 != me:
							if (abs(peak.frequency[1]-me2.frequency[0]) < s.Ctol) and (abs(peak.frequency[2]-me2.frequency[1]) < s.Htol) and peak not in Diagonal:
								Orphans.append(peak)
								peak.color = 'gold'
								peak.note = 'Orphan'
								if peak.is_assigned == 1:
									peak.assign(0, '?', '?')
									peak.assign(1, '?', '?')
									peak.assign(2, '?', '?')
									try: 
										peak.label.shows_assignment = 0
									except:
										continue
		for peakid in Used:
			peak = NOESYpeaks[peakid]
			if (abs(peak.frequency[0]-me.frequency[1]) < s.Ctol):
				Diagonal.append(peak)
				peak.color = 'cyan'
				peak.note = 'Diagonal'

		for peak in Unused:
			if peak not in Orphans and peak not in Diagonal:
				peak.color = 'red'
				peak.note = 'Bad frequency'
		self.stoppable_call(self.show_summary)
		return 

	# ------------------------------------------------------------------------------
	# 
	def Generate_MAGIC_cb(self):

		s = self.get_settings()
		outdir = self.session.project.sparky_directory + '/' + self.session.project.save_path.split('/')[-2] + '/MAGIC/'
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		HMQCpeaks = sorted(s.hmqc_spectrum.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
		replacements = {'?-?':''}
		for x in range(len(HMQCpeaks)):
			peak = HMQCpeaks[x]
			if '?' not in peak.assignment:
				if s.rename == False:
					replacements[peak.assignment]=str(x+1)
				if s.rename == True:
					if peak.assignment not in s.Labels:
						replacements[peak.assignment]=str(x+1)
		
		HMCQnew = open(outdir + s.hmqc_spectrum.name + '_magic.list', 'w')
		for x in range(len(HMQCpeaks)):
			peak = HMQCpeaks[x]
			passign = peak.assignment
			note = peak.note
			if len(note) == 0:
				tkMessageBox.showinfo('Input Error', "Missing peak type")
				return
			if s.rename == True:
				if passign in s.Labels:
					if len(note.split()) == 2:
						note = '- ' + note.split()[1]
					if len(note.split()) == 1:
						note = '-'
			if passign not in s.Labels:
				passign = str(x+1)
			for src, target in replacements.iteritems():
				note = note.replace(src, target)
				passign = passign.replace(src, target)
			HMCQnew.write("  %-14s %10.3f %10.3f\t\t%s\n" %(passign, peak.frequency[0],peak.frequency[1],note))
		HMCQnew.close()

		CCHnoesy = open(outdir + s.noesy_spectrum.name + '.list','w')
		CCHnoesy.write('13C;13C;1H\n0.1;0.1;0.01\n')
		for noe in s.noesy_spectrum.peak_list():
			CCHnoesy.write("%17s %10.3f %10.3f %10.3f  %11.0f    \n" % ('?-?-?',noe.frequency[0], noe.frequency[1], noe.frequency[2], noe.data_height))
		CCHnoesy.close()

		seqauto = open(outdir + 'seq.auto','w')
		for res in s.Seq:
			seqauto.write(res+'\n')
		seqauto.close()

		Magic_start = open(outdir + 'start.txt','w')
		Magic_start.write('############ Input ###########\n')
		Magic_start.write('##2d reference spectrum:\n%s\n' % (s.hmqc_spectrum.name + '_magic.list'))
		Magic_start.write('##3d NOESY:\n%s\n' % (s.noesy_spectrum.name + '.list'))
		Magic_start.write('##pdb:\n%s\n' % (s.pdb_path.split('/')[-1]))
		Magic_start.write('##SEQ:\nseq.auto\n')
		Magic_start.write('##labeling:\n%s\n' % (s.magiclabeling))
		Magic_start.write('############ parameters ###########\n##stereospecific=2;double methyl=1;mono-methyl=0:\n')
		if s.dimethyl == True: Magic_start.write('1\n')
		if s.mono == True: Magic_start.write('0\n')
		if s.proS == True: Magic_start.write('2\n')
		if s.proR == True: Magic_start.write('2\n')
		Magic_start.write('##score threshold factor:\n1\n##Distance threshold:\n7 10\n##Score tolerance for one-by-one step (off or on):\noff\n##Area of expected conformational changes:\n')
		Magic_start.close()
		try:
			shutil.copy(s.pdb_path,outdir)
		except OSError:
			pass
		tkMessageBox.showinfo('Finished', "Input files for MAGIC\nSuccessfully generated\n%s" %(outdir))
		return
	# ------------------------------------------------------------------------------
	# 
	def Generate_FLYA_cb(self):

		infowindo = Tkinter.Toplevel(self.top)
		infowindo.title('Generate Methyl FLYA Input')
		dp = directory_field(infowindo, 'Select Directory Location', 'Browse...')
		self.dirpath = dp.variable
		dp.frame.pack(side = 'top', anchor = 'w')
		explain = ('Specify name of output diectory and root of output files (not spaces)')
		w = Tkinter.Label(infowindo, text = explain, justify = 'left')
		w.pack(side = 'top', anchor = 'w')
		dn = tkutil.entry_field(infowindo, 'Directory Name: ', '', 15)
		self.dirname = dn.variable
		dn.frame.pack(side = 'top', anchor = 'w')
		fn = tkutil.entry_field(infowindo, 'File Name: ', '', 15)
		self.filename = fn.variable
		fn.frame.pack(side = 'top', anchor = 'w')
		
		br = tkutil.button_row(infowindo,
								('Ok', self.make_methylFLYA_cb))
		br.frame.pack(side = 'top', anchor = 'w')
	def close_child(self, w):
		w.destroy()
	def make_methylFLYA_cb(self):
		outdir = self.dirpath.get()+'/' + self.dirname.get() + '/'
		name = self.filename.get()

		if not os.path.exists(outdir):
			os.makedirs(outdir)
		s = self.get_settings()
		Iids = [res.replace('I','') for res in s.flya_atoms if res.split('.')[-1][0] == 'I']
		Lids = [res.replace('L','') for res in s.flya_atoms if res.split('.')[-1][0] == 'L']
		Vids = [res.replace('V','') for res in s.flya_atoms if res.split('.')[-1][0] == 'V']
		Mids = [res.replace('M','') for res in s.flya_atoms if res.split('.')[-1][0] == 'M']
		Aids = [res.replace('A','') for res in s.flya_atoms if res.split('.')[-1][0] == 'A']
		Tids = [res.replace('T','') for res in s.flya_atoms if res.split('.')[-1][0] == 'T']
		IDtrans = {'I':Iids,'L':Lids, 'V':Vids, 'M':Mids, 'A':Aids, 'T':Tids}
		HMQCpeaks = sorted(s.hmqc_spectrum.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
		HMQCnew = open(outdir + name + '_2DHMQC.peaks', 'w')
		HMBCout = open(outdir + name + '_HMBC-HMQC.peaks', 'w')
		HMQCnew.write('# Number of dimensions 2\n#FORMAT xeasy2D\n#INAME 1 C\n#INAME 2 H\n#SPECTRUM C13HSQC C H\n')
		HMBCout.write('# Number of dimensions 4\n#FORMAT xeasy4D\n#INAME 1 C1\n#INAME 2 H1\n#INAME 3 C2\n#INAME 4 H2\n#SPECTRUM HCcCH C1 H1 C2 H2\n')
		LVIDs = {}
		Freqs = {}
		LVpeaks = []
		for x in range(len(HMQCpeaks)):
			peak = HMQCpeaks[x]
			if peak.is_assigned != 0:
				# if len(peak.resonances()[0].group.symbol) == 1:
				if len(peak.note) < 1: 
					tkMessageBox.showinfo('Input Error', " Missing Peak Type\nfor %s\nUpdate list try again" %(peak.assignment))
					return
				if len(peak.note) >= 1 and peak.note.split()[0] not in ['I','L','V','M','A','T']:
					tkMessageBox.showinfo('Input Error', "Ambiguous Peak Type\nfor %s\nUpdate list try again" %(peak.assignment))
					return
				ptype = peak.note[0]
				trans = IDtrans[peak.note[0]]
				n = x + 1
				HMQCnew.write("%4.0f%8.3f%8.3f 1 U%19.3e         0 e   0  %8s  %8s\n" %(n, peak.frequency[0],peak.frequency[1], peak.data_height, trans[0], trans[0].replace('C','Q')))
				if ptype in ['L','V'] and ';' not in peak.note: trans.remove(trans[1])
				if ptype in ['L','V'] and 'C' in peak.note:
					LVIDs[peak.assignment] = trans[0]
					Freqs[peak.assignment] = peak.frequency

					LVpeaks.append(peak)
				trans.remove(trans[0])
		for y in range(len(LVpeaks)):
			peak = LVpeaks[y]
			n = y +1 
			pair = peak.note.split(';')[-1]
			HMBCout.write("%4.0f%8.3f%8.3f%8.3f%8.3f 1 U   1.000000E+00  0.000000E+00 e 0    %8s   %8s   %8s   %8s\n" %(n, peak.frequency[0],peak.frequency[1], Freqs[pair][0], Freqs[pair][1], LVIDs[peak.assignment], LVIDs[peak.assignment].replace('C','Q'),LVIDs[pair], LVIDs[pair].replace('C','Q')))

		HMQCnew.close()
		HMBCout.close()
		CCHnoesy = open(outdir + name + '_3D_NOESY.peaks','w')
		CCHnoesy.write('# Number of dimensions 3\n#FORMAT xeasy3D\n#INAME 1 C1\n#INAME 2 C2\n#INAME 3 H1\n#SPECTRUM CCNOESY3D C1 C2 H1\n')
		n = 0
		for noe in s.noesy_spectrum.peak_list():
			n = n +1 
			CCHnoesy.write("%4.0f%8.3f%8.3f%8.3f 1 U%19.3e         0 e   0 -         -         -\n" % (n,noe.frequency[0], noe.frequency[1], noe.frequency[2], noe.data_height))
		CCHnoesy.close()

		seqauto = open(outdir + name +'.seq','w')
		for (resi,resn) in self.sequence:
			seqauto.write('%s  %s\n' %(A_dict[resn], resi))
		seqauto.close()
		try:
			shutil.copy(s.pdb_path,'./MAGIC/')
		except OSError:
			pass
		tkMessageBox.showinfo('Finished', "Input files for methylFLYA\nSuccessfully generated\n%s" %(outdir))
		return

	# ------------------------------------------------------------------------------
	# 
	def Read_MAGIC_results_cb(self):
		s = self.get_settings()
		file_opt = options = {}
		options['defaultextension'] = '.list'
		options['filetypes'] = [('Magig Output List', '.list'), ('methylFLYA tab', '.tab')]
		options['title'] = 'Select Output File'

		path = tkFileDialog.askopenfilename(**file_opt)
		if path.split('.')[-1] == 'list':
			peak_list = open(path).readlines()
			if path.split('/')[-1] == 'hmqc.list':
				spectrum = s.hmqc_spectrum
			if path.split('/')[-1] == 'hmqc_iso.list':
				spectrum = s.hmqc_spectrum
			if path.split('/')[-1] == 'cch_iso.list':
				spectrum = s.noesy_spectrum 
			if path.split('/')[-1] == 'cch.list':
				spectrum = s.noesy_spectrum 
			save_path_new = str(spectrum.save_path+'MAGIC')
			save_new = open(save_path_new, 'w')
			save_content = open(spectrum.save_path, 'r')
			for line in save_content.readlines():
				line = line.replace('name ' + spectrum.name ,'name '+ spectrum.name + '_magic')
				if (line=='<end view>\n'): 
					break
				save_new.write(line)
			save_new.write('<end view>\n<ornament>\n<end ornament>\n<end spectrum>\n')
			save_new.close()
			self.session.open_spectrum(save_path_new)
			spectrums=self.session.project.spectrum_list()
			spectrum_ref_name=spectrum.name
			for spectrum in spectrums:
				if (spectrum.name == spectrum_ref_name+'_magic'):
					for peak in peak_list:
						self.create_peak(peak, spectrum)

		if path.split('.')[-1] == 'tab':
			spectrum = s.hmqc_spectrum
			intab = [line.strip() for line in open(path).readlines() if len(line.split()) >= 7 and 'Atom' not in line]
			
			save_path_new = str(spectrum.save_path+'FLYA')
			save_new = open(save_path_new, 'w')
			save_content = open(spectrum.save_path, 'r')
			for line in save_content.readlines():
				line = line.replace('name ' + spectrum.name ,'name '+ spectrum.name + '_FLYA')
				if (line=='<end view>\n'): 
					break
				save_new.write(line)
			save_new.write('<end view>\n<ornament>\n<end ornament>\n<end spectrum>\n')
			save_new.close()
			self.session.open_spectrum(save_path_new)
			spectrums=self.session.project.spectrum_list()
			spectrum_ref_name=spectrum.name
			for spectrum in spectrums:
				if (spectrum.name == spectrum_ref_name+'_FLYA'):
					self.pars_flya(intab, spectrum)

# ------------------------------------------------------------------------------
# 
	def parse_peak_line(self, line, dim):

		fields = string.split(line, None, dim + 1)
		if len(fields) < dim + 1:
			return None

		assignment = sputil.parse_assignment(fields[0])
		lbl = None
		if assignment == None or len(assignment) != dim:
			lbl = fields[0]
			assignment = None

		frequency = []
		try:
			for a in range(dim):
				f = pyutil.string_to_float(fields[a+1])
				if f == None:
					return None
				frequency.append(f)
		except:
			return None

		if len(fields) > dim + 1:
			note = fields[dim + 1].replace("'",'').strip()
			if 'NotAss' in note:
				note = ''
			else:
				note = note
		return (assignment, lbl, frequency, note)

# ------------------------------------------------------------------------------
# 
	def create_peak(self, line, spectrum):

		pinfo = self.parse_peak_line(line, spectrum.dimension)
		if pinfo:
			assignment, lbl, frequency, note = pinfo
			peak = spectrum.place_peak(frequency)
			self.move_peak_onto_spectrum(peak)
			assigned = 0
			if assignment != None:
				for a in range(spectrum.dimension):
					group_name, atom_name = assignment[a]
					if group_name or atom_name:
						peak.assign(a, group_name, atom_name)
				assigned = 1
			if assigned:
				peak.show_assignment_label()
			if note:
				peak.note = note
				if note.count(":") == 1:
					peak.color = 'green'
					peak.label.color = 'green'
				if note.count(":") >= 2:
					peak.color = 'gold'
					peak.label.color = 'gold'
				if note[0] == "-":
					peak.color = 'blue'
					peak.label.color = 'blue'

# ------------------------------------------------------------------------------
# 
	def pars_flya(self, intab, spectrum):

		ResDict = {}
		for res in intab :
			ResDict[AAA_dict[res.split()[1]] + res.split()[2] + res.split()[0].replace('Q','H')] = [float(res.split()[3]),res.split()[5]]
		for (resi,resn) in self.sequence:
			if resn in self.labeling.get():
				if self.dimethyl.state() == True or self.mono.state() == True:
					for me in Me_dic[resn]:
						assigment = []
						frequency = []
						note = ''
						for atom in me.split('-'):
							assigment.append(atom)
							if resn + str(resi) + atom in ResDict.keys():
								frequency.append(ResDict[resn + str(resi) + atom][0])
								note = note + ResDict[resn + str(resi) + atom][1] + ','
						if len(frequency) > 1:
							peak = spectrum.place_peak(frequency)
							self.move_peak_onto_spectrum(peak)
							assigned = 0
							peak.assign(0,resn + str(resi),assigment[0])
							peak.assign(1,resn + str(resi),assigment[1])
							assigned = 1
							peak.show_assignment_label()
							peak.note = note[:-1]
							if float(note.split(',')[0]) >= 80:
								peak.color = 'green'
								peak.label.color = 'green'
							if float(note.split(',')[0]) < 80:
								peak.color = 'gold'
								peak.label.color = 'gold'

# ------------------------------------------------------------------------------
# 
	def move_peak_onto_spectrum(self, peak):

		freq = peak.frequency
		pos = sputil.alias_onto_spectrum(freq, peak.spectrum)
		if pos != freq:
			peak.position = pos
			peak.alias = pyutil.subtract_tuples(freq, pos)


# ------------------------------------------------------------------------------
# 
class directory_field:
	def __init__(self, parent, title, cache_name,browse_button = 1, save = 0, width = 20):
		self.title = title
		self.cache_name = cache_name
		self.frame = Tkinter.Frame(parent)
		e = tkutil.entry_field(self.frame, title, '', width)
		self.entry = e
		self.variable = e.variable
		e.frame.pack(side = 'left')
		if browse_button:
			self.save = save
			b = Tkinter.Button(self.frame, text = "Browse ...", command = self.file_browse_cb)
			b.pack(side = 'left')
	# --------------------------------------------------------------------------
	# Pop up a file browsing dialog and set the variable to the chosen file.
	#
	def file_browse_cb(self):
		path = tkFileDialog.askdirectory()
		if path:
			self.variable.set(path)
			self.entry.show_end()




def show_dialog(session):
	sputil.the_dialog(methyl_dialog,session).show_window(1)
