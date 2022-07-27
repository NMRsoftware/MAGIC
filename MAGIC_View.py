import Tkinter
from Tkinter import *
import types
import math as m
import pyutil
import sparky
import myseq
import sputil
import tkutil
import os
import expectedpeaks
import tkMessageBox
import tkFileDialog
import string

# ------------------------------------------------------------------------------
#
# Updated by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
#
# Last updates: June 1, 2022
#
#
# ------------------------------------------------------------------------------
#
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}  
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
 "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L",
 "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
 "TRP": "W", "TYR": "Y", "VAL": 'V' }
class spectrum_menu(tkutil.option_menu):

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

class obj:
  pass

class Assignment_Progress_dialog(tkutil.Dialog, tkutil.Stoppable):

  # ------------------------------------------------------------------------------
  #
  def __init__(self, session):

    self.session = session
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

      return  self.close_cb
    self.sequence = sequence
    self.selection_notice = None
    self.save_file = session.project.save_path.split('/')[-1].replace('.proj', '.txt')
    tkutil.Dialog.__init__(self, session.tk, 'MAGIC-View')
    explain = ('Pymol Visualization and Assignment Accuracy Statistics\n') 
    w = Tkinter.Label(self.top, text = explain, justify = 'left')
    w.pack(side = 'top', anchor = 'w')
    explain = ('PDB does NOT need H')
    w = Tkinter.Label(self.top, text = explain, justify = 'left')
    w.pack(side = 'top', anchor = 'w')

    ep = tkutil.file_field2(self.top, 'PDB file', 'Browse...', file_type=[('Protein Data Bank File', '.pdb')], default_ext='.pdb')
    self.pdb_path = ep.variable
    ep.frame.pack(side = 'top', anchor = 'w')

    lb = tkutil.entry_field(self.top, 'Labeled Methyls: ', 'ILVMAT', 10)
    lb.frame.pack(side = 'top', anchor = 'w')
    self.labeling = lb.variable

    explain = ('Specify which chain should be used')
    w = Tkinter.Label(self.top, text = explain, justify = 'left')
    w.pack(side = 'top', anchor = 'w')
    ch = tkutil.entry_field(self.top, 'Use Chain(s): ', 'A', 5)
    self.chains = ch.variable
    ch.frame.pack(side = 'top', anchor = 'w')

    explain = ('Select 2D spectra that contain assignments')
    w = Tkinter.Label(self.top, text = explain, justify = 'left')
    w.pack(side = 'top', anchor = 'w')    
    self.sc = self.spectrum_choice_table(self.top)
    self.sc.pack(side = 'top', anchor = 'w')

  # Compare Assignments Box
    self.lf2 = Tkinter.LabelFrame(self.top, text='Compare Assignments')
    self.lf2.pack(fill=X, pady=5)
    self.pmframe= Tkinter.Frame(self.lf2)
    self.pmframe.pack(pady=5)
    # First Spectrum
    self.sc_ref_spectrum = spectrum_menu(session, self.pmframe, 'Ref Spectrum: ')
    self.sc_ref_spectrum.frame.pack(side ='top',padx=3)
  # Second Spectrum  
    self.sc_spectrum_2 = spectrum_menu(session, self.pmframe, 'Query Spectrum: ')
    self.sc_spectrum_2.frame.pack(side = 'top', anchor = 'w')

    sl = tkutil.scrolling_list(self.top, 'Assignment Report', 10)
    sl.frame.pack(fill = 'both', expand = 1)
    self.summary_list = sl

  # PPM Tolerances for filtering noesy and finding geminal pairs
    er = tkutil.entry_row(self.top, 'PPM tolerance: ',
                                    ('1H', '0.01', 5),
                                    ('13C', '0.10', 5))
    er.frame.pack(side = 'top', anchor = 'w',padx=2)
    self.ppm_range = er

    br = tkutil.button_row(self.top,
                           ('Stats', self.show_summary2),
                           ('Compare', self.show_summary),
                           ('Save', self.save_pymol_cb),
                           ('Close', self.close_cb))
    br.frame.pack(side = 'top', anchor = 'w')

  # ------------------------------------------------------------------------------
  #
  def spectrum_choice_table(self, parent):

    headings = ('')
    st = sputil.spectrum_table(self.session, parent, headings, self.add_spectrum, self.remove_spectrum)
    st.spectrum_to_checkbutton = {}
    st.chosen_spectra = []
    st.spectrum_epeak_menus = {}
    st.axis_order_menu = {}
    self.spectrum_table = st

    spectra = self.session.project.spectrum_list()
    for spectrum in spectra:
      if spectrum.dimension == 2:
        st.add_spectrum(spectrum)
    return st.frame
# -----------------------------------------------------------------------------
#
  def add_spectrum(self, spectrum, table, row):

    pat_name, pat_axes = expectedpeaks.recall_pattern(spectrum)
    # Make spectrum check button
    cb = tkutil.checkbutton(table.frame, spectrum.name, 0)
    choose_cb = pyutil.precompose(sputil.choose_spectrum_cb, spectrum, table.chosen_spectra)
    cb.add_callback(choose_cb)
    cb.button.grid(row = row, column = 0, sticky = 'w')
    if pat_name:
      cb.set_state(1)
    table.spectrum_to_checkbutton[spectrum] = cb
# -----------------------------------------------------------------------------
#
  def remove_spectrum(self, spectrum, table):

    cb = table.spectrum_to_checkbutton[spectrum]
    cb.set_state(0)
    cb.button.destroy()
    del table.spectrum_to_checkbutton[spectrum]

    table.spectrum_epeak_menus[spectrum].frame.destroy()
    del table.spectrum_epeak_menus[spectrum]

    table.axis_order_menu[spectrum].frame.destroy()
    del table.axis_order_menu[spectrum]
  # ------------------------------------------------------------------------------
  # 
  def get_settings(self):
  
    settings = pyutil.generic_class()
    Labels = []
    for (resi,resn) in self.sequence:
      Labels.append(resn + str(resi))
    spectra = self.spectrum_table.chosen_spectra
    settings.spectrum_list=[spectrum for spectrum in spectra]
    settings.pdb_path = self.pdb_path.get()
    settings.chains = self.chains.get().replace(',','')
    settings.labeling = self.labeling.get()
    settings.sequence = self.sequence
    settings.ref_spectrum = self.sc_ref_spectrum.spectrum()
    settings.spectrum_2 = self.sc_spectrum_2.spectrum()
    settings.Htol = float(self.ppm_range.variables[0].get())
    settings.Ctol = float(self.ppm_range.variables[1].get())
    settings.Labels = Labels

    return settings

  # ------------------------------------------------------------------------------
  
  def save_pymol_cb(self):
    file_opt = options = {}
    options['defaultextension'] = '.txt'
    options['filetypes'] = [('text file', '.txt')]
    options['title'] = 'Save Pymol/Chimera X Script'
    options['initialdir'] = self.session.project.save_path
    options['initialfile'] = self.save_file

    path = tkFileDialog.asksaveasfilename(**file_opt)
    if path:
      self.Generate_Pymol_cb(path)

# -----------------------------------------------------------------------------
#
  def Generate_Pymol_cb(self,path):

    s = self.get_settings()
 
    ## Read sequence file and generate list of available sequence indexes for NH and CH residues. 
    ## This will be used to determine missing assignments and generate selection statements.

    pdb_name = s.pdb_path.split('/')[-1].split('.')[0]
    methyls,backbone, aromatic = [], [], []
    showmet = 'show :'
    showaro = 'show :'
    sdict = {}
    for (resi,resn) in self.sequence:
      sdict[resi] = resn+str(resi)
      if resn in ['I', 'L', 'V', 'M', 'A', 'T']:
          methyls.append(str(resi))
      if resn in ['Y', 'F']:
          aromatic.append(str(resi))
      if resn != 'P':
        backbone.append(str(resi))

    cmx_header = "open %s\nshow nucleic\nhide protein|solvent|H\nsurface hide\nstyle (protein|nucleic|solvent) & @@draw_mode=0 stick\ncartoon style modeHelix tube sides 20\ngraphics silhouettes tru\ncartoon style radius 1.5\nset bgColor white\n#style :ala,leu,ile,val,met,thr,phe,tyr ball\ncolor all gray(150)\n" %(s.pdb_path)
    pml_header  = "load %s\nset ray_opaque_background, 0\nset depth_cue, off\nbg_color white\nset label_color, black\nset ray_shadows, 0\
    \nset orthoscopic, on\nhide everything, all\nshow cartoon, %s and chain %s and resi %s-%s\n" %(s.pdb_path, pdb_name, s.chains, backbone[0], backbone[-1])
    
    plm_out = open(path.replace('.txt', '_pymol.pml'),'w')
    cmx_out = open(path.replace('.txt','_chimera.cxc'),'w')

    plm_out.write(pml_header)
    cmx_out.write(cmx_header)
    plm_out.write("color gray70, %s and chain %s and resi %s-%s\n" %( pdb_name, s.chains, backbone[0], backbone[-1]))
    mn = 1
    for spectrum in s.spectrum_list:
      mn+=1
      cmx_out.write('combine #1 modelId %s name %s\ncartoon style radius 1.5\n' %(mn, spectrum.name))
      plm_out.write('create %s, %s and chain %s\n'%(spectrum.name, pdb_name, s.chains))
      assigned = []
      tentative = []
      used = []
      for label in spectrum.label_list():
        peak = label.peak
        if '?' not in peak.assignment:
          if (peak.resonances()[0].group.number, peak.resonances()[0].group.name[0]) in self.sequence:
            if 'green' in peak.color.lower() or 'green' in label.color.lower():
                if int(peak.resonances()[0].group.number) not in assigned:
                  assigned.append(int(peak.resonances()[0].group.number))
                  used.append(str(peak.resonances()[0].group.number))
            if 'gold' in peak.color.lower() or 'gold' in label.color.lower():
                if int(peak.resonances()[0].group.number) not in tentative:
                  tentative.append(int(peak.resonances()[0].group.number))
                  used.append(str(peak.resonances()[0].group.number))
            if 'yellow' in peak.color.lower() or 'yellow' in label.color.lower():
                if int(peak.resonances()[0].group.number) not in tentative:
                  tentative.append(int(peak.resonances()[0].group.number))
                  used.append(str(peak.resonances()[0].group.number))

      assigned = sorted(assigned)
      tentative = sorted(tentative)


### Prepare output, Pymol does not support resi list containing more than 250 + values 
### so any color/selection statements must be shorter than that or pymol will crash 


      if '13C' in spectrum.nuclei:
        if spectrum.region[0][0] < 40:
          assignments = methyls
          showstr = 'show #%s:ala,val,leu,met,ile,thr\n' %mn
          pmlshow = 'show sticks, %s and resn ALA+VAL+LEU+MET+ILE+THR\nhide sticks, name N+C\nhide sticks, elem H\ncolor red, elem O\ncolor gold, elem S\n' %spectrum.name
        if spectrum.region[0][0] > 40:
          assignments = aromatic
          showstr = 'show #%s:phe,thr\n' %mn
          pmlshow = 'show sticks, %s and resn PHE+TYR\nhide sticks, name N+C\nhide sticks, elem H\ncolor red, elem O\n' %spectrum.name
        unassigned = []
        for resi in assignments:
          if resi not in used: unassigned.append(int(resi))
        if len(assigned) != 0:
          for x in range(0,len(assigned),50):
            i = x
            plmgreen = 'color smudge, %s and resi ' %spectrum.name
            cmxgreen = 'color #%s:' %mn
            for j in range(50):
              plmgreen = plmgreen + str(assigned[i]) + '+'
              cmxgreen = cmxgreen + str(assigned[i]) + ','
              i=i+1
              if i== len(assigned): break
            plm_out.write(plmgreen[:-1] + '\n')
            cmx_out.write(cmxgreen[:-1] + ' green target ac\n')
        if len(tentative) != 0:
          tassign = '##### Tentative Assigments %d: ' %(len(tentative))
          for x in range(0,len(tentative),50):
            i = x
            plmgold = 'color gold, %s and resi ' %spectrum.name
            cmxgold = 'color #%s:' %mn
            for j in range(50):
              tassign = tassign + sdict[tentative[i]] + ' '
              plmgold = plmgold + str(tentative[i]) + '+'
              cmxgold = cmxgold + str(tentative[i]) + ','
              i=i+1
              if i== len(tentative): break
            plm_out.write(plmgold[:-1] + '\n')
            cmx_out.write(cmxgold[:-1] + ' gold target ac\n')
        if len(unassigned) != 0:
          missing = '##### Unassigned Residues %d: ' %(len(unassigned))
          for x in range(0,len(unassigned),50):
            i = x
            plmred = 'color red, %s and resi ' %spectrum.name
            cmxred = 'color #%s:' %mn
            for j in range(50):
              missing = missing + sdict[unassigned[i]] + ' '
              plmred = plmred + str(unassigned[i]) + '+'
              cmxred = cmxred + str(unassigned[i]) + ','
              i=i+1
              if i== len(unassigned): break
            plm_out.write(plmred[:-1] + '\n')
            cmx_out.write(cmxred[:-1] + ' red target ac\n')
        cmx_out.write(showstr)
        plm_out.write(pmlshow)
        if len(tentative) != 0:
          plm_out.write(tassign+ '\n')
          cmx_out.write(tassign+ '\n')
        if len(unassigned) != 0:
          plm_out.write(missing+ '\n\n')
          cmx_out.write(missing+ '\n\n')

      if '15N' in spectrum.nuclei:
        unassigned = []
        for resi in backbone:
          if resi not in used: unassigned.append(int(resi))
        if len(assigned) != 0:
          cmx_greens = 'color #' + str(mn) + ':' + str(assigned[0])
          ngreens = 'color smudge, %s and resi %s' %(spectrum.name, str(assigned[0]))
          for g in range(len(assigned))[:-1]:
            idx = assigned[g]
            idx2 = assigned[g+1]
            if idx2 != idx +1:
              cmx_greens = cmx_greens + '-' + str(idx) + ',' + str(idx2)
              ngreens = ngreens + '-' + str(idx) + '+' + str(idx2)
          if assigned[-2] +1  == assigned[-1]:
              ngreens = ngreens + '-' + str(assigned[-1])
              cmx_greens = cmx_greens  + '-' + str(assigned[-1])
          plm_out.write(ngreens + '\n')
          cmx_out.write(cmx_greens + ' green target c\n')
        if len(tentative) != 0:
          tassign = '##### Tentative Assigments %d: ' %(len(tentative))
          ngolds = 'color gold, %s and resi %s' %(spectrum.name, str(tentative[0]))
          cmx_golds = 'color #' + str(mn) + ':' + str(tentative[0])
          for t in range(len(tentative))[:-1]:
            tassign = tassign + sdict[tentative[t]] + ' '
            idx = tentative[t]
            idx2 = tentative[t+1]
            if idx2 != idx +1:
              ngolds = ngolds + '-' + str(idx) + '+' + str(idx2)
              cmx_golds = cmx_golds + '-' + str(idx) + ',' + str(idx2)
          if tentative[-2] +1  == tentative[-1]:
              ngolds = ngolds + '-' + str(tentative[-1])
              cmx_golds = cmx_golds + '-' + str(tentative[-1])
              tassign = tassign + sdict[tentative[-1]]
          plm_out.write(ngolds+ '\n')
          cmx_out.write(cmx_golds + ' gold target c\n')
        if len(unassigned) != 0:
          nreds = 'color red,  %s and resi %s' %(spectrum.name, str(unassigned[0]))
          cmx_reds = 'color #' + str(mn) + ':' + str(unassigned[0])
          missing = '##### Unassigned NH %d: ' %(len(unassigned))
          for r in range(len(unassigned))[:-1]:
            missing = missing + sdict[unassigned[r]] + ' '
            idx = unassigned[r]
            idx2 = unassigned[r+1]
            if idx2 != idx +1:
              nreds = nreds + '-' + str(idx) + '+' + str(idx2)
              cmx_reds = cmx_reds + '-' + str(idx) + ',' + str(idx2)
          if unassigned[-2] +1  == unassigned[-1]:
            nreds = nreds + '-' + str(unassigned[-1])
            cmx_reds = cmx_reds + '-' + str(unassigned[-1])
          missing = missing + sdict[unassigned[-1]]
          plm_out.write(nreds+ '\n')
          cmx_out.write(cmx_reds + ' red target c\n')
        if len(tentative) != 0:
          plm_out.write(tassign+ '\n')
          cmx_out.write(tassign+ '\n')
        if len(unassigned) != 0:
          plm_out.write(missing+ '\n\n')
          cmx_out.write(missing+ '\n\n')
    cmx_out.write('hide #1 target ac\n')
    plm_out.write("color gray70, %s and chain %s and resi %s-%s\ncolor gray70, resn PRO\n" %( pdb_name, s.chains, backbone[0], backbone[-1]))
    plm_out.write("hide everything, %s\n" %pdb_name)
    cmx_out.close()
    plm_out.close()

  def show_summary(self):
    s = self.get_settings()
    self.summary_list.clear()
    self.stoppable_loop('shifts', 100)
    sA,sC,sD,sE,sF,sG,sH,sI,sS,sK,sL,sM,sN,sP,sQ,sR,sS,sT,sV,sW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    uA,uC,uD,uE,uF,uG,uH,uI,uS,uK,uL,uM,uN,uP,uQ,uR,uS,uT,uV,uW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    aA,aC,aD,aE,aF,aG,aH,aI,aS,aK,aL,aM,aN,aP,aQ,aR,aS,aT,aV,aW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    cA,cC,cD,cE,cF,cG,cH,cI,cS,cK,cL,cM,cN,cP,cQ,cR,cS,cT,cV,cW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    isA,isC,isD,isE,isF,isG,isH,isI,isS,isK,isL,isM,isN,isP,isQ,isR,isS,isT,isV,isW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    imA,imC,imD,imE,imF,imG,imH,imI,imS,imK,imL,imM,imN,imP,imQ,imR,imS,imT,imV,imW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    # print s.Labels
    aTotal = 0
    good, bad = self.Check_Assignments_cb()
    for assignment in s.Labels:
      if assignment[0] == 'I':sI+=1
      if assignment[0] == 'L':sL+=1
      if assignment[0] == 'V':sV+=1
      if assignment[0] == 'M':sM+=1
      if assignment[0] == 'A':sA+=1
      if assignment[0] == 'T':sT+=1
    Correct, Incorrect,Used = [], [], []
    for peak in s.spectrum_2.peak_list():
      if '?' not in peak.assignment:
        restype = peak.resonances()[0].group.symbol
        if restype == '-':restype = peak.assignment[0]
        if peak.note.count(":") == 1:
          if restype == 'I': 
            uI += 1
            if peak in bad: isI+=1
          if restype == 'L': 
            uL += 1
            if peak in bad: isL+=1
          if restype == 'V': 
            uV += 1
            if peak in bad: isV+=1
          if restype == 'M': 
            uM += 1
            if peak in bad: isM+=1
          if restype == 'A': 
            uA += 1
            if peak in bad: isA+=1
          if restype == 'T': 
            uT += 1
            if peak in bad: isT+=1
        if peak.note.count(":") >= 2:
          aTotal = aTotal +1
          if 'I'in restype: 
            aI+=1
            if peak in bad: imI+=1
          if 'L'in restype: 
            aL+= 1
            if peak in bad: imL+=1
          if 'V'in restype: 
            aV+= 1
            if peak in bad: imV+=1
          if 'M'in restype: 
            aM+= 1
            if peak in bad: imM+=1
          if 'A'in restype: 
            aA+= 1
            if peak in bad: imA+=1
          if 'T'in restype: 
            aT+= 1
            if peak in bad: imT+=1

    for peak in good:
      if 'I'in peak.assignment[0]: cI+=1
      if 'L'in peak.assignment[0]: cL+=1
      if 'V'in peak.assignment[0]: cV+=1
      if 'M'in peak.assignment[0]: cM+=1
      if 'A'in peak.assignment[0]: cA+=1
      if 'T'in peak.assignment[0]: cT+=1
    print len(good)
    uTotal = uI+uL+uV+uM+uA+uT
    cTotal = cI+cL+cV+cM+cA+cT
    sTotal = sI+2*sL+2*sV+sM+sA+sT
    icTotal = isI+isL+isV+isM+isA+isT + imI+imL+imV+imM+imA+imT
    total = uI+uL+uV+uM+uA+uT+aI+aL+aV+aM+aA+aT
    self.summary_list.append(('Assigned: %3d (%2.1f%%)' %(total, 100*float(total)/sTotal)))
    self.summary_list.append(('Correct: %3d (%2.1f%%)' %(cTotal, 100*float(cTotal)/sTotal)))
    self.summary_list.append(('Incorrect: %3d(%2.1f%%)' %(icTotal, 100*float(icTotal)/sTotal)))
    self.summary_list.append('              Assignment Opt               Incorrect')
    self.summary_list.append('   Expected  Single  Multiple  Correct  Single  Multiple ')
    self.summary_list.append('I    %3d       %3d     %3d      %3d      %3d     %3d' %(  sI, uI, aI, cI, isI, imI))
    self.summary_list.append('L    %3d       %3d     %3d      %3d      %3d     %3d' %(2*sL, uL, aL, cL, isL, imL))
    self.summary_list.append('V    %3d       %3d     %3d      %3d      %3d     %3d' %(2*sV, uV, aV, cV, isV, imV))
    self.summary_list.append('M    %3d       %3d     %3d      %3d      %3d     %3d' %(  sM, uM, aM, cM, isM, imM))
    self.summary_list.append('A    %3d       %3d     %3d      %3d      %3d     %3d' %(  sA, uA, aA, cA, isA, imA))
    self.summary_list.append('T    %3d       %3d     %3d      %3d      %3d     %3d' %(  sT, uT, aI, cT, isT, imT))
  def show_summary2(self):
    s = self.get_settings()
    self.summary_list.clear()
    self.stoppable_loop('shifts', 100)
    sA,sC,sD,sE,sF,sG,sH,sI,sS,sK,sL,sM,sN,sP,sQ,sR,sS,sT,sV,sW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    uA,uC,uD,uE,uF,uG,uH,uI,uS,uK,uL,uM,uN,uP,uQ,uR,uS,uT,uV,uW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    aA,aC,aD,aE,aF,aG,aH,aI,aS,aK,aL,aM,aN,aP,aQ,aR,aS,aT,aV,aW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    cA,cC,cD,cE,cF,cG,cH,cI,cS,cK,cL,cM,cN,cP,cQ,cR,cS,cT,cV,cW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    isA,isC,isD,isE,isF,isG,isH,isI,isS,isK,isL,isM,isN,isP,isQ,isR,isS,isT,isV,isW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    imA,imC,imD,imE,imF,imG,imH,imI,imS,imK,imL,imM,imN,imP,imQ,imR,imS,imT,imV,imW = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    # print s.Labels
    aTotal = 0
    for assignment in s.Labels:
      if assignment[0] == 'I':sI+=1
      if assignment[0] == 'L':sL+=1
      if assignment[0] == 'V':sV+=1
      if assignment[0] == 'M':sM+=1
      if assignment[0] == 'A':sA+=1
      if assignment[0] == 'T':sT+=1
    Correct, Incorrect,Used = [], [], []
    for peak in s.ref_spectrum.peak_list():
      if '?' not in peak.assignment:
        restype = peak.resonances()[0].group.symbol
        if restype == '-':restype = peak.assignment[0]
        if peak.note.count(":") == 1:
          if restype == 'I': uI += 1
          if restype == 'L': uL += 1
          if restype == 'V': uV += 1
          if restype == 'M': uM += 1
          if restype == 'A': uA += 1
          if restype == 'T': uT += 1
        if peak.note.count(":") >= 2:
          aTotal = aTotal +1
          if 'I'in restype: aI+=1
          if 'L'in restype: aL+= 1
          if 'V'in restype: aV+= 1
          if 'M'in restype: aM+= 1
          if 'A'in restype: aA+= 1
          if 'T'in restype: aT+= 1
    uTotal = uI+uL+uV+uM+uA+uT
    cTotal = cI+cL+cV+cM+cA+cT
    sTotal = sI+2*sL+2*sV+sM+sA+sT
    total = uI+uL+uV+uM+uA+uT+aI+aL+aV+aM+aA+aT
    self.summary_list.append(('Assigned: %3d (%2.1f%%)' %(total, 100*float(total)/sTotal)))
    self.summary_list.append('              Assignment Opt   ')
    self.summary_list.append('   Expected  Single  Multiple  ')
    self.summary_list.append('I    %3d       %3d     %3d   ' %(  sI, uI, aI))
    self.summary_list.append('L    %3d       %3d     %3d   ' %(2*sL, uL, aL))
    self.summary_list.append('V    %3d       %3d     %3d   ' %(2*sV, uV, aV))
    self.summary_list.append('M    %3d       %3d     %3d   ' %(  sM, uM, aM))
    self.summary_list.append('A    %3d       %3d     %3d   ' %(  sA, uA, aA))
    self.summary_list.append('T    %3d       %3d     %3d   ' %(  sT, uT, aI))
  # -----------------------------------------------------------------------------
  #
  def Check_Assignments_cb(self):
    s = self.get_settings()
    RefPeaks = sorted(s.ref_spectrum.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
    CompPeaks = sorted(s.spectrum_2.peak_list(), key = lambda x: (x.assignment, x.frequency[0]))
    if s.ref_spectrum.name == s.spectrum_2.name:
      tkMessageBox.showinfo('Input Error', "Please select two different spectra")
      return
    if s.ref_spectrum.name != s.spectrum_2.name:
      gc = 0 
      bc = 0 
      good, bad, used = [], [], []
      for peak in CompPeaks:
        if '?' not in peak.assignment:
          print peak.assignment
          Matches = [rpeak for rpeak in RefPeaks if (abs(rpeak.frequency[0] -peak.frequency[0]) < s.Ctol) and (abs(rpeak.frequency[1] -rpeak.frequency[1]) < s.Htol) and rpeak not in used]
          diff = []
          if len(Matches) >= 1:
            for rpeak in Matches: 
              if '?' not in rpeak.assignment:
                if rpeak.resonances()[0].group.name == peak.resonances()[0].group.name: diff.append(-1)
                if rpeak.resonances()[0].group.name != peak.resonances()[0].group.name:
                  diff.append(abs(peak.frequency[0] - rpeak.frequency[0]) + abs(peak.frequency[1] - rpeak.frequency[1]))
              if '?' in rpeak.assignment:
                diff.append(abs(peak.frequency[0] - rpeak.frequency[0]) + abs(peak.frequency[1] - rpeak.frequency[1]))
            peak2 = Matches[diff.index(min(diff))]
            if '?' not in peak2.assignment:
              if peak.note.count(":") == 1:
                if peak.resonances()[0].group.name == peak2.resonances()[0].group.name:
                  good.append(peak)
                  gc += 1
                  print 'single Correct %d' %gc
                if peak.resonances()[0].group.name != peak2.resonances()[0].group.name:
                  bad.append(peak)
                  bc += 1
                  print 'single Incorrect %d' %bc
                  peak.color = 'red'
                used.append(peak2)
              if peak.note.count(":") >= 2:
                note = peak.note.split('{')[1]
                options = []
                for opt in note.split(','):
                  options.append(opt.split('C')[0].replace(' ',''))
                if peak2.resonances()[0].group.name in options:
                  peak.assign(0, peak2.resonances()[0].group.name, peak2.resonances()[0].atom.name)
                  peak.assign(1, peak2.resonances()[1].group.name, peak2.resonances()[1].atom.name)
                  good.append(peak)
                  gc += 1
                  print 'multiple Correct %d' %gc
                if peak2.resonances()[0].group.name not in options:
                  bad.append(peak)
                  bc += 1
                  print 'multiple Incorrect %d' %bc
                  peak.color = 'red'
                used.append(peak2)
      print gc
    return good, bad

  # ------------------------------------------------------------------------------

def show_dialog(session):
  sputil.the_dialog(Assignment_Progress_dialog,session).show_window(1)

