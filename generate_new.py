import types
import re
import numpy as np
import math as m
import datetime as time
import sys
import os
import multiprocessing as mp
import cPickle
import shutil
import resource

if str(sys.argv[1])=='-h':
  print 'python generate.py [2D peak list file name]'
  print '                   [construct sequence file name, fasta format]'
  print '                   [labeling, e.g. nAILMTV]'
  print '                   [starting sequence number, e.g. 48]'
  print '                   [rename: y or n]'
  print '                   [short mixing time CCH peak list] if available'
  exit()

class obj:
  pass
result=obj()    
 
labeling=str(sys.argv[3])
HMQC_file=open('./'+str(sys.argv[1]),'r')
HMQC=HMQC_file.readlines()
HMQC_file.close()
Seq_file=open('./'+str(sys.argv[2]),'r')
SEQ=Seq_file.read()
Seq_file.close()
seq_starter=int(sys.argv[4])
rename_status=str(sys.argv[5])

flag_geminal=0
try:
  CCH_file=open('./'+str(sys.argv[6]),'r')
  CCH=CCH_file.readlines()
  CCH_file.close()
  flag_geminal=1
except IndexError:pass
SEQ=SEQ.replace('\n','')
SEQ=SEQ.replace(' ','')
n,A,I,L,M,T,V=0,0,0,0,0,0,0
seq_output=open('./seq.auto','w')
for i in range(len(SEQ)):
  if 'n' in labeling:
    seq_output.write(str(SEQ[i])+str(seq_starter+i)+'\n')
    n+=1
    
  elif (SEQ[i]=='A' and SEQ[i] in labeling):
    seq_output.write('A'+str(seq_starter+i)+'\n')
    A+=1
  elif (SEQ[i]=='I' and SEQ[i] in labeling):
    seq_output.write('I'+str(seq_starter+i)+'\n')
    I+=1
  elif (SEQ[i]=='L' and SEQ[i] in labeling):
    seq_output.write('L'+str(seq_starter+i)+'\n')
    L+=1
  elif (SEQ[i]=='M' and SEQ[i] in labeling):
    seq_output.write('M'+str(seq_starter+i)+'\n')
    M+=1
  elif (SEQ[i]=='T' and SEQ[i] in labeling):
    seq_output.write('T'+str(seq_starter+i)+'\n')
    T+=1
  elif (SEQ[i]=='V' and SEQ[i] in labeling):
    seq_output.write('V'+str(seq_starter+i)+'\n')
    V+=1
seq_output.close()
print '###################################'
print 'Number of methyls: '+str(A+I+2*L+M+T+2*V)
print 'A: ',A,'\nI: ',I,'\nL: ',L,'\nM: ',M,'\nT: ',T,'\nV: ',V
seq_output.close()

sdAh=0.28
sdAc=1.79
muAh=1.36
muAc=18.97

sdIh=0.29
sdIc=1.67
muIh=0.68
muIc=13.40

sdLh=0.28
sdLc=1.70
muLh=0.74
muLc=24.365

sdVh=0.28
sdVc=1.54
muVh=0.815
muVc=21.395

sdMh=0.40
sdMc=1.70
muMh=1.89
muMc=17.11

sdTh=0.22
sdTc=1.11
muTh=1.14
muTc=21.55

peak_type=np.zeros((len(HMQC),6))
for i in range(len(HMQC)):
  line=HMQC[i]
  Wh=float(line.split()[2])
  Wc=float(line.split()[1])
  PA,PI,PL,PV,PM,PT=0,0,0,0,0,0
  if 'A' in labeling:PA=round(m.exp(-((m.pow((Wc-muAc),2)/(2*m.pow(sdAc,2)))+(m.pow((Wh-muAh),2)/(2*m.pow(sdAh,2))))),5)
  if 'I' in labeling:PI=round(m.exp(-((m.pow((Wc-muIc),2)/(2*m.pow(sdIc,2)))+(m.pow((Wh-muIh),2)/(2*m.pow(sdIh,2))))),5)
  if 'L' in labeling:PL=round(m.exp(-((m.pow((Wc-muLc),2)/(2*m.pow(sdLc,2)))+(m.pow((Wh-muLh),2)/(2*m.pow(sdLh,2))))),5)
  if 'V' in labeling:PV=round(m.exp(-((m.pow((Wc-muVc),2)/(2*m.pow(sdVc,2)))+(m.pow((Wh-muVh),2)/(2*m.pow(sdVh,2))))),5)
  if 'M' in labeling:PM=round(m.exp(-((m.pow((Wc-muMc),2)/(2*m.pow(sdMc,2)))+(m.pow((Wh-muMh),2)/(2*m.pow(sdMh,2))))),5)
  if 'T' in labeling:PT=round(m.exp(-((m.pow((Wc-muTc),2)/(2*m.pow(sdTc,2)))+(m.pow((Wh-muTh),2)/(2*m.pow(sdTh,2))))),5) 
  peak_type[i,0]=PA
  peak_type[i,1]=PI
  peak_type[i,2]=PL
  peak_type[i,3]=PV
  peak_type[i,4]=PM
  peak_type[i,5]=PT
HMQC_newfile=open('./new_'+str(sys.argv[1]),'w')
tot_0=0
tot_1=0
tot_2=0
tot_3=0
tot_types=0
name=0
for i in range(peak_type.shape[0]):
  name+=1
  tot=peak_type[i,:].sum()
  text=''
  if tot==0:pass
  else:
   for j in range(6):
    if (((j==1 or j==4) and (peak_type[i,j]/float(tot))>=0.2) or
        ((j==2 or j==3) and (peak_type[i,j]/float(tot))>=0.1) or
        ((j==0 or j==5) and (peak_type[i,j]/float(tot))>=0.01)
        ):
    #if (peak_type[i,j])>0.01:
      if (j==0 and 'A' in labeling):text+='A'
      elif (j==1 and 'I' in labeling):text+='I'
      elif (j==2 and 'L' in labeling):text+='L'
      elif (j==3 and 'V' in labeling):text+='V'
      elif (j==4 and 'M' in labeling):text+='M'
      elif (j==5 and 'T' in labeling):text+='T'
    else:pass
  tot_types+=len(text)
  if text=='':text='!!!!!'
  else:pass
  #if HMQC[i].split()[0][0] in text:pass
  #else:print HMQC[i].split()[0],text
  if rename_status=='y':peak_name=str(name)
  else:peak_name=HMQC[i].split()[0]
  spaces=len('  '+peak_name)
  HMQC_newfile.write('  '+peak_name+(14-spaces)*' '+HMQC[i].split()[1]+'   '+HMQC[i].split()[2]+'\t\t'+text+'\n')
#HMQC_newfile.write('tot_0='+str(tot_0)+', tot_1='+str(tot_1)+', tot_2='+str(tot_2)+', tot_3='+str(tot_3))
print 'Average number of methyl type per peak: ', round(tot_types/float(len(HMQC)),2)
HMQC_newfile.close()
############Automatic first guess for geminal pairing according to short mixing time CCH noesy#############
if flag_geminal==1:
  HMQC_newfile=open('./new_'+str(sys.argv[1]),'r')
  HMQC=HMQC_newfile.readlines()
  HMQC_newfile.close()
  os.remove('./new_'+str(sys.argv[1]))
  HMQC_newfile=open('./new_'+str(sys.argv[1]),'w')
  CHH_geminal=open('./CCH_geminal_assigned.list','w')
  w1tol,w2tol,w3tol=0.1,0.1,0.01
  for line_hmqc in HMQC:
      split_hmqc=line_hmqc.split()
      w,wnoe=float(split_hmqc[1]),float(split_hmqc[1])
      wh=float(split_hmqc[2])
      setattr(result, str(split_hmqc[0])+'_freq', (w,wh))
      noesy_peaks=[]
      for line_noe in CCH:
        split_noe=line_noe.split()	
        w1=float(split_noe[1])	
        w2=float(split_noe[2])
        w3=float(split_noe[3])
        if (w2>(w-w2tol) and 
            w2<(w+w2tol) and 
            w3>(wh-w3tol) and 
            w3<(wh+w3tol) and 
           (w1<(wnoe-w1tol) or w1>(wnoe+w1tol))
           ): noesy_peaks.append(line_noe)
      setattr(result, str(split_hmqc[0])+'_3dpeaks', noesy_peaks) 
  for linei in HMQC:
    HMQC_newfile.write(linei[:-1])
    linei_split=linei.split()
    #spaces=len('  '+str(linei_split[0]))
    #HMQC_newfile.write('  '+str(linei_split[0])+(14-spaces)*' '+str(linei_split[1])+'   '+str(linei_split[2])+'   '+'\t\t')
    try:noe_peaki_list=getattr(result, str(linei_split[0])+'_3dpeaks')
    except AttributeError: pass
    peak_geminal=''
    if noe_peaki_list:
      counter=0
      for noe_peaki in noe_peaki_list:
        noe_peaki_split=noe_peaki.split()
        for linef in HMQC:
          linef_split=linef.split()
          if ((float(noe_peaki_split[1])>(float(linef_split[1])-w1tol)) and
              (float(noe_peaki_split[1])<(float(linef_split[1])+w1tol))):
            try:noe_peakf_list=getattr(result, str(linef_split[0])+'_3dpeaks')
            except AttributeError: pass
            if noe_peakf_list:
              for noe_peakf in noe_peakf_list:
                noe_peakf_split=noe_peakf.split()
                if ((float(noe_peakf_split[1])>(float(linei_split[1])-w1tol)) and
                    (float(noe_peakf_split[1])<(float(linei_split[1])+w1tol))):
                  peak_geminal+=';'+str(linef_split[0])
                  CHH_geminal.write(str(linef_split[0][:-2])+'-'+str(linei_split[0])+'   '+str(linef_split[1])+
                                    '   '+str(linei_split[1])+'   '+str(linei_split[2])+'\n')
    HMQC_newfile.write('\t'+peak_geminal+'\n')
  HMQC_newfile.close()
  CHH_geminal.close()
print 'New HMQC peak list generated, please review it for accuracy.'