import types
import re
import numpy as np
import math
import datetime as time
import sys
import os
import multiprocessing as mp
import cPickle
import shutil
import resource
import psutil
from random import shuffle
from operator import itemgetter, attrgetter, methodcaller

class obj:pass
result=obj()    
# ---------------------------------------------------------------------
# Extract distances between of interest atoms
#
def distances(pdb,Cutoff,lowCut,flag):
  A,M,T,ICD1,ICG2,LCD1,LCD2,VCG1,VCG2=0,0,0,0,0,0,0,0,0
  for i in labeling.split(';'):
    if i[0]=='A':A=1
    elif i[0]=='M':M=1
    elif i[0]=='T':T=1
    elif i[0]=='I':
      for j in i.split(','):
        if j=='CD1':ICD1=1
        elif j=='CG2':ICG2=1
    elif i[0]=='L':
      for j in i.split(','):
        if j=='CD1':LCD1=1
        elif j=='CD2':LCD2=1
    elif i[0]=='V':
      for j in i.split(','):
        if j=='CG1':VCG1=1
        elif j=='CG2':VCG2=1   
  Atomlist=[]
  Met_alone=[]
  convert=[('ALA', 'A'),
           ('ARG', 'R'),
           ('ASP', 'D'),
           ('ASN', 'N'),
           ('CYS', 'C'),
           ('GLU', 'E'),
           ('GLN', 'Q'),
           ('GLY', 'G'),
           ('HIS', 'H'),
           ('ILE', 'I'),
           ('LEU', 'L'),
           ('LYS', 'K'),
           ('MET', 'M'),
           ('PHE', 'F'),
           ('PRO', 'P'),
           ('SER', 'S'),
           ('THR', 'T'),
           ('TRP', 'W'),
           ('TYR', 'Y'),
           ('VAL', 'V')]
  for line in pdb:
    split = line.split()
    if (split[0]=='ATOM' and (split[2]=='H' or 
                             (A==1 and split[3]=='ALA' and  split[2]=='CB') or 
                             (ICG2==1 and split[3]=='ILE' and split[2]=='CG2') or 
                             (ICD1==1 and split[3]=='ILE' and  split[2]=='CD1') or
                             ((LCD2==1 or LCD1==1) and split[3]=='LEU' and split[2]=='CD2') or
                             ((LCD2==1 or LCD1==1) and split[3]=='LEU' and split[2]=='CD1') or
                             (M==1 and split[3]=='MET' and  split[2]=='CE') or
                             (T==1 and split[3]=='THR' and  split[2]=='CG2') or  
                             ((VCG1==1 or VCG2==1) and split[3]=='VAL' and split[2]=='CG1') or 
                             ((VCG1==1 or VCG2==1) and split[3]=='VAL' and split[2]=='CG2')                     
      )): Atomlist.append(line)
  Methyl_list=[]
  distances_CHCH=[]
  if peak_geminal.sum()==0:flag_no_pairs=1
  else:flag_no_pairs=0
  Atomlist2=[]
  for i in range(len(Atomlist)):
    split1=Atomlist[i].split()
    if ((split1[3]=='LEU' and split1[2]=='CD1') or 
        (split1[3]=='VAL' and split1[2]=='CG1')
        ):
      split2=Atomlist[i+1].split()
      xc=round((float(split1[6])+float(split2[6]))/float(2),3)
      yc=round((float(split1[7])+float(split2[7]))/float(2),3)
      zc=round((float(split1[8])+float(split2[8]))/float(2),3)
      Atomlist2.append([split1[0],split1[1],split1[2][:-1],split1[3],split1[4],split1[5],xc,yc,zc])
    elif ((split1[3]=='LEU' and split1[2]=='CD2') or 
          (split1[3]=='VAL' and split1[2]=='CG2')
          ):pass
    else:Atomlist2.append(split1)
  
  for linei in Atomlist2:
    xi=float(linei[6])
    yi=float(linei[7])
    zi=float(linei[8])
    for aa in convert:
      if linei[3]==aa[0]: Res_name=aa[1]
    if linei[2]!='H':  
      Atom_name=Res_name+linei[5]+linei[2]
      if not Atom_name in Methyl_list:Methyl_list.append(Atom_name)     
    if (flag_geminal==1 and ((Res_name=='L' and linei[2]=='CD') or 
                             (Res_name=='V' and linei[2]=='CG'))
                             ):
      coef=1.2
      if flag=='run':distances_CHCH.append([str(Res_name+linei[5]+linei[2]),
                                            str(Res_name+linei[5]+linei[2]),coef])                               
    for linef in Atomlist2:
      xf=float(linef[6])
      yf=float(linef[7])
      zf=float(linef[8])
      if (xf!=xi or yf!=yi or zf!=zi):
        d=round(math.pow(math.pow(xf-xi,2)+math.pow(yf-yi,2)+math.pow(zf-zi,2),0.5),1) 
        if (linei[2]!='H' and linef[2]!='H'):   
            if (d<=lowCut and linei[5]!=linef[5]):coef=1
            elif d>lowCut and d<=Cutoff:coef=round((Cutoff-d)/float(Cutoff-lowCut),2)
            else:coef=0
            for line in convert:
              if linei[3]==line[0]:Res_name_i=line[1]
              if linef[3]==line[0]:Res_name_f=line[1]         
            if flag=='run':
              distances_CHCH.append([str(Res_name_i+linei[5]+linei[2]),
                                     str(Res_name_f+linef[5]+linef[2]),coef])
              #if coef!=0:methyl_file.write(str(distances_CHCH[-1])+'\n')
            elif flag=='histo':distances_CHCH.append([str(Res_name_i+linei[5]+linei[2]),
                                                      str(Res_name_f+linef[5]+linef[2]),d])
  empty=[]
  matrix_CHCH,matrix_geminal=matrix_it(distances_CHCH, Methyl_list, 'pdb',empty,empty,empty)
  matrix_CHCH_cluster=0

  return matrix_geminal, matrix_CHCH, Methyl_list, distances_CHCH
# ---------------------------------------------------------------------
# Extract 3d noesy peaks for each 2d peaks
#
def Extract_3DPeaks(object, peak_list_noesy, peak_list_2d_CH):    
  wtol_noeC = float(tolerance[0])
  wtol_refC = float(tolerance[-2])
  wtol_refH = float(tolerance[-1])
  if len(nuclei)==4:wtol_noeH = float(tolerance[1])
  Ref_nuclei=nuclei[-2]
  peak_list_2d_ref = peak_list_2d_CH
  DistNOE=[]
  tot_intensity=0
  Max_noe=0
  liste_NOE_intensity=[]
  liste_NOE=[]
  for line_noe in peak_list_noesy[2:]:
    split_noe=line_noe.split()   
    tot_intensity+=float(split_noe[-1])    
    liste_NOE_intensity.append(float(split_noe[-1]))    
    if float(split_noe[-1])>Max_noe:Max_noe=float(split_noe[-1])
  mean_intensity_per_strip=tot_intensity/float(len(peak_list_2d_CH))
  mean_intensity=tot_intensity/float(len(peak_list_noesy)-2)
  array_NOE=np.array(liste_NOE_intensity)
  NOE_std=np.std(array_NOE)
  NOE_per_Strip=round((len(peak_list_noesy)-2)/float(len(peak_list_2d_CH)),2)
  DistNOE.append(mean_intensity_per_strip)
  DistNOE.append(Max_noe)
  DistNOE.append(len(peak_list_noesy)-2)
  DistNOE.append(mean_intensity)
  DistNOE.append(NOE_per_Strip)
  file_log.write('NOESY data: '+str(NOESY_name)+'\n')
  file_log.write('Number of NOEs: '+str(len(peak_list_noesy)-2)+'\n')
  file_log.write('NOEs per strip: '+str(NOE_per_Strip)+'\n')  
  ####### Start sorting NOE according to 2D frequencies
  for line_hmqc in peak_list_2d_ref:
      split_hmqc=line_hmqc.split()
      w=float(split_hmqc[1])
      wh=float(split_hmqc[2])
      setattr(object, str(split_hmqc[0])+'_freq', (w,wh))
      noesy_peaks=[]
      for line_noe in peak_list_noesy[2:]:
        split_noe=line_noe.split()
        wnoeC=float(split_noe[1])	
        wrefC=float(split_noe[-3])
        wrefH=float(split_noe[-2])
        if len(nuclei)==4:
          wnoeH=float(split_noe[2])
          if (wrefC>(w-wtol_refC) and 
              wrefC<(w+wtol_refC) and 
              wrefH>(wh-wtol_refH) and 
              wrefH<(wh+wtol_refH) and 
             (wnoeC<(w-wtol_noeC) or wnoeC>(w+wtol_noeC) 
                                  or wnoeH<(wh-wtol_noeH) 
                                  or wnoeH>(wh+wtol_noeH))
             ):noesy_peaks.append(line_noe)        
        elif nuclei[0]=='1H':
          if (wrefC>(w-wtol_refC) and 
              wrefC<(w+wtol_refC) and 
              wrefH>(wh-wtol_refH) and 
              wrefH<(wh+wtol_refH) and 
             (wnoeC<(wh-wtol_noeC) or wnoeC>(wh+wtol_noeC))
             ):noesy_peaks.append(line_noe)        
        else:
          if (wrefC>(w-wtol_refC) and 
              wrefC<(w+wtol_refC) and 
              wrefH>(wh-wtol_refH) and 
              wrefH<(wh+wtol_refH) and 
             (wnoeC<(w-wtol_noeC) or wnoeC>(w+wtol_noeC))
             ):noesy_peaks.append(line_noe)
      setattr(object, str(split_hmqc[0])+'_3dpeaks', noesy_peaks)   
  return DistNOE
# ----------------------------------------------------------------------------------------
# Build up the matrix
#
def matrix_it(links, element_list, flag, factors,geminal_mark,overlap_topo):
  N=len(element_list)
  matrix=np.zeros((N,N))
  if flag=='peak':
    matrix_scoring=np.zeros((N,N))
    peak_geminal=np.zeros((N,N))  
    for entry in geminal_mark:peak_geminal[entry[0],entry[1]]=1
  if flag=='pdb':matrix_geminal=np.zeros((N,N))
  for line in links:
    i=element_list.index(line[0])
    j=element_list.index(line[1])
    #if flag=='peak':matrix[i,i]=1    
    if flag=='pdb':
      matrix[i,j]=line[2]
      if i==j:matrix_geminal[i,j]=1
    if flag=='peak':
      chi=float(line[6])
      if ((float(line[2])>0) or factors[2]>0):
        if float(line[2])>0.2:intensity_factor=1
        else:intensity_factor=0.5
        if (line[4]==0):matrix[i,j]=float(line[7])*intensity_factor*math.pow(math.pow(1+float(line[5]),2)*float(factors[0]/(20*chi*math.pow(float(line[3]),2))),0.5)
        elif (line[4]==1):matrix[i,j]=float(line[7])*intensity_factor*math.pow(math.pow(1+float(line[5]),2)*float(factors[1]/(20*chi*math.pow(float(line[3]),2))),0.5)
        elif (line[4]==2 or line[4]==3):matrix[i,j]=float(line[7])*math.pow(math.pow(1+float(line[5]),2)*float(factors[2]*intensity_factor/(20*chi*math.pow(float(line[3]),2))),0.5)
      else:matrix[i,j]=0
  if flag=='peak': 
    for a in range(matrix.size):
      i=int(a)/int(N)
      j=int(a)%int(N)
      if (matrix[i,j]<1 and matrix[j,i]<1):matrix[i,j],matrix[j,i]=0,0
      if (matrix[i,j]<1 and matrix[j,i]>2 and len(overlap_topo[j])==0):matrix[i,j]=1
  
  if flag=='pdb':return matrix,matrix_geminal
  if flag=='peak':
    matrix_ri=np.zeros((N,N))
    for line in links:
      i=element_list.index(line[0])
      j=element_list.index(line[1])
      chi=float(line[6])
      matrix_ri[i,j]=float(line[2])
      matrix_scoring[i,j]=float(line[7])*math.pow(math.pow(1+float(line[5]),2)/(20*chi*math.pow(float(line[3]),2)),0.5)
    matrix_scoring=matrix_scoring*matrix_ri
    return matrix,matrix_scoring,peak_geminal
# ---------------------------------------------------------------------------
#
def noe2matrix(result,factors,flag_cchlist):
    Noe_correlations_clustering_ALL=[]
    HMQC_peak_list=[]
    if flag_cchlist==0:noe_network=open('./'+directory.split('.')[0]+'/Peak_connectivity_'+str(NOESY_name), 'w')
    elif flag_cchlist==1:noe_network=open('./'+directory.split('.')[0]+'/Peak_connectivity_ISO_'+str(NOESY_name), 'w')
    overlap_list=[]
    overlap_topo=[]
    assignment_locked_peak=[]
    biblio_crosspeaks={}
    for i in range(len(HMQC)):
      checked_peak=HMQC[i]
      checked_peak_split=checked_peak.split()
      try:
        if checked_peak_split[3]=='-':assignment_locked_peak.append(i)
      except IndexError:pass
      overlaping_peaks=[]
      for j in range(len(HMQC)):
        peak_around=HMQC[j]
        peak_around_split=peak_around.split()
        if (float(peak_around_split[1])>(float(checked_peak_split[1])-0.1) and 
            float(peak_around_split[1])<(float(checked_peak_split[1])+0.1) and
            float(peak_around_split[2])>(float(checked_peak_split[2])-0.01) and 
            float(peak_around_split[2])<(float(checked_peak_split[2])+0.01) and
            checked_peak_split[0]!=peak_around_split[0]):
          overlaping_peaks.append(j)
          if not checked_peak_split[0] in overlap_list: overlap_list.append(checked_peak_split[0])
          if not peak_around_split[0] in overlap_list: overlap_list.append(peak_around_split[0])
      overlap_topo.append(overlaping_peaks)

    geminal_mark=[]
    for i in range(len(HMQC)):
      if (len(HMQC[i].split())==5 and len(overlap_topo[i])==0):
        geminal_i_list=HMQC[i].split()[4].split(';')
        for geminals in geminal_i_list[1:]:
          for j in range(len(HMQC)):
            if (HMQC[j].split()[0]==geminals and len(overlap_topo[j])==0):geminal_mark.append((i,j,len(geminal_i_list[1:])))
    

    for linei in HMQC:
       Noe_correlations_clustering=[]
       Noe_correlations_scoring=[]

       height_i2mean=1
       score_lowNOE=0
       score_tot=0
       linei_split=linei.split()
       HMQC_peak_list.append(linei_split[0])
       
       try:
           noe_peaki_list=getattr(result, str(linei_split[0])+'_3dpeaks')
           total_height_i=0
           imax=0
           for peak in noe_peaki_list:
             peak_split=peak.split()
             total_height_i+=float(peak_split[-1])
             if float(peak_split[4])>imax:imax=float(peak_split[-1])
           if total_height_i>0:rimax=float(imax/total_height_i)
           height_i2mean=round(float(total_height_i)/DistNOE[0], 2)
       except AttributeError: pass
       if not str(linei_split[0]) in biblio_crosspeaks.keys():biblio_crosspeaks[str(linei_split[0])]={}     
       if noe_peaki_list:
         counter=0     
         for noe_peaki in noe_peaki_list:
             noe_peaki_split=noe_peaki.split()          
             ri=round(float(noe_peaki_split[-1])/total_height_i, 2)
             ri2max=round(float(noe_peaki_split[-1])/DistNOE[1],3)    
             peak_neighbors_clustering=[]
             peak_neighbors_scoring=[]
             shared_peak_neighbors=0
             for linef in HMQC:
                 linef_split=linef.split()
                 
                 ##### FOR 3D
                 if (len(nuclei)==3 and nuclei[0]=='13C'):
                   if ((float(noe_peaki_split[1])>(float(linef_split[1])-float(tolerance[0]))) and 
                       (float(noe_peaki_split[1])<(float(linef_split[1])+float(tolerance[0])))):
                     shared_peak_neighbors=0
                     sym=sym_coef
                     ppm_matching=[0.035,0.035]
                     ppm_matching[0]=float(noe_peaki_split[1])-float(linef_split[1])
                     try:
                       noe_peakf_list=getattr(result, str(linef_split[0])+'_3dpeaks')
                     except AttributeError: pass
                     if noe_peakf_list:
                       for noe_peakf in noe_peakf_list:
                         noe_peakf_split=noe_peakf.split()
                         for noe_peakii in noe_peaki_list:
                           noe_peakii_split=noe_peakii.split()
                           if ((float(noe_peakii_split[1])>(float(noe_peakf_split[1])-float(tolerance[0]))) and 
                               (float(noe_peakii_split[1])<(float(noe_peakf_split[1])+float(tolerance[0])))
                                ):shared_peak_neighbors+=1
                         if ((float(noe_peakf_split[1])>(float(linei_split[1])-float(tolerance[0]))) and 
                             (float(noe_peakf_split[1])<(float(linei_split[1])+float(tolerance[0]))) and
                             (not str(linef_split[0]) in peak_neighbors_clustering)):
                             ppm_matching[1]=float(noe_peakf_split[1])-float(linei_split[1])
                             sym=1
                     peak_neighbors_clustering.append([str(linef_split[0]),shared_peak_neighbors,ppm_matching,sym])
                     if sym==1:peak_neighbors_scoring.append((str(linef_split[0]),shared_peak_neighbors,ppm_matching)) # for noe matrix used for scoring                        

                 elif len(nuclei)==3 and nuclei[0]=='1H':
                   if ((float(noe_peaki_split[1])>(float(linef_split[2])-float(tolerance[0]))) and 
                       (float(noe_peaki_split[1])<(float(linef_split[2])+float(tolerance[0])))):
                     shared_peak_neighbors=0
                     sym=sym_coef
                     ppm_matching=[0.035,0.035]
                     ppm_matching[0]=float(noe_peaki_split[1])-float(linef_split[2])
                     try:
                       noe_peakf_list=getattr(result, str(linef_split[0])+'_3dpeaks')
                     except AttributeError: pass
                     if noe_peakf_list:
                       for noe_peakf in noe_peakf_list:
                         noe_peakf_split=noe_peakf.split()
                         for noe_peakii in noe_peaki_list:
                           noe_peakii_split=noe_peakii.split()
                           if ((float(noe_peakii_split[1])>(float(noe_peakf_split[2])-float(tolerance[0]))) and 
                               (float(noe_peakii_split[1])<(float(noe_peakf_split[2])+float(tolerance[0])))
                                ):shared_peak_neighbors+=1
                         if ((float(noe_peakf_split[1])>(float(linei_split[2])-float(tolerance[0]))) and 
                             (float(noe_peakf_split[1])<(float(linei_split[2])+float(tolerance[0]))) and
                             (not str(linef_split[0]) in peak_neighbors_clustering)):
                             ppm_matching[1]=float(noe_peakf_split[1])-float(linei_split[2])
                             sym=1
                     peak_neighbors_clustering.append([str(linef_split[0]),shared_peak_neighbors,ppm_matching,sym])
                     if sym==1:peak_neighbors_scoring.append((str(linef_split[0]),shared_peak_neighbors,ppm_matching)) # for noe matrix used for scoring                        
                 
                 
                 ### FOR 4D
                 else:
                   if ((float(noe_peaki_split[1])>(float(linef_split[1])-float(tolerance[0]))) and 
                       (float(noe_peaki_split[1])<(float(linef_split[1])+float(tolerance[0]))) and
                       (float(noe_peaki_split[2])>(float(linef_split[2])-float(tolerance[1]))) and 
                       (float(noe_peaki_split[2])<(float(linef_split[2])+float(tolerance[1])))                       
                       ):
                     shared_peak_neighbors=0
                     sym=sym_coef
                     ppm_matching=[0.035,0.035]
                     ppm_matching[0]=float(noe_peaki_split[1])-float(linef_split[1])
                     try:
                       noe_peakf_list=getattr(result, str(linef_split[0])+'_3dpeaks')
                     except AttributeError: pass
                     if noe_peakf_list:
                       for noe_peakf in noe_peakf_list:
                         noe_peakf_split=noe_peakf.split()
                         for noe_peakii in noe_peaki_list:
                           noe_peakii_split=noe_peakii.split()
                           if ((float(noe_peakii_split[1])>(float(noe_peakf_split[1])-float(tolerance[0]))) and 
                               (float(noe_peakii_split[1])<(float(noe_peakf_split[1])+float(tolerance[0]))) and
                               (float(noe_peakii_split[2])>(float(noe_peakf_split[2])-float(tolerance[1]))) and 
                               (float(noe_peakii_split[2])<(float(noe_peakf_split[2])+float(tolerance[1])))
                                ):shared_peak_neighbors+=1
                         if ((float(noe_peakf_split[1])>(float(linei_split[1])-float(tolerance[0]))) and 
                             (float(noe_peakf_split[1])<(float(linei_split[1])+float(tolerance[0]))) and
                             (float(noe_peakf_split[2])>(float(linei_split[2])-float(tolerance[1]))) and 
                             (float(noe_peakf_split[2])<(float(linei_split[2])+float(tolerance[1]))) and
                             (not str(linef_split[0]) in peak_neighbors_clustering)
                             ):
                             ppm_matching[1]=float(noe_peakf_split[1])-float(linei_split[1])
                             sym=1
                     peak_neighbors_clustering.append([str(linef_split[0]),shared_peak_neighbors,ppm_matching,sym])
                     if sym==1:peak_neighbors_scoring.append((str(linef_split[0]),shared_peak_neighbors,ppm_matching)) # for noe matrix used for scoring    
             n=0
             m=len(peak_neighbors_clustering)
             for i in range(len(peak_neighbors_clustering)):
               if peak_neighbors_clustering[i][3]==1:n+=1
             if len(peak_neighbors_clustering)>0:
               biblio_crosspeaks[str(linei_split[0])][str(counter)]=([peak_neighbors_clustering[i][0] for i in range(len(peak_neighbors_clustering))],ri)
               counter+=1
             for i in range(len(peak_neighbors_clustering)):
               if (not peak_neighbors_clustering[i][0] in overlap_list and not linei_split[0] in overlap_list):a=0
               elif (not peak_neighbors_clustering[i][0] in overlap_list and linei_split[0] in overlap_list):a=1
               elif (peak_neighbors_clustering[i][0] in overlap_list and not linei_split[0] in overlap_list):a=2
               elif (peak_neighbors_clustering[i][0] in overlap_list and linei_split[0] in overlap_list):a=3
               d1=float(peak_neighbors_clustering[i][2][0])
               d2=float(peak_neighbors_clustering[i][2][1])
               chi=round(math.pow(math.pow(d1,2)+math.pow(d2,2),0.5),4)
               if chi<0.05:chi=0.05
               if peak_neighbors_clustering[i][3]==1:N=n
               else:N=m
               
               sym=peak_neighbors_clustering[i][3]
               ################## CONNECTION UPGRADING ###################
               if (sym==0 and ((n==0 and m>1 and int(peak_neighbors_clustering[i][1])>=2) or 
                               (n==0 and m==1 and int(peak_neighbors_clustering[i][1])>=1) or
                               (n>=1 and int(peak_neighbors_clustering[i][1])>=4))
                   ):
                 sym=0.5
                 score=round(float(sym)*math.pow(math.pow(1+float(peak_neighbors_clustering[i][1]),2)*float(factors[0]/(20*chi*math.pow(float(N),2))),0.5),3)
               ###########################################################
               Noe_correlations_clustering.append([linei_split[0],peak_neighbors_clustering[i][0],
                                                  ri,N,a,peak_neighbors_clustering[i][1],chi,sym])
               if flag_cchlist==0 or flag_cchlist==1:
                 score=round(float(sym)*math.pow(math.pow(1+float(peak_neighbors_clustering[i][1]),2)*float(factors[0]/(20*chi*math.pow(float(N),2))),0.5),3)
                 noe_network.write(str(Noe_correlations_clustering[-1])+'  '+str(score)+'\n')
                                       
       for element in Noe_correlations_clustering:Noe_correlations_clustering_ALL.append(element)
                        
    noe_peak_matrix_clustering,noe_peak_matrix_scoring,peak_geminal=matrix_it(Noe_correlations_clustering_ALL, HMQC_peak_list, 'peak',factors,geminal_mark,overlap_topo)    
    if flag_cchlist==0 or flag_cchlist==1:noe_network.close()
    return noe_peak_matrix_scoring, noe_peak_matrix_clustering, HMQC_peak_list,overlap_topo,biblio_crosspeaks,assignment_locked_peak,peak_geminal

# ---------------------------------------------------------------------------
#
def build_assignment(index):
              highest_score_inloop=highest_score
              flag_ram=0
              assignment_archive=[]
              assignment_collection={}
              index_files=0
              Time=time.datetime.now()
              report_follow=open('./'+directory.split('.')[0]+'/'+str(path_name)+'/'+str(logFile_name), 'a')
              report_follow.write(str(Time).split('.')[0]+'\n')
              report_follow.close()              
              file_log=open('./'+directory.split('.')[0]+'/log', 'a')
              file_size=int(1000+float(file_size_max*size/float(5+size)))              
              file_log.write('file_size: '+str(file_size)+'\n')
              file_log.close()
              for iii in index:             
                free_RAM=round(float(psutil.virtual_memory().available)/float(1024*1024*1024),1)              
                file_log=open('./'+directory.split('.')[0]+'/log', 'a')
                file_log.write('free_RAM: '+str(free_RAM)+' GB'+'\n')          
                file_log.write(str(index.index(iii)+1)+'/'+str(len(index))+'\n')
                file_log.close()            
                file_name=list_of_files[iii]
                file=open('./'+directory.split('.')[0]+'/run/temp/'+str(file_name), 'r')
                assignment_archive_FINAL=cPickle.load(file)
                file.close()         
                if len(assignment_archive)==0:
                  index_peaks=assignment_archive_FINAL[0]+involved_peak
                  matrix_noe=peak_matrix_Scoring[:,index_peaks][index_peaks,:]
                  assignment_archive.append(index_peaks)
                  matrix_noe_geminal=peak_geminal[:,index_peaks][index_peaks,:]

                for ii in range(len(assignment_archive_FINAL[1:])):
                  for element in possible_peak_assignments:
                    index_methyls=assignment_archive_FINAL[ii+1][0]+element                                    
                    flag_Stop=0
                    for i in index_methyls:
                      if (index_methyls.count(i)>2 or
                          (index_methyls.count(i)==2 and ((metrics[2][i][0]!='L' and metrics[2][i][0]!='V') or 
                                                          flag_geminal==2))
                         ):
                         flag_Stop=1
                         break         
                    if flag_Stop==1:continue         
                    
                    matrix_assignment=new_metrics[0][:,index_methyls][index_methyls,:]

                    #####Geminal scaling factor####################
                    if (flag_geminal==1 and matrix_noe_geminal.sum()>0):
                      matrix_assignment_geminal=new_metrics[1][:,index_methyls][index_methyls,:]
                      testing_geminal=matrix_noe_geminal*matrix_assignment_geminal
                      if (geminal_hold==1 and np.any(testing_geminal-matrix_noe_geminal)==True):continue
                      elif geminal_hold==0:matrix_assignment=matrix_assignment+(testing_geminal/float(5))
                      else:pass
                    ####Geminal scaling factor

                    testing=matrix_assignment*matrix_noe
                    tot=round(testing.sum(),2)

                    if tot>highest_score_inloop:highest_score_inloop=tot                
                    cutoff=highest_score_inloop-FILTER

                    if tot>=cutoff:assignment_archive.append([index_methyls,tot])
                    else:pass
                      
                                       
                  free_RAM=round(float(psutil.virtual_memory().available)/float(1024*1024*1024),1)              
                  if ((free_RAM<=0.7*free_RAM_start and len(assignment_archive)>file_size) or 
                      (free_RAM<=0.5*free_RAM_start and len(assignment_archive)>10000)): 
                    cutoff=highest_score_inloop-FILTER
                    tot=len(assignment_archive)-1
                    deleted=0
                    assignment_archive_sorted=[assignment_archive[0]]
                    for i in range(len(assignment_archive)-1):
                      if (assignment_archive[i+1][1]>=cutoff):
                        assignment_archive_sorted.append(assignment_archive[i+1])                 
                        for j in range(len(assignment_archive[i+1][0])):
                          peak=HMQC_peak_list[int(assignment_archive[0][j])]
                          methyl=metrics[2][int(assignment_archive[i+1][0][j])]
                          if not peak in assignment_collection.keys():
                            assignment_collection[peak]={}
                            assignment_collection[peak][methyl]=assignment_archive[i+1][1]
                          else:
                            if not methyl in assignment_collection[peak].keys():
                              assignment_collection[peak][methyl]=assignment_archive[i+1][1]
                            elif (methyl in assignment_collection[peak].keys() and 
                              assignment_archive[i+1][1]>assignment_collection[peak][methyl]):
                              assignment_collection[peak][methyl]=assignment_archive[i+1][1]
                      else:deleted=deleted+1
                    report_follow=open('./'+directory.split('.')[0]+'/'+str(path_name)+'/'+str(logFile_name), 'a')
                    report_follow.write('----total#'+str(tot-deleted)+'\n')
                    report_follow.close()
                    assignment_archive=list(assignment_archive_sorted)
                    assignment_archive_sorted=[]
                    number_of_file=int((len(assignment_archive)-1))/int(file_size)
                    number_of_file_rest=(len(assignment_archive)-1)%int(file_size)          
                    for i in range(number_of_file):
                        assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[i*file_size+1:(i+1)*file_size+1]
                        file=open('./'+directory.split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                        cPickle.dump(assignment_archive_to_save,file,-1)
                        file.close()
                        index_files+=1
                    if number_of_file_rest!=0:
                        assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[number_of_file*file_size+1:]
                        file=open('./'+directory.split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                        cPickle.dump(assignment_archive_to_save,file,-1)
                        file.close()
                        index_files+=1
                    assignment_archive=[assignment_archive[0]] 
                    
              if len(assignment_archive)>1:                 
                  cutoff=highest_score_inloop-FILTER
                  tot=len(assignment_archive)-1
                  deleted=0
                  assignment_archive_sorted=[assignment_archive[0]]
                  for i in range(len(assignment_archive)-1):
                    if (assignment_archive[i+1][1]>=cutoff):
                      assignment_archive_sorted.append(assignment_archive[i+1])                 
                      for j in range(len(assignment_archive[i+1][0])):
                        peak=HMQC_peak_list[int(assignment_archive[0][j])]
                        methyl=metrics[2][int(assignment_archive[i+1][0][j])]
                        if not peak in assignment_collection.keys():
                          assignment_collection[peak]={}
                          assignment_collection[peak][methyl]=assignment_archive[i+1][1]
                        else:
                          if not methyl in assignment_collection[peak].keys():
                            assignment_collection[peak][methyl]=assignment_archive[i+1][1]
                          elif (methyl in assignment_collection[peak].keys() and 
                            assignment_archive[i+1][1]>assignment_collection[peak][methyl]):
                            assignment_collection[peak][methyl]=assignment_archive[i+1][1]
                    else:deleted=deleted+1
                  report_follow=open('./'+directory.split('.')[0]+'/'+str(path_name)+'/'+str(logFile_name), 'a')
                  report_follow.write('----total#'+str(tot-deleted)+'\n')
                  report_follow.close()
                  assignment_archive=list(assignment_archive_sorted)
                  assignment_archive_sorted=[]                 
                  if (len(assignment_archive)-1)>file_size:
                    number_of_file=(len(assignment_archive)-1)/int(file_size)
                    number_of_file_rest=(len(assignment_archive)-1)%int(file_size)
                    for i in range(number_of_file):
                      assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[i*file_size+1:(i+1)*file_size+1]
                      file=open('./'+directory.split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                      cPickle.dump(assignment_archive_to_save,file,-1)
                      file.close()
                      index_files+=1
                    if number_of_file_rest!=0:
                      assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[number_of_file*file_size+1:]
                      file=open('./'+directory.split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                      cPickle.dump(assignment_archive_to_save,file,-1)
                      file.close()
                      index_files+=1
                    assignment_archive=[]
                  else:
                    file=open('./'+directory.split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                    cPickle.dump(assignment_archive,file,-1)
                    file.close()                    
                    assignment_archive=[]
              max=highest_score_inloop
              return assignment_collection,max

def build_assignment_peak(index):
        count=0
        global P
        global q
        global mean_MetType
        global total_Assignments
        global FILTER_scale
        global geminal_hold
        assignment_archive=[]        
        for index_number in index:
          assignment_old=assignment_archive_old[index_number]     
          if total_Assignments==[]:total_Assignments=[[],[]]
          index_matrix_noe=total_Assignments[0]+assignment_old[0]
          matrix_noe=peak_matrix_Scoring[:,index_matrix_noe][index_matrix_noe,:]
          matrix_noe_geminal=peak_geminal[:,index_matrix_noe][index_matrix_noe,:]
          
          for assignment in total_Assignments[1:]:          
            index_matrix_methyls=assignment+assignment_old[1]
            flag_Stop=0
            for i in index_matrix_methyls:
              if (index_matrix_methyls.count(i)>2 or 
                  (index_matrix_methyls.count(i)==2 and 
                   ((metrics[2][i][0]!='L' and metrics[2][i][0]!='V') or 
                    flag_geminal==2))
                  ):
                  flag_Stop=1
                  break   
            if flag_Stop==1:continue
            matrix_assignment=new_metrics[0][:,index_matrix_methyls][index_matrix_methyls,:]
            
            #####IN: Test allowed assignment for all peaks #######
            #if P<=P10:
            if P<=P_high:
              try:
                flag_stop=0
                for i in range(len(index_matrix_noe)):
                  peak=HMQC_peak_list[int(index_matrix_noe[i])]
                  methyl=metrics[2][int(index_matrix_methyls[i])]
                  if (peak in archive_assignment_result.keys() and
                     len(overlap_topo[int(index_matrix_noe[i])])==0):
                    if not methyl in archive_assignment_result[peak].keys():
                      flag_stop=1
                      break
                    else:pass
                  else:pass
                if flag_stop==1:continue
              except NameError:pass
            else:pass
            #####OUT: Test allowed assignment for all peaks #######

            #####Geminal scaling factor
            if (flag_geminal==1 and matrix_noe_geminal.sum()>0):
              matrix_assignment_geminal=new_metrics[1][:,index_matrix_methyls][index_matrix_methyls,:]
              testing_geminal=matrix_noe_geminal*matrix_assignment_geminal
              if (geminal_hold==1 and np.any(testing_geminal-matrix_noe_geminal)==True):continue
              elif geminal_hold==0:matrix_assignment=matrix_assignment+(testing_geminal/float(5))
              else:pass            
            #####Geminal scaling factor
            testing=matrix_assignment*matrix_noe
            tot=testing.sum()
            assignment_archive.append([index_matrix_noe,index_matrix_methyls,tot])
            count+=1
                           
        if len(assignment_archive)>1:          
          assignment_archive_brut=list(assignment_archive)
          scorelist=np.zeros((1,len(assignment_archive_brut)))
          for i in range(len(assignment_archive_brut)):scorelist[0,i]=assignment_archive_brut[i][2]
          max=scorelist[0,np.argmax(scorelist)]
          if max<highest_score_small:max=highest_score_small
          number_peaks=len(index_matrix_noe)
          if number_peaks<=2:q=0           
          elif number_peaks==3:q=0.01
          elif number_peaks==4:q=0.7
          elif number_peaks==5:q=0.8
          elif number_peaks>=6:q=0.9
          if flag_fin!=1:q*=0.9
          if mean_MetType>1:q*=0.9
          if FILTER_scale>1:scale_cluster=FILTER_scale
          else:scale_cluster=1
          cutoff=highest_score_small-highest_score_small*(1-q)*scale_cluster
          
          tot=len(assignment_archive_brut)
          deleted=0
          assignment_archive=[]
          for i in range(len(assignment_archive_brut)):
                if (assignment_archive_brut[i][2]>=cutoff):
                  assignment_archive.append(assignment_archive_brut[i])
                else:deleted=deleted+1
          return assignment_archive,max
        else: 
          max=0
          return assignment_archive,max

# ---------------------------------------------------------------------------
# Build up OUTPUT: hmqc+cch peaklists, histogram, checking properties
#
def outputing(highest_score,peak_matrix_Scoring):
    global P   
    file_log=open('./'+directory.split('.')[0]+'/log', 'a')
    file_log.write('----> Building up OUTPUT files...'+'\n')
    file_log.close()

    if P=='iso':
      hmqc_result_list=open('./'+directory+'/Output/hmqc_iso.list', 'w')
      cch_result_list=open('./'+directory+'/Output/cch_iso.list', 'w')
    else:
      hmqc_result_list=open('./'+directory+'/Output/hmqc.list', 'w')
      cch_result_list=open('./'+directory+'/Output/cch.list', 'w')
          
    hmqc_result_list.write('Assignment\tw1\tw2\tNote\n\n')
    cch_result_list.write('Assignment\tw1\tw2\tw3\tNote\n\n')

############## TO IMPROVE ###############
    peaks=[]
    methyls=[]
    list_peak=archive_assignment_result.keys()
    list_peak.sort()   
    for peak_name in list_peak:
      peaks.append(HMQC_peak_list.index(peak_name))      
      score=[]
      for i in range(len(archive_assignment_result[peak_name].keys())):
        methyl_name=archive_assignment_result[peak_name].keys()[i]
        score.append((archive_assignment_result[peak_name][methyl_name],methyl_name))
      score.sort()
      score=score[::-1]
      score_max=float(score[0][0])
      flag=0
      for i in range(len(score)):
        if i==0:continue
        elif float(score[i][0])==score_max and i>0:
          flag=1
          break
        else:break        
      if flag==0:methyls.append(metrics[2].index(score[0][1]))
      else:
        select_methyl=''
        for i in range(len(score)):
          if (score[i][0]==score_max and 
              (methyls.count(metrics[2].index(score[i][1]))==0 or
               (methyls.count(metrics[2].index(score[i][1]))==1 and flag_geminal!=2 
                                                               and (score[i][1][0]=='V' or score[i][1][0]=='L')))
             ):
             select_methyl=score[i][1]
             break
        if select_methyl=='':
          print 'error',peak_name
          methyls.append(metrics[2].index(score[0][1]))
        else:methyls.append(metrics[2].index(select_methyl))
############## TO IMPROVE ###############    
    
    matrix_assignment=new_metrics[0][:,methyls][methyls,:]
    matrix_assignment=np.where(matrix_assignment>1.0,1.0,matrix_assignment)
    factors=(1,1,0.25)
    sym_coef=0.1
    peak_matrix_Scoring2,peak_matrix_clustering2,glop,glip,glep,glyp,glap=noe2matrix(result,factors,3)
    matrix_peaks_CS=peak_matrix_clustering2[peaks,:][:,peaks] ###peak_matrix_clustering2 instead of peak_matrix_clustering if factor is changed
    matrix_peaks_Cluster=peak_matrix_clustering[peaks,:][:,peaks]
    #matrix_peaks=peak_matrix_Scoring2[peaks,:][:,peaks]
    matrix_peaks=peak_matrix_Scoring[peaks,:][:,peaks]
    matrix_score=matrix_assignment*matrix_peaks
    
    #matrix_peaks_NOEdensity=np.dot(matrix_peaks_Cluster,matrix_peaks_Cluster)
    
    correct=0
    total=0
    metrics2histo=distances(pdb,100,100,'histo')
    distance2histo,glop=matrix_it(metrics2histo[3], metrics[2], 'pdb',empty,empty,empty)
    matrix_methyls=distance2histo[methyls,:][:,methyls]    
    
    hmqc_new_list=[]
    
    for i in range(len(peaks)):
          peak=HMQC_peak_list[int(peaks[i])]
          methyl=metrics[2][int(methyls[i])]  
          try:
            if re.search('\w(\d+)\w',str(methyl)).group(1)==re.search('(\d+)\w', peak).group(1):correct+=1
          except AttributeError:pass
          total+=1
          freq=getattr(result, str(peak)+'_freq')
          wC=freq[0]
          wH=freq[1]
          if methyl[0]=='I':Name=str(methyl)+'-HD1'
          elif methyl[0]=='M':Name=str(methyl)+'-HE'
          elif methyl[0]=='L' and flag_L_nuclei=='CD1':Name=str(methyl)+'1-HD1'
          elif methyl[0]=='L' and flag_L_nuclei=='CD2':Name=str(methyl)+'2-HD2'
          elif (methyl[0]=='L' and flag_L_nuclei=='nd' and methyl in hmqc_new_list):Name=str(methyl)+'2-HD2'
          elif (methyl[0]=='L' and flag_L_nuclei=='nd' and not methyl in hmqc_new_list):Name=str(methyl)+'1-HD1'
          elif methyl[0]=='V' and flag_V_nuclei=='CG1':Name=str(methyl)+'1-HG1'
          elif methyl[0]=='V' and flag_V_nuclei=='CG2':Name=str(methyl)+'2-HG2'
          elif (methyl[0]=='V' and flag_V_nuclei=='nd' and methyl in hmqc_new_list):Name=str(methyl)+'2-HG2'
          elif (methyl[0]=='V' and flag_V_nuclei=='nd' and not methyl in hmqc_new_list):Name=str(methyl)+'1-HG1'          
          elif methyl[0]=='A':Name=str(methyl)+'-HB'
          elif methyl[0]=='T':Name=str(methyl)+'-HG2'         
          hmqc_new_list.append(methyl)
          hmqc_result_list.write(Name+'\t'+str(wC)+'\t'+str(wH))
          already_assigned=[]
          completeness=np.zeros((1,len(biblio_crosspeaks[peak].keys())))
          for j in range(matrix_score[i,:].size):
                if matrix_score[i,j]>0:
                  peakNOE=HMQC_peak_list[int(peaks[j])]
                  methylNOE=metrics[2][int(methyls[j])]
                  wCnoe=getattr(result, str(peakNOE)+'_freq')[0]
                  if matrix_peaks_Cluster[i,j]>0:a=1
                  else:a=0
                  cch_result_list.write(str(methylNOE)+'-'+Name+'\t'+str(wCnoe)+'\t'+str(wC)+'\t'+str(wH))
                  cch_result_list.write('  '+str(round(matrix_peaks_CS[i,j],3))+'  d='+str(matrix_methyls[i,j])
                                        +'  '+str(a)+'\n')               
                  for NOE_index in biblio_crosspeaks[peak].keys():
                    for donor in biblio_crosspeaks[peak][NOE_index][0]:
                      if (donor==peakNOE and not k in already_assigned):already_assigned.append(k)                    
                  for NOE_index in biblio_crosspeaks[peak].keys():
                    if (peakNOE in biblio_crosspeaks[peak][NOE_index][0] and
                        completeness[0,int(NOE_index)]==0):
                      completeness[0,int(NOE_index)]=1
                      break
          if completeness.size>0:NOE_completeness=round(completeness.sum()/float(completeness.size),2)
          else: NOE_completeness='no NOE'
          
          Methyl_type=HMQC[int(peaks[i])].split()[3]        
          total_Density=round(matrix2[i,:].sum(),1)         
          total_Score=round(matrix_peaks[i,:].sum(),1)          
          Satisfied_Score=round(matrix_score[i,:].sum(),1)
          
          if float(total_Score)!=0:Score_completeness=round((float(Satisfied_Score)/float(total_Score)),1)
          #else:print 'total_Score=0 ??!',peak,methyl
          else:Score_completeness='nd'
          
          hmqc_result_list.write('  '+str(Methyl_type)+
                                 #'  '+str(round(matrix_peaks_CS[i,:].sum(),2))+
                                 '  '+str(total_Density)+
                                 '  '+str(NOE_completeness)+
                                 '  '+str(Score_completeness)+
                                 '  '+str(archive_assignment_result[peak]))
          hmqc_result_list.write('\n')    
                           
    cch_result_list.close()
    already_used=methyls[:]
    for i in range(len(HMQC_peak_list)):
      if not i in peaks:
        peak=HMQC_peak_list[i]
        possible_assignment=[]
        allowed_assignment=[]
        code_methyl_type=HMQC[i].split()[3]
        if 'A' in code_methyl_type:possible_assignment+=Ala_methyls
        if 'I' in code_methyl_type:possible_assignment+=Ile_methyls
        if 'L' in code_methyl_type:possible_assignment+=Leu_methyls
        if 'V' in code_methyl_type:possible_assignment+=Val_methyls
        if 'M' in code_methyl_type:possible_assignment+=Met_methyls
        if 'T' in code_methyl_type:possible_assignment+=Thr_methyls
        for methyl_index in possible_assignment:
          if (methyls.count(methyl_index)==0 or 
             (methyls.count(methyl_index)==1 and
              flag_geminal!=2 and
              (methyl_index in Leu_methyls or methyl_index in Val_methyls))
             ):allowed_assignment.append(metrics[2][methyl_index])
        freq=getattr(result, str(peak)+'_freq')
        wC=freq[0]
        wH=freq[1]
        Name2=HMQC[i].split()[0]
        Methyl_type=HMQC[i].split()[3]
        hmqc_result_list.write(Name2+'\t'+str(wC)+'\t'+str(wH))
        hmqc_result_list.write('   '+str(Methyl_type)+'\tNotAss'+str(allowed_assignment)+'\n')        
    hmqc_result_list.close()
    
    if P=='iso':
      file=open('./'+directory+'/Output/cch_iso.list', 'r')
      noelist=file.readlines()
      file.close()
      file2=open('./'+directory+'/Output/hmqc_iso.list', 'r')
      hmqclist=file2.readlines()
      file2.close()
      results=open('./'+directory+'/Output/mapping_iso.pml','w')
    else:
      file=open('./'+directory+'/Output/cch.list', 'r')
      noelist=file.readlines()
      file.close()
      file2=open('./'+directory+'/Output/hmqc.list', 'r')
      hmqclist=file2.readlines()
      file2.close()
      results=open('./'+directory+'/Output/mapping.pml','w')
                
    results.write('hide everything\n'+
                  'show ribbon\n'+
                  'set ribbon_width,5\n'+
                  'bg_color white\n'+
                  'color gray30\n'+
                  'set dash_width,4\n'+
                  'opaque_background set to off\n'+
                  'create Me, ')
    for i in labeling.split(';'):
            if i[0]=='A':results.write('resn ala and name cb')
            elif i[0]=='M':results.write('resn met and name ce')
            elif i[0]=='T':results.write('resn thr and name CG2')
            elif i[0]=='I':
              for j in i.split(',')[1:]:
                if j=='CD1':results.write('resn ile and name CD1')
                elif j=='CG2':results.write('resn ile and name CG2')
                #if j==i.split(',')[-1]:break
                results.write(' or ')
              continue
            elif i[0]=='L':
              for j in i.split(',')[1:]:
                if j=='CD1':results.write('resn leu and name CD1')
                elif j=='CD2':results.write('resn leu and name CD2')
                #if j==i.split(',')[-1]:break
                results.write(' or ')
              continue
            elif i[0]=='V':
              for j in i.split(',')[1:]:
                if j=='CG1':results.write('resn val and name cg1')
                elif j=='CG2':results.write('resn val and name cg2')
                #if j==i.split(',')[-1]:break
                results.write(' or ')
              continue
            if i==labeling.split(';')[-1]:break
            results.write(' or ')
    results.write('\nshow nb_sphere, object Me\ncolor gray70, object Me\n')
    M=3
    for line in noelist[2:]:
            split=line.split()
            if float(split[6])==0:continue
            else:pass
            split0=split[0].split('-')
            res1i=int(re.search('[A-Z]([0-9]+)[A-Z]', split0[0]).group(1))
            name1i=re.search('[0-9]+(\w+)', split0[0]).group(1)
            res2i=int(re.search('[A-Z]([0-9]+)[A-Z]', split0[1]).group(1))
            name2i=re.search('[0-9]+(\w+)', split0[1]).group(1)
            sm=float(split[4])
            for line_2 in noelist[2:]:
              split=line_2.split()
              split0f=split[0].split('-')
              res1f=int(re.search('[A-Z]([0-9]+)[A-Z]', split0f[0]).group(1))
              name1f=re.search('[0-9]+(\w+)', split0f[0]).group(1)
              res2f=int(re.search('[A-Z]([0-9]+)[A-Z]', split0f[1]).group(1))
              name2f=re.search('[0-9]+(\w+)', split0f[1]).group(1)  
              if (res1f==res2i and name1f==name2i and res1i==res2f and name1i==name2f):
                sf=float(split[4])
                sm=(sm+sf)/float(2)
                break	
            M4=float(M)
            M3=float(3*M/float(4))
            M2=float(2*M/float(4))
            M1=float(M/float(4))
            if (sm>M4):r,g,b=0,0,1
            elif (sm>M3):r,g,b=0,round(1-(sm-M3)/(M4-M3), 2),1
            elif (sm>M2):r,g,b=0,1,round((sm-M2)/(M3-M2), 2)
            elif (sm>M1):r,g,b=round(1-(sm-M1)/(M2-M1), 2),1,0
            else:r,g,b=1,round(sm/M1, 2),0
            if (split0[0][0]=='L' or split0[0][0]=='V'):star1='*'
            else:star1=''
            if (split0[1][0]=='L' or split0[1][0]=='V'):star2='*'
            else:star2=''
            #star1,star2='',''
            results.write('distance '+str(res1i)+'-'+str(name1i)+'_'+str(res2i)+'-'+str(name2i)+', object Me and resi '+str(res1i)+' and name '+str(name1i)+star1+', object Me and resi '+str(res2i)+' and name '+str(name2i)+star2+'\n')
            results.write('set_color rgb'+str(res1i)+'-'+str(name1i)+'_'+str(res2i)+'-'+str(name2i)+', ['+str(r)+', '+str(g)+', '+str(b)+']\n')
            results.write('color rgb'+str(res1i)+'-'+str(name1i)+'_'+str(res2i)+'-'+str(name2i)+', '+str(res1i)+'-'+str(name1i)+'_'+str(res2i)+'-'+str(name2i)+'\n')
    M=1
    for line in hmqclist[2:]:
        split=line.split()
        if split[4][0]!='N':
            split0=split[0].split('-')
            res1i=int(re.search('[A-Z]([0-9]+)[A-Z]', split0[0]).group(1))
            name1i=re.search('[0-9]+(\w+)', split0[0]).group(1)
            try:
              score=float(split[5])
              M4=float(M)
              M3=float(3*M/float(4))
              M2=float(2*M/float(4))
              M1=float(M/float(4))
              if (score>M4):r,g,b=0,0,1
              elif (score>M3):r,g,b=0,round(1-(score-M3)/(M4-M3), 2),1
              elif (score>M2):r,g,b=0,1,round((score-M2)/(M3-M2), 2)
              elif (score>M1):r,g,b=round(1-(score-M1)/(M2-M1), 2),1,0
              else:r,g,b=1,round(score/M1, 2),0
              if (split0[0][0]=='L' or split0[0][0]=='V'):star1='*'
              else:star1=''
              results.write('set_color rgb-'+str(res1i)+'-'+str(name1i)+', ['+str(r)+', '+str(g)+', '+str(b)+']\n')
              results.write('color rgb-'+str(res1i)+'-'+str(name1i)+', resi '+str(res1i)+' and name '+str(name1i)+star1+'\n')
            except ValueError:pass 
        else:pass
    results.write('hide labels\nset sphere_scale,0.5\nshow sphere, object Me')
    results.close() 
##########################################################################################
##########################################################################################
##########################################################################################
################################### MAIN STREAM ##########################################
#
Time_start=time.datetime.now()
try:
  if sys.argv[2]:directory=str(sys.argv[2])
except IndexError:directory=str(Time_start).split('.')[0]
directory = directory.replace(' ','_').replace(':','-')
os.makedirs('./'+directory)
os.makedirs('./'+directory+'/Input')
os.makedirs('./'+directory+'/Output')
try: 
  shutil.copy('./'+str(sys.argv[1]),'./'+directory+'/Input/'+str(sys.argv[1]))
except OSError:pass
try:
  shutil.copy('./'+str(sys.argv[0]),'./'+directory+'/Input/'+str(sys.argv[0]))
except OSError:pass

file_log=open('./'+directory+'/log', 'w')
file_log.write(directory+'\n')

file=open(str(sys.argv[1]),'r')
start_file=file.readlines()
file.close()

NOESY_list=start_file[4][:-1].split(';')
HMQC_name=start_file[2][:-1]
pdb_name=start_file[6][:-1]
Seq_name=start_file[8][:-1]
labeling=start_file[10]
flag_geminal=float(start_file[13][:-1])
cutoff_factor=float(start_file[15][:-1])
lowCut,max_distance=float(start_file[17].split()[0]),float(start_file[17].split()[1])
score_tol_end=str(start_file[19])

flag_L_nuclei,flag_V_nuclei='nd','nd'
for residue in labeling.split(';'):
  if (residue.split(',')[0]=='L' and len(residue.split(','))==2):flag_L_nuclei=residue.split(',')[1]
  elif (residue.split(',')[0]=='V' and len(residue.split(','))==2):flag_V_nuclei=residue.split(',')[1]

FILTER_scale=cutoff_factor

geminal_hold=0  #### On going test
try:
  shutil.copy('./'+str(HMQC_name),'./'+directory+'/Input/'+str(HMQC_name))
except OSError:pass
try:
  shutil.copy('./'+str(pdb_name),'./'+directory+'/Input/'+str(pdb_name))
except OSError:pass
try:
  shutil.copy('./'+str(Seq_name),'./'+directory+'/Input/'+str(Seq_name))
except OSError:pass

low_P,P_high=1,40
cpu=mp.cpu_count()
Tot_RAM=round(float(psutil.virtual_memory().total)/float(1024*1024*1024),1)
free_RAM_start=round(float(psutil.virtual_memory().available)/float(1024*1024*1024),1)
file_log.write('cpu: '+str(cpu)+' and total RAM: '+str(Tot_RAM)+' GB'+'\n')
file_size_max=int(1000000*(0.5*free_RAM_start-0.2*free_RAM_start)/float(2*cpu*0.06))
file_size=10000

HMQC_file=open('./'+str(HMQC_name),'r')
HMQC=HMQC_file.readlines()
HMQC_file.close()
pdb_file=open('./'+str(pdb_name),'r')
pdb=pdb_file.readlines()
pdb_file.close()
Seq_file=open('./'+str(Seq_name),'r')
SEQ=Seq_file.readlines()
Seq_file.close()

for NOESY_name in NOESY_list:
  NOESY_file=open('./'+str(NOESY_name),'r')
  NOESY=NOESY_file.readlines()
  NOESY_file.close()
  try:
    shutil.copy('./'+str(NOESY_name),'./'+directory+'/Input/'+str(NOESY_name))
  except OSError:pass
  nuclei=NOESY[0].split(';')
  tolerance=NOESY[1].split(';')
  if len(nuclei)==3:sym_coef=0      
  else:sym_coef=1
  DistNOE=Extract_3DPeaks(result,NOESY,HMQC)					
  factors=(1,0,0)
  peak_matrix_Scoring_i,peak_matrix_clustering_i,HMQC_peak_list,overlap_topo,biblio_crosspeaks,assignment_locked_peak,peak_geminal=noe2matrix(result,factors,0)
  C_matrix=open('./'+directory+'/C_matrix_'+str(NOESY_name),'w')
  np.savetxt(C_matrix,peak_matrix_clustering_i,fmt='%.3f')
  C_matrix.close()
  if NOESY_name==NOESY_list[0]:
    peak_matrix_Scoring,peak_matrix_clustering=peak_matrix_Scoring_i,peak_matrix_clustering_i
  else:
    peak_matrix_Scoring,peak_matrix_clustering=peak_matrix_Scoring+peak_matrix_Scoring_i,peak_matrix_clustering+peak_matrix_clustering_i
peak_matrix_clustering=(peak_matrix_clustering/float(len(NOESY_list)))+np.identity(len(HMQC_peak_list))

score_matrix=open('./'+directory+'/score_matrix','w')
np.savetxt(score_matrix,peak_matrix_Scoring,fmt='%.3f')
score_matrix.close()
C_matrix=open('./'+directory+'/C_matrix_TOT','w')
np.savetxt(C_matrix,peak_matrix_clustering,fmt='%.3f')
C_matrix.close()      

if len(HMQC_peak_list)>100:scaler=len(HMQC_peak_list)/float(100)
else:scaler=1

FILTER_start=np.mean(peak_matrix_Scoring.sum(axis=0))#*len(HMQC_peak_list)

file_log.write('Mean maximal peak score and scaler: '+str(round(FILTER_start,3))+'; '+str(scaler)+'\n')
calcul=peak_matrix_clustering[np.where(peak_matrix_clustering>0)]
main_confident=np.mean(calcul)
file_log.write('Average peak connection confidence score: '+str(round(main_confident,3))+'\n')

metrics=distances(pdb, max_distance,lowCut,'run')

methyl_to_add=[]
for line in SEQ:
  resID=re.search('\w(\d+)', line).group(1)
  resTYPE=line[0]
  flag=0
  for methyl in metrics[2]:
    methylID=re.search('\w(\d+)\w', methyl).group(1)
    if (methylID==resID and resTYPE==methyl[0]):
      flag=1
      break
  if flag==0:
    if resTYPE[0]=='A':name='CB'
    elif resTYPE[0]=='I':name='CD1'
    elif resTYPE[0]=='L':name='CD'
    elif resTYPE[0]=='M':name='CE'
    elif resTYPE[0]=='T':name='CG2'
    elif resTYPE[0]=='V':name='CG'
    methyl_to_add.append(resTYPE[0]+resID+name)
for i in methyl_to_add:
  metrics[2].append(i)
  if ((i[0]=='L' or i[0]=='V') and flag_geminal==1):
    if peak_geminal.sum()==0:metrics[3].append([i,i,1.2])
    else:metrics[3].append([i,i,1.0])
empty=[]
new_metrics=matrix_it(metrics[3],metrics[2],'pdb',empty,empty,empty)

methyl_file=open('./'+directory+'/Methyl_connectivity','w')
for line in metrics[3]:
  if line[2]!=0:methyl_file.write(str(line)+'\n')
methyl_file.close()

file_log.write('###################################'+'\n')
A,I,L,M,T,V=0,0,0,0,0,0
for i in range(len(metrics[2])):
  if metrics[2][i][0]=='A':A+=1
  elif metrics[2][i][0]=='I':I+=1
  elif metrics[2][i][0]=='L':L+=1
  elif metrics[2][i][0]=='M':M+=1
  elif metrics[2][i][0]=='T':T+=1
  elif metrics[2][i][0]=='V':V+=1
if flag_geminal!=2:number_of_methyls_from_pdb=A+I+2*L+M+T+2*V
else:number_of_methyls_from_pdb=A+I+L+M+T+V
file_log.write('Number of methyls: '+str(number_of_methyls_from_pdb)+'\n')
file_log.write('A: '+str(A)+'\nI: '+str(I)+'\nL: '+str(2*L)+'\nM: '+str(M)+'\nT: '+str(T)+'\nV: '+str(2*V)+'\n')

file_log.write('###################################'+'\n')
file_log.write('Number of peaks: '+str(len(HMQC))+'\n')
file_log.write('If stated for unique methyl type only:'+'\n')
a,il,l,m,t,v=0,0,0,0,0,0
total=0
for i in range(len(HMQC)):
  sum=0
  for j in HMQC[i].split()[3]:sum+=1
  total+=sum
  if sum==1:
    if 'A' in HMQC[i].split()[3]:a+=1
    elif 'I' in HMQC[i].split()[3]:il+=1
    elif 'L' in HMQC[i].split()[3]:l+=1
    elif 'M' in HMQC[i].split()[3]:m+=1
    elif 'T' in HMQC[i].split()[3]:t+=1
    elif 'V' in HMQC[i].split()[3]:v+=1
mean_MetType=round(total/float(len(HMQC)),1)
file_log.write('A: '+str(a)+'\nI: '+str(il)+'\nL: '+str(l)+'\nM: '+str(m)+'\nT: '+str(t)+'\nV: '+str(v)+'\n')
file_log.write('Methyl types per peak:'+str(mean_MetType)+'\n')
file_log.write('overlap_topo:'+'\n')
for i in range(len(overlap_topo)):
  if len(overlap_topo[i])>0:
    for j in overlap_topo[i]:file_log.write(str(HMQC_peak_list[i])+' '+str(HMQC_peak_list[int(j)])+'\n')
file_log.write('assignment_locked_peak: '+str(len(assignment_locked_peak))+'\n')
for i in assignment_locked_peak:
  file_log.write(str(HMQC_peak_list[i])+'\n')
file_log.write('\n\nCalculation starts:\n\n')
  
if (a>A or il>I or l/2>L or v/2>V or t>T or m>M):
  file_log.write('More peaks than methyls, please review the hmqc list'+'\n')
  print 'More peaks than methyls, please review the hmqc list'
  sys.exit()
else:pass
file_log.close()

report=open('./'+directory+'/Peak_list', 'w')
for i in range(len(HMQC_peak_list)): report.write(str(i)+' '+str(HMQC_peak_list[i])+'\n')
report.close()

report=open('./'+directory+'/Methyl_list', 'w')
for i in range(len(metrics[2])): report.write(str(i)+' '+str(metrics[2][i])+'\n')
report.close()

assignment_locked_methyl=[]
for i in assignment_locked_peak:
  name=HMQC_peak_list[i]
  assignment_name=name[0]+re.search('\w(\d+)', name).group(1)
  flag_assigned=0
  for j in range(len(metrics[2])):
    if metrics[2][j][0]+re.search('\w(\d+)', metrics[2][j]).group(1)==assignment_name:
      assignment_locked_methyl.append(j)
      flag_assigned=1
  if flag_assigned==0:
    file_log.write('issue with locked assignment'+'\n')
    file_log.write(str(assignment_name)+'\n')
    sys.exit()

##########################################################################################
################# Build up small redundant clusters ######################################
#
if not os.path.exists('./'+directory+'/run/Local'):os.makedirs('./'+directory+'/run/Local')
if not os.path.exists('./'+directory+'/run/Local/temp'):os.makedirs('./'+directory+'/run/Local/temp')
if not os.path.exists('./'+directory+'/run/Local/archive'):os.makedirs('./'+directory+'/run/Local/archive')

selected_peaks_clusters=obj()
N=peak_matrix_clustering.shape[0]
matrix2=np.dot(peak_matrix_clustering,peak_matrix_clustering)

for i in range(N):
  for j in range(N):
    if peak_matrix_clustering[i,j]==0:matrix2[i,j]/=2
  matrix2[i,i]=0

P_high=round(np.amax(matrix2),0)
P=float(P_high)

Density_matrix=open('./'+directory+'/Density_matrix','w')
np.savetxt(Density_matrix,matrix2,fmt='%.3f')
Density_matrix.close()

for i in range(N):
  neighbors_with_sharing=np.argwhere(matrix2[i,:]>P_high) # 2 => 1 shared neighbors, 3=>2, etc.
  neighbors_with_sharing=[neighbors_with_sharing[j,0] for j in range(neighbors_with_sharing.size)] 
  kept_neighbors_with_sharing=[]
  for j in neighbors_with_sharing:
    if peak_matrix_clustering[i,j]>0:kept_neighbors_with_sharing.append(j) 
  if not i in kept_neighbors_with_sharing:kept_neighbors_with_sharing.append(i) 
  kept_neighbors_with_sharing=np.array(kept_neighbors_with_sharing).reshape((len(kept_neighbors_with_sharing),1))
  setattr(selected_peaks_clusters, str(i), kept_neighbors_with_sharing)

list_peak=dir(selected_peaks_clusters)
list_peak.remove('__doc__')
list_peak.remove('__module__')
list_peak.sort(key=float)

for name_peak_index in range(len(list_peak)):
  peak_index=int(list_peak[name_peak_index])        
  counter=np.argwhere(matrix2[peak_index,:]>=5)
  cluster_size=int(str(len(counter))[0])
  if cluster_size<3:cluster_size=int(3)
  if P<5:cluster_size=int(3)

set_of_methyls=metrics[2]      
Ile_methyls=[]
Met_methyls=[]
Ala_methyls=[]
Thr_methyls=[]
VL_methyls=[]
Val_methyls=[]
Leu_methyls=[]
for i in range(len(set_of_methyls)):
        assignment=metrics[2][i]
        if (assignment[0]=='I' and not i in assignment_locked_methyl):Ile_methyls.append(i)
        if (assignment[0]=='M' and not i in assignment_locked_methyl):Met_methyls.append(i)
        if (assignment[0]=='A' and not i in assignment_locked_methyl):Ala_methyls.append(i)
        if (assignment[0]=='T' and not i in assignment_locked_methyl):Thr_methyls.append(i)
        if (assignment[0]=='V' and not i in assignment_locked_methyl):Val_methyls.append(i)
        elif (assignment[0]=='V' and assignment_locked_methyl.count(i)==1 and flag_geminal!=2):Val_methyls.append(i)
        if (assignment[0]=='L' and not i in assignment_locked_methyl):Leu_methyls.append(i)    
        elif (assignment[0]=='L' and assignment_locked_methyl.count(i)==1 and flag_geminal!=2):Leu_methyls.append(i)     
for name_peak in list_peak:   
      set_of_peaks=getattr(selected_peaks_clusters, name_peak)   
      previous_assignments=[[],[]]
      n=0
      possible_assignment=[]
      All_possible_assignment=[]
      for i in set_of_peaks:
        possible_assignment=[]
        code_methyl_type=HMQC[i[0]].split()[3]
        if 'A' in code_methyl_type:possible_assignment+=Ala_methyls
        if 'I' in code_methyl_type:possible_assignment+=Ile_methyls
        if 'L' in code_methyl_type:possible_assignment+=Leu_methyls
        if 'V' in code_methyl_type:possible_assignment+=Val_methyls      
        if 'M' in code_methyl_type:possible_assignment+=Met_methyls
        if 'T' in code_methyl_type:possible_assignment+=Thr_methyls
        All_possible_assignment.append(possible_assignment)
      previous_assignments=[[],[]]
      while n<len(set_of_peaks):
        total_Assignments=[]
        if set_of_peaks[n][0] in assignment_locked_peak:
          n+=1
          continue
        total_Assignments.append(previous_assignments[0]+[set_of_peaks[n][0]])
        for assignment in previous_assignments[1:]:
          for i in All_possible_assignment[n]:
            if assignment.count(i)==0:total_Assignments.append(assignment+[i])
            elif (flag_geminal!=2 and 
                  assignment.count(i)==1 and 
                  (metrics[2][i][0]=='L' or metrics[2][i][0]=='V')):
              total_Assignments.append(assignment+[i])
            else:continue
        previous_assignments=total_Assignments
        n+=1      
      assignment_archive=[]
      count=0
      locked_noe=[]
      locked_methyl=[]
      for i in range(len(assignment_locked_peak)):
        if assignment_locked_peak[i] in set_of_peaks:
          locked_noe=locked_noe+[assignment_locked_peak[i]]
          locked_methyl=locked_methyl+[assignment_locked_methyl[i]]
      if total_Assignments==[]:total_Assignments=[[],[]]
      index_matrix_noe=locked_noe+total_Assignments[0]
      matrix_noe=peak_matrix_Scoring[:,index_matrix_noe][index_matrix_noe,:]
      
      for assignment in total_Assignments[1:]:
        index_matrix_methyls=locked_methyl+assignment
        matrix_assignment=new_metrics[0][:,index_matrix_methyls][index_matrix_methyls,:]
        testing=matrix_assignment*matrix_noe
        tot=testing.sum()
        #if set_of_peaks.size==1:continue
        archive=[]
        archive.append(index_matrix_noe)
        archive.append(index_matrix_methyls)
        archive.append(tot)
        assignment_archive.append(archive)
        count=count+1             

      report=open('./'+directory+'/run/Local/c#'+str(name_peak)+'_P='+str(P), 'w')
           
      archive_assignment_cluster={}
      for i in range(len(assignment_archive)):
          for j in range(len(assignment_archive[i][0])):
            peak=HMQC_peak_list[assignment_archive[i][0][j]]
            methyl=metrics[2][assignment_archive[i][1][j]]
            if not peak in archive_assignment_cluster.keys():
              archive_assignment_cluster[peak]={}
              archive_assignment_cluster[peak][methyl]=assignment_archive[i][2]
            else:
              if not methyl in archive_assignment_cluster[peak].keys():
                archive_assignment_cluster[peak][methyl]=round(assignment_archive[i][2],3)
              elif (methyl in archive_assignment_cluster[peak].keys() and
                    archive_assignment_cluster[peak][methyl]<assignment_archive[i][2]):
                archive_assignment_cluster[peak][methyl]=round(assignment_archive[i][2],3)
              else:pass   
      setattr(result, str(name_peak)+'_possible_peak_assignments',archive_assignment_cluster)             

      score_index=np.zeros((2,len(assignment_archive)))      
      for i in range(len(assignment_archive)):
        assignment=assignment_archive[i]
        score_index[0,i]=assignment[2]
        score_index[1,i]=str(i)      
      score_sorting=score_index[0,:]
      sort=np.argsort(score_sorting, axis=None)
      sort=sort[::-1] #mirror flip on both dimensions
      for i in sort[:]:
        assignment=assignment_archive[i]
        report.write('#'+str(int(score_index[1,i]))+': '+str(assignment[2])+' ### ')
        for j in range(len(assignment[0])):
          peak=HMQC_peak_list[assignment[0][j]]
          methyl=metrics[2][assignment[1][j]]
          report.write(str(peak)+'=='+str(methyl)+' # ')
        report.write('\n')
      report.close()
      file=open('./'+directory+'/run/Local/archive/archive_'+str(name_peak), 'w')
      cPickle.dump(assignment_archive,file,-1)
      file.close()
      assignment_archive=[]      
######################## iterative RUNs ####################
#

P10=round(P_high*0.1,0)
P5=round(P_high*0.2,0)
P3=round(P_high*0.33,0)
P2=round(P_high*0.5,0)
filling_id=np.argwhere(matrix2==0)
matrix2_edited=np.array(matrix2)
matrix2_edited[filling_id[0]]=2
matrix2_edited[filling_id[1]]=3
matrix2_edited[filling_id[2]]=5
matrix2_edited[filling_id[3]]=P10
matrix2_edited[filling_id[4]]=P5
matrix2_edited[filling_id[5]]=P3
P_list=np.unique(matrix2_edited)
P_list=P_list[::-1]

selected_peaks_clusters=obj()
peak_clusters=obj()
N=peak_matrix_clustering.shape[0]
List_of_assigned_peaks=[]
already_called_cluster=[]
postpone_peaks=[]
delayed_clusters=[]
flag_score_equation=0
archive_assignment_result={}
completeness=np.zeros((len(HMQC_peak_list),2))

for i in range(len(HMQC_peak_list)):completeness[i,1]=len(biblio_crosspeaks[HMQC_peak_list[i]].keys())


for P in P_list:
    if P>2:
      flag_new_peaks=0
      for i in range(N):
        neighbors_with_sharing=np.argwhere(matrix2[i,:]>P)
        neighbors_with_sharing=[neighbors_with_sharing[j,0] for j in range(neighbors_with_sharing.size)] 
        kept_neighbors_with_sharing=[]
        for j in neighbors_with_sharing:
          if peak_matrix_clustering[i,j]>0:kept_neighbors_with_sharing.append(j) 
        if not i in kept_neighbors_with_sharing:kept_neighbors_with_sharing.append(i) 
        kept_neighbors_with_sharing=np.array(kept_neighbors_with_sharing).reshape((len(kept_neighbors_with_sharing),1))
        setattr(selected_peaks_clusters, str(i), kept_neighbors_with_sharing)

        
      list_peak=dir(selected_peaks_clusters)
      list_peak.remove('__doc__')
      list_peak.remove('__module__')
      list_peak.sort(key=float)
      for name_peak_index in range(len(list_peak)):
        highest_score_small=0
        peak_index=int(list_peak[name_peak_index])        
        
        ####Define cluster size       
        counter=np.argwhere(matrix2[peak_index,:]>=3)
        cluster_size=int(str(len(counter))[0])
        #if cluster_size>5:cluster_size=int(5)
        if P>P10 and cluster_size<int(5):cluster_size=int(5)
        if P<P10 and cluster_size>int(5):cluster_size=int(5)
        if (cluster_size<4 or P<3):cluster_size=int(4)
        if P==1:cluster_size=int(3)
         
        if peak_index in already_called_cluster:continue
        elif peak_index in delayed_clusters:continue
        #elif len(overlap_topo[peak_index])>0:continue
        #######IN: Skip peak assignment if overlapping already assigned peak #########
        flag_continue=0
        for overpeak in overlap_topo[peak_index]:
              if (overpeak in already_called_cluster or overpeak in assignment_locked_peak):flag_continue=1
        if flag_continue==1:continue
        #######OUT: Skip peak assignment if overlapping already assigned peak ########       
        set_of_peaks=getattr(selected_peaks_clusters, str(peak_index))
        file=open('./'+directory+'/run/Local/archive/archive_'+str(peak_index), 'r')
        assignment_archive_old=cPickle.load(file)
        file.close()
        
        if len(assignment_archive_old)==0:continue
                         
        try:
          if (len(assignment_archive_old[0][0])>=cluster_size and P<P3):continue
          if (len(assignment_archive_old[0][0])==1 and P<3):continue						####EDIT 12/11/2019
          new_peaks=[]
          for peak in set_of_peaks:
            if not peak in assignment_archive_old[0][0]:
            #if ((not peak in assignment_archive_old[0][0]) and 
            #    (len(overlap_topo[peak[0]])==0)
            #   ):new_peaks.append(peak[0])
              flag_continue=0
              for overpeak in overlap_topo[peak[0]]:
                if (overpeak in assignment_archive_old[0][0] or overpeak in new_peaks):flag_continue=1
              if flag_continue==0:new_peaks.append(peak[0])          
          setattr(result, str(peak_index)+'_new_peaks',new_peaks)
        except IndexError:continue
        
        if (len(new_peaks)==0):continue
        file_log=open('./'+directory+'/log', 'a')
        file_log.write('Local: Tc= '+str(round(P,2))+' ; Clustering around peak # '+str(peak_index)+'\n')
        file_log.close()
        
        if flag_new_peaks==1:pass
        elif set_of_peaks.size>=cluster_size:flag_new_peaks=1
        else:flag_new_peaks=0

        sorting=[(peak_matrix_clustering[peak,peak_index],peak) for peak in new_peaks]
        sorting.sort()        
        sorting=sorting[::-1]
               
        for line in sorting:
          if line!=sorting[-1]:flag_fin=0			
          else:flag_fin=1
          peak=line[1]
          file=open('./'+directory+'/run/Local/archive/archive_'+str(peak_index), 'r')
          assignment_archive_old=cPickle.load(file)
          file.close()
          All_possible_assignment=[]  
          possible_assignment=[]
          if peak in assignment_locked_peak:
            for i in range(len(assignment_locked_peak)):
              if int(assignment_locked_peak[i])==int(peak):All_possible_assignment.append([assignment_locked_methyl[i]])
          else:
            code_methyl_type=HMQC[peak].split()[3]
            if 'A' in code_methyl_type:possible_assignment+=Ala_methyls
            if 'I' in code_methyl_type:possible_assignment+=Ile_methyls
            if 'L' in code_methyl_type:possible_assignment+=Leu_methyls
            if 'V' in code_methyl_type:possible_assignment+=Val_methyls      
            if 'M' in code_methyl_type:possible_assignment+=Met_methyls
            if 'T' in code_methyl_type:possible_assignment+=Thr_methyls
            All_possible_assignment.append(possible_assignment)
          total_Assignments=[[peak]]
          for i in All_possible_assignment[0]:total_Assignments.append([i])
          #################### Multi-processed assignment building #####################
          total=len(assignment_archive_old)
          if total==0:continue
          m=(int(total*len(total_Assignments)))/50000
          m_reste=(int(total*len(total_Assignments)))%50000       
          if m<cpu:m=cpu          
          if total<m:
            size=1
            m=total
          else:size=total/(m)
          reste=total%(m)
          list_of_assignment_index=[]
          if size!=0:
            for i in range(m):
              list_of_assignment_index.append([j for j in range(size*i,(i+1)*size)])
          if reste!=0:list_of_assignment_index.append([j for j in range(m*size,m*size+reste)])   
          number_of_pool=len(list_of_assignment_index)/cpu
          number_of_pool_rest=len(list_of_assignment_index)%cpu
          number_of_files=0
          for i in range(number_of_pool):
              pool=mp.Pool(processes=cpu) 
              pool_result=pool.map(build_assignment_peak, list_of_assignment_index[cpu*i:cpu*(i+1)])
              pool.close()

              for j in range(len(pool_result)):
                if  pool_result[j][1]>highest_score_small:highest_score_small=pool_result[j][1]
              assignments_table=[]
              for j in range(len(pool_result)):
                if len(pool_result[j][0])!=0:assignments_table.append(pool_result[j][0])
              pool_result=[]
              if len(assignments_table)!=0:
                file=open('./'+directory+'/run/Local/temp/temp_'+str(number_of_files), 'w')
                cPickle.dump(assignments_table,file,-1)
                file.close()
                assignments_table=[]
                number_of_files+=1 
          if number_of_pool_rest!=0:
              pool = mp.Pool(processes=cpu) 
              pool_result=pool.map(build_assignment_peak, list_of_assignment_index[number_of_pool*cpu:]) 
              pool.close()
              
              for j in range(len(pool_result)):
                if  pool_result[j][1]>highest_score_small:highest_score_small=pool_result[j][1]
              assignments_table=[]
              for j in range(len(pool_result)):
                if len(pool_result[j][0])!=0:assignments_table.append(pool_result[j][0])
              pool_result=[]

              if len(assignments_table)!=0:
                file=open('./'+directory+'/run/Local/temp/temp_'+str(number_of_files), 'w')
                cPickle.dump(assignments_table,file,-1)
                file.close()
                assignments_table=[]
                number_of_files+=1
          ##################### Multi-processed assignment building ####################  
          number_peaks=len(assignment_archive_old[0][0])+1
          if number_peaks==2:q=0          
          elif number_peaks==3:q=0.01
          elif number_peaks==4:q=0.7
          elif number_peaks==5:q=0.8
          elif number_peaks>=6:q=0.9
          if line!=sorting[-1]:q*=0.9        
          if mean_MetType>1:q*=0.9
          
          if FILTER_scale>1:scale_cluster=FILTER_scale
          else:scale_cluster=1
          cutoff=highest_score_small-highest_score_small*(1-q)*scale_cluster
          tot=0
          deleted=0
          
          if highest_score_small==0.0:cutoff=1.0

          assignment_archive=[]
          list_of_files=os.listdir('./'+directory+'/run/Local/temp/')
          for i in range(len(list_of_files)):
              file=open('./'+directory+'/run/Local/temp/'+str(list_of_files[i]), 'r')
              element=cPickle.load(file)
              file.close()
              for j in range(len(element)):
                for k in range(len(element[j])):
                  tot+=1
                  if element[j][k][2]>=cutoff:assignment_archive.append(element[j][k])
              os.remove('./'+directory+'/run/Local/temp/'+str(list_of_files[i]))  
          deleted=tot-len(assignment_archive)+1
        
          archive_assignment_cluster={}
          for i in range(len(assignment_archive)):
            for j in range(len(assignment_archive[i][0])):
              peak=HMQC_peak_list[assignment_archive[i][0][j]]
              methyl=metrics[2][assignment_archive[i][1][j]]
              if not peak in archive_assignment_cluster.keys():
                archive_assignment_cluster[peak]={}
                archive_assignment_cluster[peak][methyl]=assignment_archive[i][2]
              else:
                if not methyl in archive_assignment_cluster[peak].keys():
                  archive_assignment_cluster[peak][methyl]=round(assignment_archive[i][2],3)
                elif (methyl in archive_assignment_cluster[peak].keys() and
                      archive_assignment_cluster[peak][methyl]<assignment_archive[i][2]):
                  archive_assignment_cluster[peak][methyl]=round(assignment_archive[i][2],3)
                else:pass
          setattr(result, str(list_peak[name_peak_index])+'_possible_peak_assignments',archive_assignment_cluster)
          
          report=open('./'+directory+'/run/Local/c#'+str(peak_index)+'_P='+str(round(P,2)), 'w')        
          score_index=np.zeros((2,len(assignment_archive)))   
          for i in range(len(assignment_archive)):
            assignment=assignment_archive[i]
            score_index[0,i]=assignment[2]
            score_index[1,i]=str(i)      
          score_sorting=score_index[0,:]
          sort=np.argsort(score_sorting, axis=None)
          sort=sort[::-1] #mirror flip on both dimensions
          report.write('#cutoff= '+str(cutoff)+'\n')
          for i in sort[:]:
            assignment=assignment_archive[int(score_index[1,i])]
            report.write('#'+str(int(score_index[1,i]))+': '+str(assignment[2])+' ### ')
            for j in range(len(assignment[0])):
              peak=HMQC_peak_list[assignment[0][j]]
              methyl=metrics[2][assignment[1][j]]
              report.write(str(peak)+'=='+str(methyl)+' # ')
            report.write('\n')
          report.close()
          file=open('./'+directory+'/run/Local/archive/archive_'+str(peak_index), 'w')
          cPickle.dump(assignment_archive,file,-1)
          file.close()
          assignment_archive=[] 
#####################################################################################################      
###################### Mix all clusters together #################################
    else:pass
    #if P>=2:
    if P>=1:
      if (flag_new_peaks==1 or P==3 or P==5 or P==P5 or P==P10 or P==P3 or P==2 or P==1):pass
      else:continue

      logFile_name='logFile_'+str(round(P,2))
      path_name='/run/Global'
      if not os.path.exists('./'+directory+'/run/Global'):os.makedirs('./'+directory+'/run/Global')
      if not os.path.exists('./'+directory+'/run/temp'):
        
        possible_peak_assignments_total={}######ATTENTION  
        
        os.makedirs('./'+directory+'/run/temp')
        assignment_archive_FINAL=[assignment_locked_peak,[assignment_locked_methyl,]]
        file=open('./'+directory+'/run/temp/100000000#', 'w')
        cPickle.dump(assignment_archive_FINAL,file,-1)
        file.close()
        highest_score=0
        flag_alone=0
      
      list_peak=dir(selected_peaks_clusters)
      list_peak.remove('__doc__')
      list_peak.remove('__module__')
      list_peak.sort(key=float)
      sort_index=np.zeros((2,len(list_peak)))          
      for i in range(len(list_peak)):
            file=open('./'+directory+'/run/Local/archive/archive_'+str(i), 'r')
            assignment_archive_clusterX=cPickle.load(file)
            file.close()
            sort_index[0,i]=list_peak[i]
            sort_index[1,i]=len(assignment_archive_clusterX)
      score_sorting=sort_index[1,:]
      sort=np.argsort(score_sorting, axis=None)
      flag_print_result=0
      
      for cluster_index in sort:
            cluster_name=int(sort_index[0,cluster_index])
            set_of_peaks=getattr(selected_peaks_clusters, str(cluster_name))
            possible_peak_assignments=getattr(result, str(cluster_name)+'_possible_peak_assignments')
            number_of_peaks=len(possible_peak_assignments.keys())
 
            ####Define cluster size       
            counter=np.argwhere(matrix2[cluster_name,:]>=3)
            cluster_size=int(str(len(counter))[0])
            #if cluster_size>5:cluster_size=int(5)
            if P>P10 and cluster_size<int(5):cluster_size=int(5)
            if P<P10 and cluster_size>int(5):cluster_size=int(5)
            if (cluster_size<4 or P<3):cluster_size=int(4)
            if P==2:cluster_size=int(3)
                       
            if (number_of_peaks<cluster_size 
                or cluster_name in already_called_cluster
                or (cluster_name in delayed_clusters and P!=P5 and P!=P10 and P!=5 and P!=3 and P!=2 and P!=1)
                ):continue                     
            else:pass

            #######IN: Skip peak assignment if overlapping already assigned peak #########
            flag_continue=0
            for overpeak in overlap_topo[cluster_name]:
              if (overpeak in already_called_cluster or overpeak in assignment_locked_peak):flag_continue=1
            if flag_continue==1:continue
            #######OUT: Skip peak assignment if overlapping already assigned peak ########
            
            file=open('./'+directory+'/run/Local/archive/archive_'+str(cluster_name), 'r')
            assignment_archive_clusterX=cPickle.load(file)
            file.close()

            list_of_files=os.listdir('./'+directory+'/run/temp/')
            file=open('./'+directory+'/run/temp/'+str(list_of_files[0]), 'r')
            assignment_previous=cPickle.load(file)
            file.close()
            
            try:involved_peak=assignment_archive_clusterX[0][0]
            except IndexError:continue
            index_to_keep=[]
            for i in range(len(involved_peak)):
              if not involved_peak[i] in assignment_previous[0]:index_to_keep.append(i)          
            
            if (len(index_to_keep)==len(involved_peak) and 
                len(involved_peak)>=cluster_size):flag_newArea=1
            else:flag_newArea=0
            
            list2stock=[]															#EDIT
            involved_peak=[involved_peak[i] for i in index_to_keep]
            possible_peak_assignments=[]
            if len(involved_peak)==0:
              already_called_cluster.append(cluster_name)
              if cluster_name in delayed_clusters:delayed_clusters.remove(cluster_name)
              continue            
            for line in assignment_archive_clusterX:
              #####IN: Test allowed assignment for all peaks #######
              #if P<=P10:
              if P<=P_high:
                try:
                  flag_stop=0
                  for i in range(len(line[1])):
                    peak=HMQC_peak_list[int(assignment_archive_clusterX[0][0][i])]
                    methyl=metrics[2][int(line[1][i])]
                    if (peak in archive_assignment_result.keys() and
                       len(overlap_topo[int(assignment_archive_clusterX[0][0][i])])==0):
                      if not methyl in archive_assignment_result[peak].keys():
                        flag_stop=1
                        break
                      else:pass
                    else:pass
                  if flag_stop==1:continue
                except NameError:pass
              else:pass
              #####OUT: Test allowed assignment for all peaks #######               
              line_to_add=[line[1][i] for i in index_to_keep]
              #####IN: Test whether peak assignment is not already used #######
              flag_skip_line=0
              for i in line_to_add:
                methyl=metrics[2][int(i)]
                for peak in archive_assignment_result.keys():
                  if methyl in archive_assignment_result[peak].keys():
                    if (len(archive_assignment_result[peak].keys())==1 and flag_geminal==2):
                      flag_skip_line=1
                      break
                    if (len(archive_assignment_result[peak].keys())==1 and methyl[0]!='L' and methyl[0]!='V'):
                      flag_skip_line=1
                      break
                    else:
                      excluding_assignment=archive_assignment_result[peak].keys()
                      count=0
                      for methyl_ex in excluding_assignment:
                        if (methyl_ex[0]=='L' or methyl_ex[0]=='V'):count+=-1     
                      for peak2 in archive_assignment_result.keys():
                        if excluding_assignment==archive_assignment_result[peak2].keys():count+=1
                        else:pass
                        if count==len(excluding_assignment):
                          flag_skip_line=1
                          break
                    if flag_skip_line==1:break  
                if flag_skip_line==1:break
              if flag_skip_line==1:continue    
              #####OUT: Test whether peak assignment is not already used #######
              if not line_to_add in possible_peak_assignments:							#EDIT
                possible_peak_assignments.append(line_to_add)							#EDIT
                list2stock.append([involved_peak,line_to_add,line[2]])					#EDIT

            #EDIT
            file=open('./'+directory+'/run/Local/archive/archive_'+str(cluster_name), 'w')
            cPickle.dump(list2stock,file,-1)
            file.close()
            #EDIT          
                          
            if (len(possible_peak_assignments)==0):
                if not cluster_name in delayed_clusters:delayed_clusters.append(cluster_name)
                continue            
            else:pass

            sum_local=peak_matrix_Scoring[involved_peak,:][:,List_of_assigned_peaks].sum()
            sum_global=peak_matrix_Scoring[List_of_assigned_peaks,:][:,involved_peak].sum()
            ratio=round(float((len(List_of_assigned_peaks)+len(involved_peak)+len(assignment_locked_peak))/float(len(list_peak))),2)
            
            peak2use=List_of_assigned_peaks+involved_peak
            TOTscore_mean=np.mean(peak_matrix_Scoring[peak2use,:][:,peak2use].sum(axis=0))

            FILTER=round(scaler*FILTER_scale*TOTscore_mean*(1-ratio*0.9),3) #(1-ratio*0.85)
            #FILTER_gap=round(scaler*TOTscore_mean*(1-ratio*0.9),3) #*(1-ratio*0.85)   
            keep2test=[] 
            for peak in list_peak:
              if not int(peak) in involved_peak:keep2test.append(int(peak))
              else:pass
            sum_out=peak_matrix_Scoring[involved_peak,:][:,keep2test].sum()
            sum_tot=sum_out+sum_local+sum_global
            
            if len(involved_peak)==1:max_len=10		#5
            elif len(involved_peak)==2:max_len=20	#5x4
            elif len(involved_peak)==3:max_len=60	#5x4x3
            elif len(involved_peak)==4:max_len=120	#5x4x3x2
            else:max_len=120*(1.5**(len(involved_peak)-4))        
            
            #if P<=3:max_len*=2										#edit du 10/01/2020
            if flag_newArea==1 and P>5:max_len*=3							
            max_len*=mean_MetType**(len(involved_peak))		        #edit du 10/01/2020
            if (sum_local+sum_global)>FILTER:max_len*=((sum_local+sum_global)/FILTER)                    
            #if sum_tot>FILTER:max_len*=(sum_tot/FILTER)        
            max_len=int(round(max_len,0)) 
                        
            file_test=open('./'+directory+'/test', 'a')
            file_test.write('P: '+str(round(P,2))+(7-len(str(round(P,2))))*' '+
                            '\tclus: '+str(cluster_name)+(5-len(str(cluster_name)))*' '+
                            '\t#peak: '+str(len(involved_peak))+(3-len(str(len(involved_peak))))*' '+
                            '\tsum_out: '+str(round(sum_out,3))+(7-len(str(round(sum_out,3))))*' '+
                            '\tsum_g+l= '+str(round(sum_global+sum_local,3))+(7-len(str(round(sum_global+sum_local,3))))*' '+
                            '\ttol= '+str(round(FILTER,3))+(7-len(str(round(FILTER,3))))*' '+
                            '\tass= '+str(len(possible_peak_assignments))+
                            '\t/ '+str(max_len)+(7-len(str(len(possible_peak_assignments))+'\t/'+str(max_len)))*'  '+
                            '\tflag= '+str(flag_newArea))
            
            if ((P>P3 and len(possible_peak_assignments)<=max_len)
                or len(List_of_assigned_peaks)==0
                or len(possible_peak_assignments)==1
                or (sum_tot>=FILTER and len(possible_peak_assignments)<=max_len)
                or (P<=3 and len(possible_peak_assignments)<=max_len)
                ):
               file_test.write(' *\n')
               file_test.close()
               pass
            else:  
               if not cluster_name in delayed_clusters:delayed_clusters.append(cluster_name)
               file_test.write('\n')
               file_test.close()
               continue                   

            file_log=open('./'+directory+'/log', 'a')
            file_log.write('Global: Tc= '+str(round(P,1))+' ; Cluster merging # '+str(cluster_name)+'\n')
            file_log.close()      
            report_follow=open('./'+directory+'/run/Global/logFile_'+str(round(P,2)), 'a')
            report_follow.write('##########  Cluster : '+str(cluster_name)+'\n')
            report_follow.write('involved_peak: '+str(involved_peak)+'\n'+'number of possible assignments: '+
                             str(len(possible_peak_assignments))+'\n'+str(possible_peak_assignments)+'\n')
            report_follow.close()

            #####################IN:Multi-processed assignment building###################
            list_of_files=os.listdir('./'+directory+'/run/temp/')
            archive_assignment_result_previous=archive_assignment_result
            archive_assignment_result={}
            shuffle(list_of_files)
            size=len(list_of_files)/cpu
            size_reste=len(list_of_files)%cpu 
            list_of_file_index=[]
            if size!=0:
              for i in range(cpu):
                if i<size_reste:list_of_file_index.append([j for j in range(size*i,(i+1)*size)]+[int(len(list_of_files)-i-1)])
                else:list_of_file_index.append([j for j in range(size*i,(i+1)*size)])
              pool=mp.Pool(processes=cpu) 
              pool_result=pool.map(build_assignment, list_of_file_index)
              pool.close() 
            elif size==0:
              for i in range(len(list_of_files)):list_of_file_index.append([i])
              pool=mp.Pool(processes=len(list_of_file_index)) 
              pool_result=pool.map(build_assignment, list_of_file_index)
              pool.close()       
            for j in range(len(pool_result)):
              if  pool_result[j][1]>highest_score:highest_score=pool_result[j][1]                 
            for element in pool_result:
              archive=element[0]
              for peak in archive.keys():
                if not peak in archive_assignment_result.keys():archive_assignment_result[peak]={}
                for methyl in archive[peak].keys():
                  if not methyl in archive_assignment_result[peak].keys():
                    archive_assignment_result[peak][methyl]=archive[peak][methyl]
                  if (methyl in archive_assignment_result[peak].keys() and
                      archive[peak][methyl]>archive_assignment_result[peak][methyl]):
                    archive_assignment_result[peak][methyl]=archive[peak][methyl]                              
            #####################OUT:Multi-processed assignment building##################       
            list_of_files=os.listdir('./'+directory+'/run/temp/')
            flag_not_null=0
            for file_name in list_of_files:
              peak_id=re.search('^([0-9]+)#', file_name).group(1)
              if int(peak_id)==cluster_name:
                flag_not_null=1
                flag_print_result=1
                break
            if flag_not_null==1:
              for file_name in list_of_files:
                peak_id=re.search('^([0-9]+)#', file_name).group(1)
                if int(peak_id)!=cluster_name:os.remove('./'+directory+'/run/temp/'+str(file_name))
              for i in involved_peak: 
                if not i in List_of_assigned_peaks:List_of_assigned_peaks.append(i)              
            else:archive_assignment_result=archive_assignment_result_previous
        
            already_called_cluster.append(cluster_name)
            if cluster_name in delayed_clusters:delayed_clusters.remove(cluster_name) 
            report_follow=open('./'+directory+'/run/Global/logFile_'+str(round(P,2)), 'a')
            for peak in archive_assignment_result.keys():
              for methyl in archive_assignment_result[peak]:
                try:
                  if re.search('\w(\d+)\w',str(methyl)).group(1)==re.search('(\d+)\w', peak).group(1):
                    if archive_assignment_result[peak][methyl]==highest_score:report_follow.write('R* ')
                    else:report_follow.write('R  ')
                    break
                except AttributeError:pass
              report_follow.write(str(peak)+'  '+str(archive_assignment_result[peak])+'\n')     
            report_follow.close()
               
            if flag_not_null==1:
              outputing(highest_score,peak_matrix_Scoring)  
              Time_end=time.datetime.now()
              list_of_peaks=archive_assignment_result.keys()    
              print 'Tc= '+str(round(P,2))+'; '+str(Time_end-Time_start).split('.')[0].split('.')[0],' => ',str(round(100*len(list_of_peaks)/float(len(HMQC_peak_list)),1)),'%'                 
            setattr(result,'archive_assignment_result',archive_assignment_result)

##########################################################################################
############ Assignment of isolated peak #################################################
P='iso'
if score_tol_end[1]=='f':pass
else:
  #### Define not yet assigned peaks ####
  isolated_peaks_list=[]
  for i in range(len(HMQC_peak_list)):
    if not i in List_of_assigned_peaks:isolated_peaks_list.append(i)
    else:pass

  path_name='/run/Global'
  logFile_name='logFile_iso'
  list_of_files=os.listdir('./'+directory+'/run/temp/')
  
  FILTER_iso=float(FILTER)
  
  #print highest_score
  #print 'before',peak_matrix_Scoring.sum()
  ####### Recalculate score to take in account low confident connection
  #file_log=open('./'+directory+'/log', 'w')
  #for NOESY_name in NOESY_list:
  #  NOESY_file=open('./'+str(NOESY_name),'r')
  #  NOESY=NOESY_file.readlines()
  #  NOESY_file.close()
  #  nuclei=NOESY[0].split(';')
  #  tolerance=NOESY[1].split(';') 
  #  if len(nuclei)==3:sym_coef=0      
  #  else:sym_coef=1
  #  DistNOE=Extract_3DPeaks(result,NOESY,HMQC)					
  #  factors=(1,0,0)
  #  peak_matrix_Scoring_i,peak_matrix_clustering_i,HMQC_peak_list,overlap_topo,biblio_crosspeaks,assignment_locked_peak,peak_geminal=noe2matrix(result,factors,1)
  #  if NOESY_name==NOESY_list[0]:
  #    peak_matrix_Scoring,peak_matrix_clustering=peak_matrix_Scoring_i,peak_matrix_clustering_i
  #  else:
  #    peak_matrix_Scoring,peak_matrix_clustering=peak_matrix_Scoring+peak_matrix_Scoring_i,peak_matrix_clustering+peak_matrix_clustering_i
  #peak_matrix_clustering=(peak_matrix_clustering/float(len(NOESY_list)))+np.identity(len(HMQC_peak_list))
  #file_log.close()		#moyenne pour clustering mais non pour scoring ???
  ######END
  #print 'after',peak_matrix_Scoring.sum()
         
  still_not_assigned=list(isolated_peaks_list)
  no_crosspeaks_peak=[]
  cycler=0
  turn=0
  a=0
  list_peak=dir(selected_peaks_clusters)
  list_peak.remove('__doc__')
  list_peak.remove('__module__')
  list_peak.sort(key=float)

  while cycler==0:
    a+=1
    turn+=1
    sort_index=np.zeros((2,len(still_not_assigned)))
    for i in range(len(still_not_assigned)):
        sort_index[0,i]=still_not_assigned[i]
        sort_index[1,i]=np.amax(peak_matrix_clustering[still_not_assigned[i],List_of_assigned_peaks])	
    score_sorting=sort_index[1,:]
    sort=np.argsort(score_sorting, axis=None)
    sort=sort[::-1]
    #### Accretion to final_assignment ####
    for index in sort:
      cluster_name=int(sort_index[0,index])
      file_log=open('./'+directory+'/log', 'a')
      file_log.write('Isolated peak # '+str(cluster_name)+'\n')
      
      previous_assignments=[[],[]]
      n=0
      possible_assignment=[]
      All_possible_assignment=[]
      possible_assignment=[]
      code_methyl_type=HMQC[cluster_name].split()[3]
      if 'A' in code_methyl_type:possible_assignment+=Ala_methyls
      if 'I' in code_methyl_type:possible_assignment+=Ile_methyls
      if 'L' in code_methyl_type:possible_assignment+=Leu_methyls
      if 'V' in code_methyl_type:possible_assignment+=Val_methyls      
      if 'M' in code_methyl_type:possible_assignment+=Met_methyls
      if 'T' in code_methyl_type:possible_assignment+=Thr_methyls
      possible_peak_assignments=[[i] for i in possible_assignment]
      involved_peak=[cluster_name]     
      possible_peak_assignments2keep=[]
      for line in possible_peak_assignments:
          #####IN: Test whether peak assignment is not already used #######
          flag_skip_line=0
          methyl=metrics[2][int(line[0])]
          for peak in archive_assignment_result.keys():
            if methyl in archive_assignment_result[peak].keys():
              if (len(archive_assignment_result[peak].keys())==1 and flag_geminal==2):
                flag_skip_line=1
                break
              if (len(archive_assignment_result[peak].keys())==1 and methyl[0]!='L' and methyl[0]!='V'):
                flag_skip_line=1
                break
              else:
                excluding_assignment=archive_assignment_result[peak].keys()
                count=0
                for methyl_ex in excluding_assignment:
                  if (methyl_ex[0]=='L' or methyl_ex[0]=='V'):count+=-1 
                for peak2 in archive_assignment_result.keys():
                  if excluding_assignment==archive_assignment_result[peak2].keys():count+=1
                  else:pass
                  if count==len(excluding_assignment):
                    flag_skip_line=1
                    break
            if flag_skip_line==1:break  
          if flag_skip_line==1:continue
          #####OUT: Test whether peak assignment is not already used #######
          if not line in possible_peak_assignments2keep:possible_peak_assignments2keep.append(line)
      possible_peak_assignments=list(possible_peak_assignments2keep)
      report_follow=open('./'+directory+'/'+path_name+'/'+logFile_name, 'a')
      report_follow.write('##########  Peak alone : '+str(cluster_name)+'\n')
      report_follow.write('assignment table length: '+str(len(possible_peak_assignments))+'\n'+str(possible_peak_assignments)+'\n')
      report_follow.close() 

      if peak_matrix_Scoring[involved_peak[0],:].sum()==0:
          file_log.write('no crosspeaks'+'\n')
          no_crosspeaks_peak.append([cluster_name,possible_peak_assignments])
          continue         

      sum_local=peak_matrix_Scoring[involved_peak,:][:,List_of_assigned_peaks].sum()
      sum_global=peak_matrix_Scoring[List_of_assigned_peaks,:][:,involved_peak].sum()
      ratio=round(float((len(List_of_assigned_peaks)+1+len(assignment_locked_peak))/float(len(list_peak))),2)      
      peak2use=List_of_assigned_peaks+involved_peak
      TOTscore_mean=np.mean(peak_matrix_Scoring[peak2use,:][:,peak2use].sum(axis=0))
      FILTER=round(scaler*FILTER_scale*TOTscore_mean*(1-ratio*0.9),3) #(1-ratio*0.85)
      #FILTER_gap=round(scaler*TOTscore_mean*(1-ratio*0.9),3) #*(1-ratio*0.85)
      
      #To keep tolerance higher enough to accomodate less accurate connections
      if FILTER<round(TOTscore_mean/float(3),3) and FILTER_iso<round(TOTscore_mean/float(3),3):FILTER=FILTER_iso
      elif FILTER<round(TOTscore_mean/float(3),3) and FILTER_iso>round(TOTscore_mean/float(3),3):FILTER=round(TOTscore_mean/float(3),3)
      else:pass
          
      keep2test=[] 
      for peak in list_peak:
        if not int(peak) in involved_peak:keep2test.append(int(peak))
        else:pass
      sum_out=peak_matrix_Scoring[involved_peak,:][:,keep2test].sum()
      
      file_test=open('./'+directory+'/test', 'a')
      file_test.write('ISO   '+
                      '\tclus: '+str(cluster_name)+(5-len(str(cluster_name)))*' '+
                      '\t#peak: '+str(len(involved_peak))+(3-len(str(len(involved_peak))))*' '+
                      '\tsum_out: '+str(round(sum_out,3))+(7-len(str(round(sum_out,3))))*' '+
                      '\tsum_g+l= '+str(round(sum_global+sum_local,3))+(7-len(str(round(sum_global+sum_local,3))))*' '+
                      '\ttol= '+str(round(FILTER,3))+(7-len(str(round(FILTER,3))))*' '+
                      '\tass= '+str(len(possible_peak_assignments))) 
      if (((sum_global+sum_local)<FILTER and turn==1)
          or ((sum_global+sum_local)<FILTER and sum_out<FILTER and turn==2)
          or ((sum_global+sum_local+sum_out)<FILTER and turn==3)
          #((sum_global+sum_local+sum_out)<(FILTER/float(2)) and turn==4)  #peut etre annule si trop long
          ):
          file_log.write('too ambiguous\n')
          file_test.write('\n')
          continue
      file_test.write('*\n')
      file_test.close()
      ##################### Multi-processed assignment building ####################         
      archive_assignment_result_previous=archive_assignment_result
      archive_assignment_result={}
      list_of_files=os.listdir('./'+directory+'/run/temp/')
      size=len(list_of_files)/cpu
      size_reste=len(list_of_files)%cpu 
      list_of_file_index=[]
      if size!=0:
              for i in range(cpu):
                if i<size_reste:list_of_file_index.append([j for j in range(size*i,(i+1)*size)]+[int(len(list_of_files)-i-1)])
                else:list_of_file_index.append([j for j in range(size*i,(i+1)*size)])
              pool=mp.Pool(processes=cpu) 
              pool_result=pool.map(build_assignment, list_of_file_index)
              pool.close() 
      elif size==0:
              for i in range(len(list_of_files)):list_of_file_index.append([i])
              pool=mp.Pool(processes=len(list_of_file_index)) 
              pool_result=pool.map(build_assignment, list_of_file_index)
              pool.close()     
      for j in range(len(pool_result)):
        if  pool_result[j][1]>highest_score:highest_score=pool_result[j][1]
      for element in pool_result:
        archive=element[0]
        for peak in archive.keys():
          if not peak in archive_assignment_result.keys():archive_assignment_result[peak]={}
          for methyl in archive[peak].keys():
            if not methyl in archive_assignment_result[peak].keys():
              archive_assignment_result[peak][methyl]=archive[peak][methyl]
            if (methyl in archive_assignment_result[peak].keys()
              and archive[peak][methyl]>archive_assignment_result[peak][methyl]):
              archive_assignment_result[peak][methyl]=archive[peak][methyl]                
      ##################### Multi-processed assignment building #################### 
      list_of_files=os.listdir('./'+directory+'/run/temp/')
      flag_not_null=0
      for file_name in list_of_files:
              peak_id=re.search('^([0-9]+)#', file_name).group(1)
              if int(peak_id)==cluster_name:
                flag_not_null=1
                flag_print_result=1
                break
      if flag_not_null==1:
              for file_name in list_of_files:
                peak_id=re.search('^([0-9]+)#', file_name).group(1)
                if int(peak_id)!=cluster_name:os.remove('./'+directory+'/run/temp/'+str(file_name))
      List_of_assigned_peaks.append(cluster_name) 
      still_not_assigned.remove(cluster_name)
      if flag_not_null==0:archive_assignment_result=archive_assignment_result_previous                     
      list_of_peaks=archive_assignment_result.keys()
      list_of_peaks.sort()
      Time_end=time.datetime.now()      
      print 'Tc= '+str(P)+'; '+str(Time_end-Time_start).split('.')[0].split('.')[0], ' => ',str(round(100*len(list_of_peaks)/float(len(HMQC_peak_list)),1)),'%'
      
      report_follow=open('./'+directory+'/'+path_name+'/'+logFile_name, 'a')
      for peak in archive_assignment_result.keys():
          for methyl in archive_assignment_result[peak]:
            try:
              if re.search('\w(\d+)\w',str(methyl)).group(1)==re.search('(\d+)\w', peak).group(1):
                if archive_assignment_result[peak][methyl]==highest_score:report_follow.write('R* ')
                else:report_follow.write('R  ')
                break
            except AttributeError:pass
          report_follow.write(str(peak)+'  '+str(archive_assignment_result[peak])+'\n')     
      report_follow.close()
       
      turn=0
      break
    #if turn==4:cycler=1
    if turn==3:cycler=1
    else:no_crosspeaks_peak=[]
    outputing(highest_score,peak_matrix_Scoring)

########STATS FILE#######################################################################
try:file=open('./'+directory+'/Output/hmqc_iso.list', 'r')
except IOError:file=open('./'+directory+'/Output/hmqc.list', 'r')
histo_ambiguity={}
for peakline in file.readlines()[2:]:
  if peakline.count('NotAss')==0:
    number_of_methyls=peakline.count(':')
    if not str(number_of_methyls) in histo_ambiguity.keys():
      histo_ambiguity[str(number_of_methyls)]=1
    else:histo_ambiguity[str(number_of_methyls)]+=1
  else:
    number_of_methyls=peakline.count(',')+1
    if not str(number_of_methyls) in histo_ambiguity.keys():
      histo_ambiguity[str(number_of_methyls)]=1
    else:histo_ambiguity[str(number_of_methyls)]+=1
file.close()
tot_ass=0
tot_met=0
sorted_histo_keys=[int(i) for i in histo_ambiguity.keys()]
sorted_histo_keys.sort()
file_log=open('./'+directory+'/log', 'a')
file_log.write('\n\n\nSome statistics...\n')
for key in sorted_histo_keys:
  val=histo_ambiguity[str(key)]
  file_log.write(str(key)+': '+str(val)+'\n')
  tot_ass+=float(key)*float(val)
  tot_met+=float(val)
Time_end=time.datetime.now()
file_log.write('\nAssigned methyls: '+str(round(100*len(List_of_assigned_peaks)/float(number_of_methyls_from_pdb),1))+' %'
               '\nAssigned peaks: '+str(round(100*len(List_of_assigned_peaks)/float(len(HMQC_peak_list)),1))+' %'
               '\nAverage alternative assignments: '+str(round(tot_ass/tot_met,1))+' per peak'
               '\nCalculation time: '+str(Time_end-Time_start).split('.')[0].split('.')[0])
file_log.close
print 'Assignment complete'
