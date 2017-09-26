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

class obj:
  pass
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
                             (LCD2==1 and split[3]=='LEU' and split[2]=='CD2') or
                             (LCD1==1 and split[3]=='LEU' and split[2]=='CD1') or
                             (M==1 and split[3]=='MET' and  split[2]=='CE') or
                             (T==1 and split[3]=='THR' and  split[2]=='CG2') or  
                             (VCG1==1 and split[3]=='VAL' and split[2]=='CG1') or 
                             (VCG2==1 and split[3]=='VAL' and split[2]=='CG2')                      
      )): Atomlist.append(line)
  Methyl_list=[]
  distances_CHCH=[]
  if peak_geminal.sum()==0:flag_no_pairs=1
  else:flag_no_pairs=0
  
  for linei in Atomlist:
    spliti=linei.split()
    xi=float(spliti[6])
    yi=float(spliti[7])
    zi=float(spliti[8])
    for aa in convert:
      if spliti[3]==aa[0]: Res_name=aa[1]
    if spliti[2]!='H': Methyl_list.append(Res_name+spliti[5]+spliti[2])  
    for linef in Atomlist:
      splitf=linef.split()
      xf=float(splitf[6])
      yf=float(splitf[7])
      zf=float(splitf[8])
      if (xf!=xi or yf!=yi or zf!=zi):
        d=round(math.pow(math.pow(xf-xi,2)+math.pow(yf-yi,2)+math.pow(zf-zi,2),0.5),1)
        if (spliti[2]!='H' and splitf[2]!='H'): # CCH
            
            ##############Relax cutoff and lowCut if conformation changes#########
            Cutoff2=Cutoff
            lowCut2=lowCut
            flag_confchange=0
            for area in areas:
              if (int(spliti[5])>=int(area[0]) and int(spliti[5])<=int(area[1])):flag_confchange=1
            if flag_confchange==1:
              Cutoff2=Cutoff+2  
              lowCut2=lowCut+2
            #######################################################################
              
            if (d<=lowCut2 and spliti[5]!=splitf[5]):coef=1
            elif d>lowCut2 and d<=Cutoff2:coef=(Cutoff2-d)/float(Cutoff2-lowCut2)
            elif (flag_geminal==1 and spliti[5]==splitf[5]):
              if flag_no_pairs==1:coef=1.2
              else: coef=1
            elif (flag_geminal==0 and spliti[5]==splitf[5]):coef=0
            else:coef=0
            for line in convert:
              if spliti[3]==line[0]:Res_name_i=line[1]
              if splitf[3]==line[0]:Res_name_f=line[1]         
            if flag=='run':
              distances_CHCH.append([str(Res_name_i+spliti[5]+spliti[2]),
                                     str(Res_name_f+splitf[5]+splitf[2]),coef])
              #if coef!=0:methyl_file.write(str(distances_CHCH[-1])+'\n')
            elif flag=='histo':distances_CHCH.append([str(Res_name_i+spliti[5]+spliti[2]),
                                                      str(Res_name_f+splitf[5]+splitf[2]),d])
  #if flag=='run':methyl_file=open('./'+str(Time_start).split('.')[0]+'/Methyl_connectivity','w')
  distances_CHCH_new=[]
  for line in distances_CHCH:
    if (line[0][0]=='L' or line[0][0]=='V'):line[0]=line[0][:-1]
    if (line[1][0]=='L' or line[1][0]=='V'):line[1]=line[1][:-1]
    distances_CHCH_new.append(line)
    #if (flag=='run' and float(line[2])!=0):methyl_file.write(str(distances_CHCH_new[-1])+'\n')
  Methyl_list_new=[]
  for line in Methyl_list:
    if ((line[0]=='L' or line[0]=='V') and not line[:-1] in Methyl_list_new):Methyl_list_new.append(line[:-1])
    elif line[:-1] in Methyl_list_new:pass
    else:Methyl_list_new.append(line)  
  #if flag=='run':methyl_file.close()
  empty=[]
  matrix_CHCH,matrix_geminal=matrix_it(distances_CHCH_new, Methyl_list_new, 'pdb',empty,empty)
  matrix_CHCH_cluster=0

  return matrix_geminal, matrix_CHCH, Methyl_list_new, distances_CHCH_new
# ---------------------------------------------------------------------
# Extract 3d noesy peaks for each 2d peaks
#
def Extract_3DPeaks(object, peak_list_3d, peak_list_2d_CH, nuclei, tolerance):
  w1tol = tolerance[0]
  w2tol = tolerance[1]
  w3tol = tolerance[2]
  Ref_nuclei=nuclei[1]
  if Ref_nuclei=='13C': peak_list_2d_ref = peak_list_2d_CH
  DistNOE=[]
  tot_intensity=0
  Max_noe=0
  liste_NOE=[]
  for line_noe in peak_list_3d:
    split_noe=line_noe.split()
    tot_intensity+=float(split_noe[4])
    liste_NOE.append(float(split_noe[4]))
    if float(split_noe[4])>Max_noe:Max_noe=float(split_noe[4])
  mean_intensity_per_strip=tot_intensity/float(len(peak_list_2d_CH))
  mean_intensity=tot_intensity/float(len(peak_list_3d))
  array_NOE=np.array(liste_NOE)
  NOE_std=np.std(array_NOE)
  DistNOE.append(mean_intensity_per_strip)
  DistNOE.append(Max_noe)
  DistNOE.append(len(peak_list_3d))
  DistNOE.append(mean_intensity)
  DistNOE.append(NOE_std)
  file_log=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/log', 'a')
  file_log.write('Number of NOEs: '+str(len(peak_list_3d))+'\n')
  file_log.write('NOEs per strip: '+str(round(len(peak_list_3d)/float(len(peak_list_2d_CH)),2))+'\n')
  file_log.close()
  for line_hmqc in peak_list_2d_ref:
      split_hmqc=line_hmqc.split()
      w=float(split_hmqc[1])
      wh=float(split_hmqc[2])
      setattr(object, str(split_hmqc[0])+'_freq', (w,wh))
      noesy_peaks=[]
      if nuclei[0]=='13C': wnoe=w
 #     elif nuclei[0]=='1H': wnoe=wh
      for line_noe in peak_list_3d:
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
      setattr(object, str(split_hmqc[0])+'_3dpeaks', noesy_peaks)
  return DistNOE
# ----------------------------------------------------------------------------------------
# Build up the matrix
#
def matrix_it(links, element_list, flag, factors,geminal_mark):
  matrix=np.zeros((len(element_list),len(element_list)))
  if flag=='peak':
    matrix_scoring=np.zeros((len(element_list),len(element_list)))
    peak_geminal=np.zeros((len(element_list),len(element_list)))  
    for entry in geminal_mark:peak_geminal[entry[0],entry[1]]=1
  if flag=='pdb':matrix_geminal=np.zeros((len(element_list),len(element_list)))
  for line in links:
    i=element_list.index(line[0])
    j=element_list.index(line[1])
    if flag=='peak':matrix[i,i]=1    
    if flag=='pdb':
      matrix[i,j]=line[2]
      if i==j:matrix_geminal[i,j]=1
    if flag=='peak':
      chi=float(line[6])
      if ((float(line[2])>0 and line[7]==1) or factors[2]>0):
        if (line[4]==0):matrix[i,j]=math.pow(math.pow(1+float(line[5]),2)*float(line[7])*float(factors[0]/(20*chi*math.pow(float(line[3]),2))),0.5)
        elif (line[4]==1):matrix[i,j]=math.pow(math.pow(1+float(line[5]),2)*float(line[7])*float(factors[1]/(20*chi*math.pow(float(line[3]),2))),0.5)
        elif (line[4]==2):matrix[i,j]=math.pow(math.pow(1+float(line[5]),2)*float(line[7])*float(factors[2]/(20*chi*math.pow(float(line[3]),2))),0.5)
      else:matrix[i,j]=0
  if flag=='pdb':return matrix,matrix_geminal
  if flag=='peak':
    matrix_ri=np.zeros((len(element_list),len(element_list)))
    for line in links:
      i=element_list.index(line[0])
      j=element_list.index(line[1])
      chi=float(line[6])
      matrix_ri[i,j]=float(line[2])
      matrix_scoring[i,j]=float(line[7])*math.pow(math.pow(1+float(line[5]),2)/(20*chi*math.pow(float(line[3]),2)),0.5)
      if matrix_scoring[i,j]<matrix_scoring[j,i]:matrix_scoring[i,j]=matrix_scoring[j,i]
      else:matrix_scoring[j,i]=matrix_scoring[i,j]
    matrix_scoring=matrix_scoring*matrix_ri
    return matrix,matrix_scoring,peak_geminal
# ---------------------------------------------------------------------------
#
def noe2matrix(result,factors,flag_cchlist):
    Noe_correlations_clustering=[]
    Noe_correlations_scoring=[]
    HMQC_peak_list=[]
    if flag_cchlist==0:noe_network=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/Peak_connectivity', 'w')
    overlap_list=[]
    overlap_topo=[]
    RiskyPeak=[]
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
      if len(HMQC[i].split())==5:
        geminal_i_list=HMQC[i].split()[4].split(';')
        for geminals in geminal_i_list[1:]:
          for j in range(len(HMQC)):
            if HMQC[j].split()[0]==geminals:geminal_mark.append((i,j,len(geminal_i_list[1:])))
        
    for linei in HMQC:
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
             total_height_i+=float(peak_split[4])
             if float(peak_split[4])>imax:imax=float(peak_split[4])
           if total_height_i>0:rimax=float(imax/total_height_i)
           height_i2mean=round(float(total_height_i)/DistNOE[0], 2)
       except AttributeError: pass
       if not str(linei_split[0]) in biblio_crosspeaks.keys():biblio_crosspeaks[str(linei_split[0])]={}     
       if noe_peaki_list:
         counter=0     
         for noe_peaki in noe_peaki_list:
             noe_peaki_split=noe_peaki.split()
             ri=round(float(noe_peaki_split[4])/total_height_i, 2)
             ri2max=round(float(noe_peaki_split[4])/DistNOE[1],3)    
             peak_neighbors_clustering=[]
             peak_neighbors_scoring=[]
             shared_peak_neighbors=0
             for linef in HMQC:
                 linef_split=linef.split()
                 if ((float(noe_peaki_split[1])>(float(linef_split[1])-ppm_range[0])) and 
                     (float(noe_peaki_split[1])<(float(linef_split[1])+ppm_range[0]))):
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
                           if ((float(noe_peakii_split[1])>(float(noe_peakf_split[1])-ppm_range[0])) and 
                               (float(noe_peakii_split[1])<(float(noe_peakf_split[1])+ppm_range[0]))
                                ):shared_peak_neighbors+=1
                         if ((float(noe_peakf_split[1])>(float(linei_split[1])-ppm_range[0])) and 
                             (float(noe_peakf_split[1])<(float(linei_split[1])+ppm_range[0])) and
                             (not str(linef_split[0]) in peak_neighbors_clustering)):
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
               else:a=2
               d1=float(peak_neighbors_clustering[i][2][0])
               d2=float(peak_neighbors_clustering[i][2][1])
               chi=round(math.pow(math.pow(d1,2)+math.pow(d2,2),0.5),4)
               if chi<0.05:chi=0.05
               if peak_neighbors_clustering[i][3]==1:N=n
               else:N=m
               
               sym=peak_neighbors_clustering[i][3]
               ################## CONNECTION UPGRADING #########################
               if (n==0 and int(peak_neighbors_clustering[i][1])>2):sym=1
               ###########################################################
               Noe_correlations_clustering.append([linei_split[0],peak_neighbors_clustering[i][0],
                                                  ri,N,a,peak_neighbors_clustering[i][1],chi,sym])
               if flag_cchlist==0:
                 score=round(math.pow(math.pow(1+float(peak_neighbors_clustering[i][1]),2)*float(sym)*float(factors[0]/(20*chi*math.pow(float(N),2))),0.5),3)
                 noe_network.write(str(Noe_correlations_clustering[-1])+'  '+str(score)+'\n') 
                 
             n=len(peak_neighbors_scoring)
             for i in range(len(peak_neighbors_scoring)):
               if (not peak_neighbors_scoring[i][0] in overlap_list and not linei_split[0] in overlap_list):a=0
               else: a=1
               d1=peak_neighbors_scoring[i][2][0]
               #chi=round(math.pow(math.pow(d1,2),0.5),4)
               #if chi<0.05:chi=0.05
               chi=0.05
               if (a==1 or n>=2):height=0
               else:
                 a,n=0,1
                 height=float(noe_peaki_split[4])
                 if flag_cchlist==2:height=float(ri/rimax)
               Noe_correlations_scoring.append([linei_split[0],peak_neighbors_scoring[i][0],height,n,a,'0',chi,'1'])

       if height_i2mean<0.1:RiskyPeak.append(1)
       else:RiskyPeak.append(0)
    noe_peak_matrix_clustering,noe_peak_matrix_scoring,peak_geminal=matrix_it(Noe_correlations_clustering, HMQC_peak_list, 'peak',factors,geminal_mark)
    
    C_matrix=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/C_matrix','w')
    np.savetxt(C_matrix,noe_peak_matrix_clustering,fmt='%.3f')
    C_matrix.close()
    
    if flag_cchlist==1:noe_peak_matrix_clustering,noe_peak_matrix_scoring,peak_geminal=matrix_it(Noe_correlations_scoring, HMQC_peak_list, 'peak',factors,geminal_mark)
    if flag_cchlist==2:noe_peak_matrix_clustering,noe_peak_matrix_scoring,peak_geminal=matrix_it(Noe_correlations_clustering, HMQC_peak_list, 'peak',factors,geminal_mark)
    if flag_cchlist==0:noe_network.close()
    return noe_peak_matrix_scoring, noe_peak_matrix_clustering, HMQC_peak_list,overlap_topo,RiskyPeak,biblio_crosspeaks,assignment_locked_peak,peak_geminal

# ---------------------------------------------------------------------------
#
def build_assignment(index):
              highest_score_inloop=highest_score
              flag_ram=0
              assignment_archive=[]
              assignment_collection={}
              index_files=0
              Time=time.datetime.now()
              report_follow=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/'+str(path_name)+'/'+str(logFile_name), 'a')
              report_follow.write(str(Time).split('.')[0]+'\n')
              report_follow.close()              
              file_log=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/log', 'a')
              file_size=int(1000+float(file_size_max*size/float(5+size)))              
              file_log.write('file_size: '+str(file_size)+'\n')
              file_log.close()
              for iii in index:             
                free_RAM=round(float(psutil.virtual_memory().available)/float(1024*1024*1024),1)              
                file_log=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/log', 'a')
                file_log.write('free_RAM: '+str(free_RAM)+' GB'+'\n')          
                file_log.write(str(index.index(iii)+1)+'/'+str(len(index))+'\n')
                file_log.close()            
                file_name=list_of_files[iii]
                file=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/run/temp/'+str(file_name), 'r')
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
                    testing=matrix_assignment*matrix_noe
                    #####Geminal scaling factor####################
                    if (flag_geminal==1 and matrix_noe_geminal.sum()>0):
                      matrix_assignment_geminal=new_metrics[1][:,index_methyls][index_methyls,:]
                      testing_geminal=matrix_noe_geminal*matrix_assignment_geminal
                      tot=round(testing.sum()*(1+0.2*(testing_geminal.sum()/float(matrix_noe_geminal.sum()))),3)
                    #####Geminal scaling factor
                    else:tot=round(testing.sum(),3)
                    if tot>highest_score_inloop:highest_score_inloop=tot                
                    cutoff=highest_score_inloop-FILTER
                    if tot>=cutoff:assignment_archive.append([index_methyls,tot])
                                       
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
                    report_follow=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/'+str(path_name)+'/'+str(logFile_name), 'a')
                    report_follow.write('----total#'+str(tot)+', deleted# '+str(deleted)+'\n')
                    report_follow.close()
                    assignment_archive=list(assignment_archive_sorted)
                    assignment_archive_sorted=[]
                    number_of_file=int((len(assignment_archive)-1))/int(file_size)
                    number_of_file_rest=(len(assignment_archive)-1)%int(file_size)          
                    for i in range(number_of_file):
                        assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[i*file_size+1:(i+1)*file_size+1]
                        file=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                        cPickle.dump(assignment_archive_to_save,file,-1)
                        file.close()
                        index_files+=1
                    if number_of_file_rest!=0:
                        assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[number_of_file*file_size+1:]
                        file=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
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
                  report_follow=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/'+str(path_name)+'/'+str(logFile_name), 'a')
                  report_follow.write('----total#'+str(tot)+', deleted# '+str(deleted)+'\n')
                  report_follow.close()
                  assignment_archive=list(assignment_archive_sorted)
                  assignment_archive_sorted=[]                 
                  if (len(assignment_archive)-1)>file_size:
                    number_of_file=(len(assignment_archive)-1)/int(file_size)
                    number_of_file_rest=(len(assignment_archive)-1)%int(file_size)
                    for i in range(number_of_file):
                      assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[i*file_size+1:(i+1)*file_size+1]
                      file=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                      cPickle.dump(assignment_archive_to_save,file,-1)
                      file.close()
                      index_files+=1
                    if number_of_file_rest!=0:
                      assignment_archive_to_save=[assignment_archive[0]]+assignment_archive[number_of_file*file_size+1:]
                      file=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                      cPickle.dump(assignment_archive_to_save,file,-1)
                      file.close()
                      index_files+=1
                    assignment_archive=[]
                  else:
                    file=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/run/temp/'+str(cluster_name)+'#'+str(index[0])+'#'+str(index_files), 'w')
                    cPickle.dump(assignment_archive,file,-1)
                    file.close()                    
                    assignment_archive=[]
              max=highest_score_inloop
              return assignment_collection,max

def build_assignment_peak(index):
        count=0
        global q
        global total_Assignments
        assignment_archive=[]        
        for index_number in index:
          assignment_old=assignment_archive_old[index_number]
          #locked_noe=[]
          #locked_methyl=[]
          #for i in range(len(assignment_locked_peak)):
          #  if (assignment_locked_peak[i] in new_peaks and 
          #     not assignment_locked_peak[i] in assignment_old[0]):
          #    locked_noe=locked_noe+[assignment_locked_peak[i]]
          #    locked_methyl=locked_methyl+[assignment_locked_methyl[i]]       
          if total_Assignments==[]:total_Assignments=[[],[]]
          #index_matrix_noe=locked_noe+total_Assignments[0]+assignment_old[0]
          index_matrix_noe=total_Assignments[0]+assignment_old[0]
          matrix_noe=peak_matrix_Scoring[:,index_matrix_noe][index_matrix_noe,:]
          matrix_noe_geminal=peak_geminal[:,index_matrix_noe][index_matrix_noe,:]
          
          for assignment in total_Assignments[1:]:          
            #index_matrix_methyls=locked_methyl+assignment+assignment_old[1]
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
            try:
              flag_stop=0
              for i in range(len(index_matrix_noe)):
                peak=HMQC_peak_list[int(index_matrix_noe[i])]
                methyl=metrics[2][int(index_matrix_methyls[i])]
                if peak in archive_assignment_result.keys():
                  if not methyl in archive_assignment_result[peak].keys():
                    flag_stop=1
                    break
                  else:pass
                else:pass
              if flag_stop==1:continue
            except NameError:pass
            #####OUT: Test allowed assignment for all peaks #######
            
            testing=matrix_assignment*matrix_noe
            tot=testing.sum()
            
            #####Geminal scaling factor
            if (flag_geminal==1 and matrix_noe_geminal.sum()>0):
              matrix_assignment_geminal=new_metrics[1][:,index_matrix_methyls][index_matrix_methyls,:]
              testing_geminal=matrix_noe_geminal*matrix_assignment_geminal
              tot=round(tot*(1+0.2*(testing_geminal.sum()/float(matrix_noe_geminal.sum()))),3)
            #####Geminal scaling factor
            else:pass
            assignment_archive.append([index_matrix_noe,index_matrix_methyls,tot])
            count=count+1   
        if len(assignment_archive)>1:          
          assignment_archive_brut=list(assignment_archive)
          scorelist=np.zeros((1,len(assignment_archive_brut)))
          for i in range(len(assignment_archive_brut)):scorelist[0,i]=assignment_archive_brut[i][2]
          max=scorelist[0,np.argmax(scorelist)]
          if max<highest_score_small:max=highest_score_small
          number_peaks=len(index_matrix_noe)
          if number_peaks<=2:q=0           
          elif number_peaks==3:q=0.5
          elif number_peaks==4:q=0.75
          elif number_peaks>=5:q=0.8
          #elif number_peaks==6:q=0.85
          #elif number_peaks>=7:q=0.9
          cutoff=max*q           
          #cutoff=max-FILTER_local
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
    file_log=open('./'+str(Time_start).split('.')[0].split('.')[0]+'/log', 'a')
    file_log.write('----> Building up OUTPUT files...'+'\n')
    file_log.close()

    if P=='iso':
      hmqc_result_list=open('./'+str(Time_start).split('.')[0]+'/Output/hmqc_iso.list', 'w')
      cch_result_list=open('./'+str(Time_start).split('.')[0]+'/Output/cch_iso.list', 'w')
    else:
      hmqc_result_list=open('./'+str(Time_start).split('.')[0]+'/Output/hmqc.list', 'w')
      cch_result_list=open('./'+str(Time_start).split('.')[0]+'/Output/cch.list', 'w')
          
    hmqc_result_list.write('Assignment\tw1\tw2\tNote\n\n')
    cch_result_list.write('Assignment\tw1\tw2\tw3\tNote\n\n')

############## TO IMPROVE ###############
    peaks=[]
    methyls=[]
    list_peak=archive_assignment_result.keys()
    list_peak.sort()   
    for peak_name in list_peak:
      peaks.append(HMQC_peak_list.index(peak_name))      
      #print 'peaks',peaks
      score=[]
      for i in range(len(archive_assignment_result[peak_name].keys())):
        methyl_name=archive_assignment_result[peak_name].keys()[i]
        score.append((archive_assignment_result[peak_name][methyl_name],methyl_name))
      #print 'score',score
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
    factors=(1,1,0.25)
    sym_coef=0.1
    peak_matrix_Scoring2,peak_matrix_clustering2,glup,glop,glip,glep,glyp,glap=noe2matrix(result,factors,3)
    matrix_peaks_CS=peak_matrix_clustering2[peaks,:][:,peaks] ###peak_matrix_clustering2 instead of peak_matrix_clustering if factor is changed
    matrix_peaks_Cluster=peak_matrix_clustering[peaks,:][:,peaks]
    matrix_peaks=peak_matrix_Scoring2[peaks,:][:,peaks]
    matrix_score=matrix_assignment*matrix_peaks
    correct=0
    total=0
    metrics2histo=distances(pdb,100,100,'histo')
    distance2histo,glop=matrix_it(metrics2histo[3], metrics[2], 'pdb',empty,empty)
    matrix_methyls=distance2histo[methyls,:][:,methyls]    
    for i in range(len(peaks)):
          peak=HMQC_peak_list[int(peaks[i])]
          methyl=metrics[2][int(methyls[i])]  
          if re.search('\w(\d+)\w',str(methyl)).group(1)==re.search('(\d+)\w', peak).group(1):correct+=1
          total+=1
          freq=getattr(result, str(peak)+'_freq')
          wC=freq[0]
          wH=freq[1]
          if methyl[0]=='I':Name=str(methyl)+'-HD1'
          elif methyl[0]=='M':Name=str(methyl)+'-HE'
          elif methyl[0]=='L':Name=str(methyl)+'*-HD*'
          elif methyl[0]=='V':Name=str(methyl)+'*-HG*'
          elif methyl[0]=='A':Name=str(methyl)+'-HB'
          elif methyl[0]=='T':Name=str(methyl)+'-HG2'
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
                        completeness[0,NOE_index]==0):
                      completeness[0,NOE_index]=1
                      break
          if completeness.size>0:NOE_completeness=round(completeness.sum()/float(completeness.size),2)
          else: NOE_completeness='no NOE'
          hmqc_result_list.write('  '+str(round(matrix_peaks_CS[i,:].sum(),2))+
                                 '  '+str(NOE_completeness)+
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
          if not methyl_index in methyls:allowed_assignment.append(metrics[2][methyl_index])
          elif (('V' or 'L') in code_methyl_type and
                methyls.count(methyl_index)==1 and
                flag_geminal!=2 and
                (methyl_index in Leu_methyls or methyl_index in Val_methyls)):allowed_assignment.append(metrics[2][methyl_index])
        freq=getattr(result, str(peak)+'_freq')
        wC=freq[0]
        wH=freq[1]
        Name2=HMQC[i].split()[0]
        hmqc_result_list.write(Name2+'\t'+str(wC)+'\t'+str(wH))
        hmqc_result_list.write('\tNotAss'+str(allowed_assignment)+'\n')        
    hmqc_result_list.close()
    
    if P=='iso':
      file=open('./'+str(Time_start).split('.')[0]+'/Output/cch_iso.list', 'r')
      noelist=file.readlines()
      file.close()
      file2=open('./'+str(Time_start).split('.')[0]+'/Output/hmqc_iso.list', 'r')
      hmqclist=file2.readlines()
      file2.close()
      results=open('./'+str(Time_start).split('.')[0]+'/Output/mapping_iso.pml','w')
    else:
      file=open('./'+str(Time_start).split('.')[0]+'/Output/cch.list', 'r')
      noelist=file.readlines()
      file.close()
      file2=open('./'+str(Time_start).split('.')[0]+'/Output/hmqc.list', 'r')
      hmqclist=file2.readlines()
      file2.close()
      results=open('./'+str(Time_start).split('.')[0]+'/Output/mapping.pml','w')
                
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
    M=2
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
        if split[3][0]!='N':
            split0=split[0].split('-')
            res1i=int(re.search('[A-Z]([0-9]+)[A-Z]', split0[0]).group(1))
            name1i=re.search('[0-9]+(\w+)', split0[0]).group(1)
            try:
              score=float(split[4])
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
os.makedirs('./'+str(Time_start).split('.')[0])
os.makedirs('./'+str(Time_start).split('.')[0]+'/Input')
os.makedirs('./'+str(Time_start).split('.')[0]+'/Output')
shutil.copy('./'+str(sys.argv[1]),'./'+str(Time_start).split('.')[0]+'/Input/'+str(sys.argv[1]))
shutil.copy('./'+str(sys.argv[0]),'./'+str(Time_start).split('.')[0]+'/Input/'+str(sys.argv[0]))

file_log=open('./'+str(Time_start).split('.')[0]+'/log', 'w')
file_log.write(str(Time_start).split('.')[0]+'\n')

file_test=open('./'+str(Time_start).split('.')[0]+'/test', 'w')
file_test.close()

file=open(str(sys.argv[1]),'r')
start_file=file.readlines()
file.close()

CCH_name=start_file[4][:-1]
HMQC_name=start_file[2][:-1]
pdb_name=start_file[6][:-1]
Seq_name=start_file[8][:-1]
labeling=start_file[10]
flag_geminal=float(start_file[13][:-1])
ppm_range=(float(start_file[15].split()[0]),float(start_file[15].split()[1]),float(start_file[15].split()[2]))
cutoff_factor=float(start_file[17][:-1])
lowCut,max_distance=float(start_file[19].split()[0]),float(start_file[19].split()[1])
score_tol_end=str(start_file[21])

FILTER_scale=cutoff_factor

areas=[]
line=str(start_file[23])
line=line.replace('\n','')
for area in line.split(';'):
  try:areas.append((area.split('-')[0],area.split('-')[1]))
  except IndexError:pass
file_log.write('set areas: '+str(areas)+'\n')

shutil.copy('./'+str(CCH_name),'./'+str(Time_start).split('.')[0]+'/Input/'+str(CCH_name))
shutil.copy('./'+str(HMQC_name),'./'+str(Time_start).split('.')[0]+'/Input/'+str(HMQC_name))
shutil.copy('./'+str(pdb_name),'./'+str(Time_start).split('.')[0]+'/Input/'+str(pdb_name))
shutil.copy('./'+str(Seq_name),'./'+str(Time_start).split('.')[0]+'/Input/'+str(Seq_name))

low_P,P_high=1,40
file_log.write('ppm_range: '+str(ppm_range)+'\n')
cpu=mp.cpu_count()
Tot_RAM=round(float(psutil.virtual_memory().total)/float(1024*1024*1024),1)
free_RAM_start=round(float(psutil.virtual_memory().available)/float(1024*1024*1024),1)
file_log.write('cpu: '+str(cpu)+' and total RAM: '+str(Tot_RAM)+' GB'+'\n')

				#purgeMAX=0.5*free_RAM_start-0.2*free_RAM_start
              	# MAXfile= 0.72 / 0.06 = 12
              	# MAXfile/cpu=12/4=3 > 1.5 for security
				#archive_limit=50000*Tot_RAM
file_size_max=int(1000000*(0.5*free_RAM_start-0.2*free_RAM_start)/float(2*cpu*0.06))
file_size=10000

CCH_file=open('./'+str(CCH_name),'r')
CCH=CCH_file.readlines()
CCH_file.close()
HMQC_file=open('./'+str(HMQC_name),'r')
HMQC=HMQC_file.readlines()
HMQC_file.close()
pdb_file=open('./'+str(pdb_name),'r')
pdb=pdb_file.readlines()
pdb_file.close()
Seq_file=open('./'+str(Seq_name),'r')
SEQ=Seq_file.readlines()
Seq_file.close()
nuclei=['13C','13C','1H']

DistNOE=Extract_3DPeaks(result,CCH,HMQC,nuclei,ppm_range)

factors=(1,1,0)
sym_coef=0              
peak_matrix_Scoring,peak_matrix_clustering,HMQC_peak_list,overlap_topo,RiskyPeak,biblio_crosspeaks,assignment_locked_peak,peak_geminal=noe2matrix(result,factors,0)

#print 'average_PER10', np.percentile(peak_matrix_Scoring[np.where(peak_matrix_Scoring>0)],10)
#if len(HMQC_peak_list)<=100:scaler=float(1)
#else:scaler=len(HMQC_peak_list)/float(100)
scaler=len(HMQC_peak_list)/float(100)
#print 'median',np.median(peak_matrix_Scoring.sum(axis=0))
#print 'mean',np.mean(peak_matrix_Scoring.sum(axis=0))
FILTER_start=np.mean(peak_matrix_Scoring.sum(axis=0))#*len(HMQC_peak_list)

file_log.write('Mean maximal peak score and scaler: '+str(FILTER_start)+'; '+str(scaler)+'\n')
calcul=peak_matrix_clustering[np.where(peak_matrix_clustering>0)]
main_confident=np.mean(calcul)
file_log.write('Average peak connection confidence score: '+str(main_confident)+'\n')

metrics=distances(pdb, max_distance,lowCut,'run')

methyl_to_add=[]
for line in SEQ:
  resID=re.search('\w(\d+)', line).group(1)
  resTYPE=line[0]
  flag=0
  for methyl in metrics[2]:
    methylID=re.search('\w(\d+)\w', methyl).group(1)
    if methylID==resID:
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
new_metrics=matrix_it(metrics[3], metrics[2], 'pdb',empty,empty)

methyl_file=open('./'+str(Time_start).split('.')[0]+'/Methyl_connectivity','w')
for line in metrics[3]:
  if line[2]!=0:methyl_file.write(str(line)+'\n')
methyl_file.close()

pdb_matrix=open('./'+str(Time_start).split('.')[0]+'/pdb_matrix','w')
np.savetxt(pdb_matrix,new_metrics[0],fmt='%.3f')
pdb_matrix.close()

file_log.write('###################################'+'\n')
A,I,L,M,T,V=0,0,0,0,0,0
for i in range(len(metrics[2])):
  if metrics[2][i][0]=='A':A+=1
  elif metrics[2][i][0]=='I':I+=1
  elif metrics[2][i][0]=='L':L+=1
  elif metrics[2][i][0]=='M':M+=1
  elif metrics[2][i][0]=='T':T+=1
  elif metrics[2][i][0]=='V':V+=1
file_log.write('Number of methyls: '+str(A+I+2*L+M+T+2*V)+'\n')
file_log.write('A: '+str(A)+'\nI: '+str(I)+'\nL: '+str(2*L)+'\nM: '+str(M)+'\nT: '+str(T)+'\nV: '+str(2*V)+'\n')

file_log.write('###################################'+'\n')
file_log.write('Number of peaks: '+str(len(HMQC))+'\n')
file_log.write('If stated for unique methyl type only:'+'\n')
a,il,l,m,t,v=0,0,0,0,0,0
for i in range(len(HMQC)):
  sum=0
  for j in HMQC[i].split()[3]:sum+=1
  if sum==1:
    if 'A' in HMQC[i].split()[3]:a+=1
    elif 'I' in HMQC[i].split()[3]:il+=1
    elif 'L' in HMQC[i].split()[3]:l+=1
    elif 'M' in HMQC[i].split()[3]:m+=1
    elif 'T' in HMQC[i].split()[3]:t+=1
    elif 'V' in HMQC[i].split()[3]:v+=1
file_log.write('A: '+str(a)+'\nI: '+str(il)+'\nL: '+str(l)+'\nM: '+str(m)+'\nT: '+str(t)+'\nV: '+str(v)+'\n')
file_log.write('overlap_topo:'+'\n')
for i in range(len(overlap_topo)):
  if len(overlap_topo[i])>0:file_log.write(str(HMQC_peak_list[i])+' '+str(HMQC_peak_list[int(overlap_topo[i][0])])+'\n')
file_log.write('assignment_locked_peak: '+str(len(assignment_locked_peak))+'\n')
for i in assignment_locked_peak:
  file_log.write(str(HMQC_peak_list[i])+'\n')
  
if (a>A or il>I or l/2>L or v/2>V or t>T or m>M):
  file_log.write('More peaks than methyls, please review the hmqc list'+'\n')
  sys.exit()
else:pass
file_log.close()

report=open('./'+str(Time_start).split('.')[0]+'/Peak_list', 'w')
for i in range(len(HMQC_peak_list)): report.write(str(i)+' '+str(HMQC_peak_list[i])+'\n')
report.close()

report=open('./'+str(Time_start).split('.')[0]+'/Methyl_list', 'w')
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
#############Analysis of network density matrix###########################################
N=peak_matrix_clustering.shape[0]
#matrix2=peak_matrix_clustering
matrix2=np.dot(peak_matrix_clustering,peak_matrix_clustering)
for i in range(N):matrix2[i,i]=0
P_high=round(np.amax(matrix2),0)


##########################################################################################
################# Build up small redundant clusters ######################################
#
if not os.path.exists('./'+str(Time_start).split('.')[0]+'/run/Local'):os.makedirs('./'+str(Time_start).split('.')[0]+'/run/Local')
if not os.path.exists('./'+str(Time_start).split('.')[0]+'/run/Local/temp'):os.makedirs('./'+str(Time_start).split('.')[0]+'/run/Local/temp')
if not os.path.exists('./'+str(Time_start).split('.')[0]+'/run/Local/archive'):os.makedirs('./'+str(Time_start).split('.')[0]+'/run/Local/archive')
P=P_high
selected_peaks_clusters=obj()
N = peak_matrix_clustering.shape[0]
#matrix2=peak_matrix_clustering
matrix2=np.dot(peak_matrix_clustering,peak_matrix_clustering)
for i in range(N):
      neighbors_with_sharing=np.argwhere(matrix2[i,:]>P) # 2 => 1 shared neighbors, 3=>2, etc.
      if not i in neighbors_with_sharing:
        neighbors_with_sharing=[neighbors_with_sharing[j,0] for j in range(neighbors_with_sharing.size)]
        neighbors_with_sharing.append(i)
        neighbors_with_sharing=np.array(neighbors_with_sharing).reshape((len(neighbors_with_sharing),1))
      setattr(selected_peaks_clusters, str(i), neighbors_with_sharing)
liste=dir(selected_peaks_clusters)
liste.remove('__doc__')
liste.remove('__module__')
liste.sort(key=float)

list_peak=dir(selected_peaks_clusters)
list_peak.remove('__doc__')
list_peak.remove('__module__')
list_peak.sort(key=float)

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
        code_methyl_type=HMQC[i].split()[3]
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
        tot_res=testing.sum(axis=0)
        archive=[]
        archive.append(index_matrix_noe)
        archive.append(index_matrix_methyls)
        archive.append(tot)
        assignment_archive.append(archive)
        count=count+1             

      report=open('./'+str(Time_start).split('.')[0]+'/run/Local/'+str(P)+'_cp#'+str(name_peak), 'w')
           
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
      #report.write('#cutoff= '+str(cutoff)+' #q= '+str(q)+'\n')
      for i in sort[:]:
        assignment=assignment_archive[i]
        report.write('#'+str(int(score_index[1,i]))+': '+str(assignment[2])+' ### ')
        for j in range(len(assignment[0])):
          peak=HMQC_peak_list[assignment[0][j]]
          methyl=metrics[2][assignment[1][j]]
          report.write(str(peak)+'=='+str(methyl)+' # ')
        report.write('\n')
      report.close()
      file=open('./'+str(Time_start).split('.')[0]+'/run/Local/archive/archive_'+str(name_peak), 'w')
      cPickle.dump(assignment_archive,file,-1)
      file.close()
      assignment_archive=[]      
######################## iterative RUNs ####################
#
P_list1=np.arange(10, P_high, 1)
P_list2=np.arange(2, 10, 0.01)
P_list=np.concatenate((P_list2,P_list1))
P_list=P_list[::-1]
selected_peaks_clusters=obj()
peak_clusters=obj()
N=peak_matrix_clustering.shape[0]
#matrix2=peak_matrix_clustering
matrix2=np.dot(peak_matrix_clustering,peak_matrix_clustering)
List_of_assigned_peaks=[]
already_called_cluster=[]
postpone_peaks=[]
delayed_clusters=[]
flag_score_equation=0
archive_assignment_result={}
completeness=np.zeros((len(HMQC_peak_list),2))
for i in range(len(HMQC_peak_list)):completeness[i,1]=len(biblio_crosspeaks[HMQC_peak_list[i]].keys())
for P in P_list:
    if P>100:cluster_size=7
    elif P>50:cluster_size=6
    elif P>10:cluster_size=5
    elif P>5:cluster_size=4
    elif P>=2:cluster_size=3
    #else:cluster_size=2 
    if P>0:        
      flag_new_peaks=0
      for i in range(N):
        neighbors_with_sharing=np.argwhere(matrix2[i,:]>=P)
        if not i in neighbors_with_sharing:
          neighbors_with_sharing=[neighbors_with_sharing[j,0] for j in range(neighbors_with_sharing.size)]
          neighbors_with_sharing.append(i)
          neighbors_with_sharing=np.array(neighbors_with_sharing).reshape((len(neighbors_with_sharing),1))
        setattr(selected_peaks_clusters, str(i), neighbors_with_sharing)  
      ############################# NEW CHECKING ######################################
      list_peak=dir(selected_peaks_clusters)
      list_peak.remove('__doc__')
      list_peak.remove('__module__')
      list_peak.sort(key=float)
      for name_peak_index in range(len(list_peak)):      
        highest_score_small=0  
        if int(list_peak[name_peak_index]) in already_called_cluster:continue
        if int(list_peak[name_peak_index]) in delayed_clusters:continue
        #######IN: Skip peak assignment if overlapping already assigned peak #########
        flag_continue=0
        for overpeak in overlap_topo[int(list_peak[name_peak_index])]:
              if (overpeak in already_called_cluster or overpeak in assignment_locked_peak):flag_continue=1
        if flag_continue==1:continue
        #######OUT: Skip peak assignment if overlapping already assigned peak ########       
        set_of_peaks=getattr(selected_peaks_clusters, str(list_peak[name_peak_index]))
        file=open('./'+str(Time_start).split('.')[0]+'/run/Local/archive/archive_'+str(list_peak[name_peak_index]), 'r')
        assignment_archive_old=cPickle.load(file)
        file.close()
        
        if len(assignment_archive_old)>100000:
          already_called_cluster.append(int(list_peak[name_peak_index]))
          continue
        
        try:                 
          new_peaks=[]
          for peak in set_of_peaks:
            if not peak in assignment_archive_old[0][0]:
              flag_continue=0
              for overpeak in overlap_topo[peak]:
                if (overpeak in assignment_archive_old[0][0] or overpeak in new_peaks):flag_continue=1
              if flag_continue==0:new_peaks.append(peak[0])          
          setattr(result, str(list_peak[name_peak_index])+'_new_peaks',new_peaks)
        except IndexError:continue
        
        if (len(new_peaks)==0):continue
        file_log=open('./'+str(Time_start).split('.')[0]+'/log', 'a')
        file_log.write('Local: Tc= '+str(P)+' ; Clustering around peak # '+str(name_peak_index)+'\n')
        file_log.close()
        
        if flag_new_peaks==1:pass
        elif set_of_peaks.size>=cluster_size:flag_new_peaks=1
        else:flag_new_peaks=0
        
        sorting=[(peak_matrix_clustering[peak,int(list_peak[name_peak_index])],peak) for peak in new_peaks]
        sorting.sort()        
        sorting=sorting[::-1]
        
        for line in sorting:
          peak=line[1]
          file=open('./'+str(Time_start).split('.')[0]+'/run/Local/archive/archive_'+str(list_peak[name_peak_index]), 'r')
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
                file=open('./'+str(Time_start).split('.')[0]+'/run/Local/temp/temp_'+str(number_of_files), 'w')
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
                file=open('./'+str(Time_start).split('.')[0]+'/run/Local/temp/temp_'+str(number_of_files), 'w')
                cPickle.dump(assignments_table,file,-1)
                file.close()
                assignments_table=[]
                number_of_files+=1
          ##################### Multi-processed assignment building ####################  
          number_peaks=len(assignment_archive_old[0][0])+1
          if number_peaks<=2:q=0           
          elif number_peaks==3:q=0.5
          elif number_peaks==4:q=0.75
          elif number_peaks>=5:q=0.8
          #elif number_peaks==6:q=0.85
          #elif number_peaks>=7:q=0.9    
          cutoff=highest_score_small*q
          #cutoff=highest_score_small-FILTER_local
          tot=0
          deleted=0
        
          assignment_archive=[]
          list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/Local/temp/')
          for i in range(len(list_of_files)):
              file=open('./'+str(Time_start).split('.')[0]+'/run/Local/temp/'+str(list_of_files[i]), 'r')
              element=cPickle.load(file)
              file.close()
              for j in range(len(element)):
                for k in range(len(element[j])):
                  tot+=1
                  if element[j][k][2]>=cutoff:assignment_archive.append(element[j][k])
              os.remove('./'+str(Time_start).split('.')[0]+'/run/Local/temp/'+str(list_of_files[i]))  
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
        
          report=open('./'+str(Time_start).split('.')[0]+'/run/Local/'+str(P)+'_cp#'+str(list_peak[name_peak_index]), 'w')        
          score_index=np.zeros((2,len(assignment_archive)))   
          for i in range(len(assignment_archive)):
            assignment=assignment_archive[i]
            score_index[0,i]=assignment[2]
            score_index[1,i]=str(i)      
          score_sorting=score_index[0,:]
          sort=np.argsort(score_sorting, axis=None)
          sort=sort[::-1] #mirror flip on both dimensions
          #report.write('#cutoff= '+str(cutoff)+' #q= '+str(q)+'\n')
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
          file=open('./'+str(Time_start).split('.')[0]+'/run/Local/archive/archive_'+str(list_peak[name_peak_index]), 'w')
          cPickle.dump(assignment_archive,file,-1)
          file.close()
          assignment_archive=[] 
#####################################################################################################      
###################### Mix all clusters together #################################
    else:pass
    if P>0:
      if (flag_new_peaks==1 or P==2 or P==10 or P==5 or P==50 or P==100):pass
      else:continue
      logFile_name='logFile_'+str(P)
      path_name='/run/Global'
      if not os.path.exists('./'+str(Time_start).split('.')[0]+'/run/Global'):os.makedirs('./'+str(Time_start).split('.')[0]+'/run/Global')
      if not os.path.exists('./'+str(Time_start).split('.')[0]+'/run/temp'):
        os.makedirs('./'+str(Time_start).split('.')[0]+'/run/temp')
        assignment_archive_FINAL=[assignment_locked_peak,[assignment_locked_methyl,]]
        file=open('./'+str(Time_start).split('.')[0]+'/run/temp/100000000#', 'w')
        cPickle.dump(assignment_archive_FINAL,file,-1)
        file.close()
        highest_score=0
        flag_alone=0
      
      if P>0: 
        list_peak=dir(selected_peaks_clusters)
        list_peak.remove('__doc__')
        list_peak.remove('__module__')
        list_peak.sort(key=float)
        counter=0
        sort_index=np.zeros((2,len(list_peak)))          
        for i in range(len(list_peak)):
            file=open('./'+str(Time_start).split('.')[0]+'/run/Local/archive/archive_'+str(i), 'r')
            assignment_archive_clusterX=cPickle.load(file)
            file.close()
            #possible_peak_assignments=getattr(result, str(list_peak[i])+'_possible_peak_assignments')
            sort_index[0,i]=list_peak[i]
            #sort_index[1,i]=len(possible_peak_assignments)
            sort_index[1,i]=len(assignment_archive_clusterX)
        score_sorting=sort_index[1,:]
        sort=np.argsort(score_sorting, axis=None)
        flag_print_result=0
      
      for cluster_index in sort:
        if P>0:
            cluster_name=int(sort_index[0,cluster_index])
            set_of_peaks=getattr(selected_peaks_clusters, str(cluster_name))
            if set_of_peaks.size<cluster_size:
              counter+=1
              continue                      
            if cluster_name in already_called_cluster:
              counter+=1
              continue

            #######IN: Skip peak assignment if overlapping already assigned peak #########
            flag_continue=0
            for overpeak in overlap_topo[cluster_name]:
              if (overpeak in already_called_cluster or overpeak in assignment_locked_peak):flag_continue=1
            if flag_continue==1:
              counter+=1
              continue
            #######OUT: Skip peak assignment if overlapping already assigned peak ########
            
            file=open('./'+str(Time_start).split('.')[0]+'/run/Local/archive/archive_'+str(cluster_name), 'r')
            assignment_archive_clusterX=cPickle.load(file)
            file.close()
 
            list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/temp/')
            file=open('./'+str(Time_start).split('.')[0]+'/run/temp/'+str(list_of_files[0]), 'r')
            assignment_previous=cPickle.load(file)
            file.close()
            
            try:involved_peak=assignment_archive_clusterX[0][0]
            except IndexError:
              counter+=1
              continue
            index_to_keep=[]
            for i in range(len(involved_peak)):
              if not involved_peak[i] in assignment_previous[0]:index_to_keep.append(i)          
            
            #print P,cluster_name,index_to_keep,involved_peak
            if len(index_to_keep)==len(involved_peak):flag_newArea=1					#EDIT
            else:flag_newArea=0														#EDIT
            
            involved_peak=[involved_peak[i] for i in index_to_keep]
            possible_peak_assignments=[]
            if len(involved_peak)==0:
              counter+=1
              already_called_cluster.append(cluster_name)
              continue            
            for line in assignment_archive_clusterX:
              #####IN: Test allowed assignment for all peaks #######
              try:
                  flag_stop=0
                  for i in range(len(line[1])):
                    peak=HMQC_peak_list[int(assignment_archive_clusterX[0][0][i])]
                    methyl=metrics[2][int(line[1][i])]
                    if peak in archive_assignment_result.keys():
                      if not methyl in archive_assignment_result[peak].keys():
                        flag_stop=1
                        break
                      else:pass
                    else:pass
                  if flag_stop==1:continue
              except NameError:pass
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
              if not line_to_add in possible_peak_assignments:possible_peak_assignments.append(line_to_add)                        
         
            if (len(possible_peak_assignments)==0):
                if not cluster_name in delayed_clusters:delayed_clusters.append(cluster_name)
                counter+=1
                continue            
            else:pass

            sum_local=peak_matrix_Scoring[involved_peak,:][:,List_of_assigned_peaks].sum()
            sum_global=peak_matrix_Scoring[List_of_assigned_peaks,:][:,involved_peak].sum()
            ratio=round(float((len(List_of_assigned_peaks)+len(involved_peak)+len(assignment_locked_peak))/float(len(list_peak))),2)
            
            mean_score=np.mean(peak_matrix_Scoring[List_of_assigned_peaks+involved_peak,:][:,List_of_assigned_peaks+involved_peak].sum(axis=0))           
            #FILTER=round(FILTER_scale*(0.03-0.02*ratio)*highest_score,3)
            #FILTER=round(FILTER_scale*len(involved_peak+List_of_assigned_peaks)*((1.05*(1-ratio))**2)*mean_score,3)           
            #FILTER=round(scaler*FILTER_scale*(FILTER_start*(1-ratio)+0.25*FILTER_start),3)
            FILTER=round(((scaler-1)*ratio+1)*FILTER_scale*(FILTER_start*(1-ratio)+0.25*FILTER_start),3)
            
            keep2test=[] 
            for peak in list_peak:
              if not int(peak) in involved_peak:keep2test.append(int(peak))
              else:pass
            sum_out=peak_matrix_Scoring[involved_peak,:][:,keep2test].sum()
            
            file_test=open('./'+str(Time_start).split('.')[0]+'/test', 'a')
            file_test.write('P: '+str(P)+(7-len(str(P)))*' '+
                            '\tclus: '+str(cluster_name)+(5-len(str(cluster_name)))*' '+
                            '\t#peak: '+str(len(involved_peak))+(3-len(str(len(involved_peak))))*' '+
                            '\tsum_out: '+str(round(sum_out,3))+(7-len(str(round(sum_out,3))))*' '+
                            '\tsum_g+l= '+str(round(sum_global+sum_local,3))+(7-len(str(round(sum_global+sum_local,3))))*' '+
                            '\ttol= '+str(round(FILTER,3))+(7-len(str(round(FILTER,3))))*' '+
                            '\tass= '+str(len(possible_peak_assignments))+(7-len(str(possible_peak_assignments)))*' '+
                            '\tflag= '+str(flag_newArea)
                            )

            if (P>5 and P<=10 and ratio<0.33):
              if len(involved_peak)==1:max_len=125
              elif len(involved_peak)==2:max_len=250
              elif len(involved_peak)==3:max_len=500
              else:max_len=750
            elif (P>=2 and P<=5 and ratio<0.5):
              if len(involved_peak)==1:max_len=250
              elif len(involved_peak)==2:max_len=500
              elif len(involved_peak)==3:max_len=1000
              else:max_len=1500
            else:         
              if len(involved_peak)==1:max_len=50
              elif len(involved_peak)==2:max_len=100
              elif len(involved_peak)==3:max_len=200
              else:max_len=300
            
            
            if (P>(P_high*0.5)
                or len(List_of_assigned_peaks)==0
                or len(possible_peak_assignments)<=5
                or (flag_newArea==1 and len(possible_peak_assignments)<=max_len and sum_out>=0.5*FILTER)
                or (flag_newArea==0 and len(possible_peak_assignments)<=max_len and 
                                                ((sum_local+sum_global)>=FILTER or
                                                 (len(involved_peak)>=2 and sum_out>=0.5*FILTER)))
                or (flag_newArea==0 and P<=3 and len(possible_peak_assignments)<=max_len 
                                             and (sum_local+sum_global)>=0.75*FILTER)
                or (P==2 and len(possible_peak_assignments)<=max_len and flag_newArea==1)
                ):
               file_test.write(' *\n')
               file_test.close()
               pass
            else:  
               if not cluster_name in delayed_clusters:delayed_clusters.append(cluster_name)
               counter+=1
               file_test.write('\n')
               file_test.close()
               continue                   

                #or (flag_newArea==1 and len(possible_peak_assignments)<=max_len and sum_out>=2*FILTER)


        file_log=open('./'+str(Time_start).split('.')[0]+'/log', 'a')
        file_log.write('Global: Tc= '+str(P)+' ; Cluster merging # '+str(cluster_name)+'\n')
        file_log.close()
        
        report_follow=open('./'+str(Time_start).split('.')[0]+'/run/Global/logFile_'+str(P), 'a')
        report_follow.write('##########  Cluster : '+str(cluster_name)+', counter: '+str(counter)+'\n')
        report_follow.write('involved_peak: '+str(involved_peak)+'\n'+'number of possible assignments: '+
                             str(len(possible_peak_assignments))+'\n'+str(possible_peak_assignments)+'\n')
        report_follow.close()

        #####################IN:Multi-processed assignment building###################
        list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/temp/')
        archive_assignment_result_previous=archive_assignment_result
        archive_assignment_result={}
        shuffle(list_of_files)
        size=len(list_of_files)/cpu
        size_reste=len(list_of_files)%cpu 
        list_of_file_index=[]
        if size!=0:
              for i in range(cpu):list_of_file_index.append([j for j in range(size*i,(i+1)*size)])
              if size_reste!=0:
                for j in range(size_reste):list_of_file_index.append([size*cpu+j])
              pool=mp.Pool(processes=cpu) 
              pool_result=pool.map(build_assignment, list_of_file_index[0:cpu])
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
        if len(list_of_file_index)>cpu:
                  pool = mp.Pool(processes=cpu) 
                  pool_result=pool.map(build_assignment, list_of_file_index[cpu:len(list_of_files)]) 
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
        list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/temp/')
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
                if int(peak_id)!=cluster_name:os.remove('./'+str(Time_start).split('.')[0]+'/run/temp/'+str(file_name))
        if flag_not_null==0:archive_assignment_result=archive_assignment_result_previous
        already_called_cluster.append(cluster_name)
        for i in involved_peak: 
              if not i in List_of_assigned_peaks:List_of_assigned_peaks.append(i)
        counter+=1
        report_follow=open('./'+str(Time_start).split('.')[0]+'/run/Global/logFile_'+str(P), 'a')
        for peak in archive_assignment_result.keys():
          for methyl in archive_assignment_result[peak]:
            if re.search('\w(\d+)\w',str(methyl)).group(1)==re.search('(\d+)\w', peak).group(1):
              report_follow.write('R  ') 
              break
          report_follow.write(str(peak)+'  '+str(archive_assignment_result[peak])+'\n')     
        report_follow.close()
               
        if flag_not_null==1:
          outputing(highest_score,peak_matrix_Scoring)  
          Time_end=time.datetime.now()
          list_of_peaks=archive_assignment_result.keys()    
          print 'Tc= '+str(P)+'; '+str(Time_end-Time_start).split('.')[0].split('.')[0],' => ',str(round(100*len(list_of_peaks)/float(len(HMQC_peak_list)),1)),'%'                 
      setattr(result,'archive_assignment_result',archive_assignment_result)

##########################################################################################
############ Assignment of isolated peak #################################################
P='iso'
#### Define not yet assigned peaks ####
isolated_peaks_list=[]
for i in range(len(HMQC_peak_list)):
  if not i in List_of_assigned_peaks:isolated_peaks_list.append(i)
  else:pass

archive_assignment_result_2save=dict(archive_assignment_result)

path_name='/run/Global'
logFile_name='logFile_iso'
list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/temp/')
           
still_not_assigned=list(isolated_peaks_list)
counter=0
no_crosspeaks_peak=[]
cycler=0
turn=0
a=0
archive_alternative_Assignments=[]

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
        sort_index[1,i]=peak_matrix_clustering[still_not_assigned[i],List_of_assigned_peaks].sum()
    score_sorting=sort_index[1,:]
    sort=np.argsort(score_sorting, axis=None)
    sort=sort[::-1]
    #### Accretion to final_assignment ####
    for index in sort:
      cluster_name=int(sort_index[0,index])
      
      file_log=open('./'+str(Time_start).split('.')[0]+'/log', 'a')
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
      report_follow=open('./'+str(Time_start).split('.')[0]+'/'+path_name+'/'+logFile_name, 'a')
      report_follow.write('##########  Peak alone : '+str(cluster_name)+'\n')
      report_follow.write('assignment table length: '+str(len(possible_peak_assignments))+'\n'+str(possible_peak_assignments)+'\n')
      report_follow.close() 

      if peak_matrix_Scoring[involved_peak[0],:].sum()==0:
          file_log.write('no crosspeaks'+'\n')
          no_crosspeaks_peak.append([cluster_name,possible_peak_assignments])
          continue         

      sum_local=peak_matrix_Scoring[involved_peak,:][:,List_of_assigned_peaks].sum()
      sum_global=peak_matrix_Scoring[List_of_assigned_peaks,:][:,involved_peak].sum()
      ratio=round((len(List_of_assigned_peaks)+len(assignment_locked_peak)+1)/float(len(list_peak)),2)

      mean_score=np.mean(peak_matrix_Scoring[List_of_assigned_peaks,:][:,List_of_assigned_peaks].sum(axis=0))           
      #FILTER=round(FILTER_scale*len(involved_peak+List_of_assigned_peaks)*((1.05*(1-ratio))**2)*mean_score,3)
      
      #FILTER=round(FILTER_scale*(0.05-0.04*ratio)*highest_score,3)
      FILTER=round(((scaler-1)*ratio+1)*FILTER_scale*(FILTER_start*(1-ratio)+0.25*FILTER_start),3)
                 
      #FILTER=round(scaler*FILTER_scale*(FILTER_start*(1-ratio)+0.25*FILTER_start),3)

      keep2test=[] 
      for peak in list_peak:
        if not int(peak) in involved_peak:keep2test.append(int(peak))
        else:pass
      sum_out=peak_matrix_Scoring[involved_peak,:][:,keep2test].sum()
      
      file_test=open('./'+str(Time_start).split('.')[0]+'/test', 'a')
      file_test.write('P: '+str(P)+(7-len(str(P)))*' '+
                      '\tclus: '+str(cluster_name)+(5-len(str(cluster_name)))*' '+
                      '\t#peak: '+str(len(involved_peak))+(3-len(str(len(involved_peak))))*' '+
                      '\tsum_out: '+str(round(sum_out,3))+(7-len(str(round(sum_out,3))))*' '+
                      '\tsum_g+l= '+str(round(sum_global+sum_local,3))+(7-len(str(round(sum_global+sum_local,3))))*' '+
                      '\ttol= '+str(round(FILTER,3))+(7-len(str(round(FILTER,3))))*' '+
                      '\tass= '+str(len(possible_peak_assignments))+'\n'
                      )
      file_test.close()

      if score_tol_end[1]=='f':FILTER=0   # ONLY HIGHEST SCORE ARE KEPT
      else:pass
      ##################### Multi-processed assignment building ####################         
      archive_assignment_result_previous=archive_assignment_result
      archive_assignment_result={}
      list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/temp/')
      size=len(list_of_files)/cpu
      size_reste=len(list_of_files)%cpu 
      list_of_file_index=[]
      if size!=0:
        for i in range(cpu):list_of_file_index.append([j for j in range(size*i,(i+1)*size)])
        if size_reste!=0:
          for j in range(size_reste):list_of_file_index.append([size*cpu+j]) 
        pool=mp.Pool(processes=cpu) 
        pool_result=pool.map(build_assignment, list_of_file_index[0:cpu])
        pool.close()
      elif size==0:
        for i in range(len(list_of_files)):list_of_file_index.append([i])  
        pool=mp.Pool(processes=cpu) 
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
      if len(list_of_file_index)>cpu:
        pool = mp.Pool(processes=cpu) 
        pool_result=pool.map(build_assignment, list_of_file_index[cpu:len(list_of_files)]) 
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
      list_of_files=os.listdir('./'+str(Time_start).split('.')[0]+'/run/temp/')
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
                if int(peak_id)!=cluster_name:os.remove('./'+str(Time_start).split('.')[0]+'/run/temp/'+str(file_name))
      List_of_assigned_peaks.append(cluster_name) 
      still_not_assigned.remove(cluster_name)
      if flag_not_null==0:archive_assignment_result=archive_assignment_result_previous                     
      list_of_peaks=archive_assignment_result.keys()
      list_of_peaks.sort()
      Time_end=time.datetime.now()      
      print 'Tc= '+str(P)+'; '+str(Time_end-Time_start).split('.')[0].split('.')[0], ' => ',str(round(100*len(list_of_peaks)/float(len(HMQC_peak_list)),1)),'%'
      
      if flag_not_null==1:
        peak=HMQC_peak_list[int(cluster_name)]
        archive_alternative_Assignments.append((peak,possible_peak_assignments))
        for peak in archive_assignment_result.keys():
          if peak in archive_assignment_result_2save.keys():
            for methyl in archive_assignment_result_2save[peak]:
              if not methyl in archive_assignment_result[peak].keys():
                archive_assignment_result[peak][methyl]=archive_assignment_result_2save[peak][methyl]
              else:pass
          else:pass
      #  for line in archive_alternative_Assignments:
      #    for methyl_index in line[1]:
      #      methyl=metrics[2][methyl_index[0]]
      #      if not methyl in archive_assignment_result[line[0]].keys():
      #        archive_assignment_result[line[0]][methyl]='nd'
      #      else:pass
      
      report_follow=open('./'+str(Time_start).split('.')[0]+'/'+path_name+'/'+logFile_name, 'a')
      for peak in archive_assignment_result.keys():
          for methyl in archive_assignment_result[peak]:
            if re.search('\w(\d+)\w',str(methyl)).group(1)==re.search('(\d+)\w', peak).group(1):
              report_follow.write('R  ') 
              break
          report_follow.write(str(peak)+'  '+str(archive_assignment_result[peak])+'\n')     
      report_follow.close()
       
      turn=0
      break
    if turn==1:cycler=1
    else:no_crosspeaks_peak=[]
    outputing(highest_score,peak_matrix_Scoring)

########STATS FILE#######################################################################
try:file=open('./'+str(Time_start).split('.')[0]+'/Output/hmqc_iso.list', 'r')
except IOError:file=open('./'+str(Time_start).split('.')[0]+'/Output/hmqc.list', 'r')
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
file_log=open('./'+str(Time_start).split('.')[0]+'/log', 'a')
for key in sorted_histo_keys:
  val=histo_ambiguity[str(key)]
  file_log.write(str(key)+': '+str(val)+'\n')
  tot_ass+=float(key)*float(val)
  tot_met+=float(val)
Time_end=time.datetime.now()
file_log.write('Mean number of alternative assignment: '+str(round(tot_ass/tot_met,1))+
           '\nCalculation time: '+str(Time_end-Time_start).split('.')[0].split('.')[0])
file_log.close
print 'Assignment complete'
