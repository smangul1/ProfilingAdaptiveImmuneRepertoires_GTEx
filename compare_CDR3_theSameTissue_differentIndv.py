import csv
import argparse
import sys
from glob import glob
import os
import numpy as np


import Bio
from Bio import SeqIO

import gzip

import itertools
from collections import Counter

import os

import jellyfish



def BrayCurtis(dict1,dict2):
    # Given a hash { 'species': count } for each sample, returns Compositional dissimilarity
    #BrayCurtis({'a': 0.2, 'b': 0.5, 'c': 0.3,},{'a': .4, 'b': 0.1, 'c':0.5 ,})
    #0.4
    #    # confirmed by vegan R package
    #http://www.inside-r.org/packages/cran/vegan/docs/vegdist
    # this formula is equivalent to formula with differnece
    
    
    if not dict1 or not dict2:
        return 1.0
    
    
    s=0
    s1=0.0
    s2=0.0
    
    sum1=sum(dict1.values())
    
    for key, value in dict1.items():
        dict1[key]=value/sum1
    
    sum2=sum(dict2.values())
    
    for key, value in dict2.items():
        dict2[key]=value/sum2
    
    
    

    
    shared=set(dict1.keys()) & set(dict2.keys())
    for i in shared:
        s+=min(dict1[i],dict2[i])
    return 1-2.0*s/(sum(dict1.values())+sum(dict2.values()))


def CDR3Shared(dict1,dict2):

    set1=set(dict1.keys())
    set2=set(dict2.keys())

    return set1.intersection(set2)




def fastq2readLength(inFile):
    readString=""
    infile = gzip.open(inFile, 'r')
    firstLine = infile.readline()
    
    for line in infile:
        if ">" in line:
            return len(readString.replace('\n',''))
        readString+=line
        



# mixcr2CDR3 - 10/11/2016
#Note : In case the output contains identical protein sequences and different nucleotide suquences of CDR3 correponding to different lines. It will sum up. Example of such situation - examples/CDR3_different_nt_thesame_aminoAcids.txt

def mixcr2CDR3(inFile,chain):
    clonotypes={}
    d={}
    
    d["IGH"]=0
    d["IGK"]=1
    d["IGL"]=2
    
    d["TRA"]=4
    d["TRB"]=5
    d["TRD"]=6
    d["TRG"]=7
    
    clonotypesTemp=set()
    
    with open(inFile, 'r') as f:
        
        if os.stat(inFile).st_size == 0:
            return {}
        
        csvFile = csv.reader(f,delimiter='\t')
        next(csvFile)
        
        
        
        for line in csvFile:
            
            
            
            
            
            if "*" not in line[32] and "_" not in line[32]:
                
                
                if (chain =="IGH" and line[32].startswith('C') and line[32].endswith('W')) or (chain !="IGH" and line[32].startswith('C') and line[32].endswith('F')):
                
                    if chain in line[5] and chain in line[7]:
                        if  line[32] in clonotypesTemp:
                            clonotypes[line[32]]+=float(line[2])
                        else:
                            clonotypesTemp.add(line[32])
                            clonotypes[line[32]]=float(line[2])




    sum1=sum(clonotypes.values())
    
    for key, value in clonotypes.items():
        clonotypes[key]=value/sum1
    

    return clonotypes



def distanceCDR3(d1,d2):
    d=[]
    if d1 and d2: #at least one
        for c in  list(itertools.product(d1.keys(), d2.keys())):
            
            distance=1-jellyfish.levenshtein_distance(unicode(c[0]),unicode(c[1]))/float(max( len(c[0]),len(c[1]) ))
            #print c[0],c[1],distance
            d.append(distance)
        return np.mean(d)
    else:
        return -1




#===================================

ap = argparse.ArgumentParser("python")

necessary_arguments = ap.add_argument_group("Necessary Inputs")
necessary_arguments.add_argument("dir",
                                 help="directory to save results of the analysis")
necessary_arguments.add_argument("dirFastq",
                                 help="dir with unmapped reads")
necessary_arguments.add_argument("outDir",
                                 help="outDir")
necessary_arguments.add_argument("tissue",
                                 help="tissue nr from 0 to 53")


#-----

args = ap.parse_args()

#-------


print "a. Testing distanceCDR3 function based on jellyfish.levenshtein_distance"
d1={}
d2={}


d1['aaa']=1
d1['aba']=1
d1['aca']=1

d2['a11']=1
d2['222']=1
d2['333']=1

print distanceCDR3(d1,d2)
print "distanceCDR3 works as expected!"



print "b. Testing mixcr2CDR3"

clonotype_IGH_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","IGH")
clonotype_IGK_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","IGK")
clonotype_IGL_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","IGL")
clonotype_TRA_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","TRA")
clonotype_TRB_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","TRB")
clonotype_TRD_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","TRD")
clonotype_TRG_1=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-PLZ5-0003.4.unmapped_after_rRNA_lostHuman.txt","TRG")


clonotype_IGH_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","IGH")
clonotype_IGK_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","IGK")
clonotype_IGL_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","IGL")
clonotype_TRA_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","TRA")
clonotype_TRB_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","TRB")
clonotype_TRD_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","TRD")
clonotype_TRG_2=mixcr2CDR3("/u/home/s/serghei/scratch/mixcr/clonotypes/data/clones_C1422.GTEX-NPJ7-0009.1.unmapped_after_rRNA_lostHuman.txt","TRG")


print clonotype_IGH_1
print clonotype_IGK_1
print clonotype_IGL_1
print clonotype_TRA_1
print clonotype_TRB_1
print clonotype_TRG_1
print clonotype_TRD_1

print "-----"
print clonotype_IGH_2
print clonotype_IGK_2
print clonotype_IGL_2
print clonotype_TRA_2
print clonotype_TRB_2
print clonotype_TRG_2
print clonotype_TRD_2


print "Compare File 1 : {'CVNTEAFF': 1.0} and File2: {'CASGPGGAYEQYF': 0.5, 'CASSHWGPGYYEQYF': 0.5}"
print distanceCDR3(clonotype_TRB_1,clonotype_TRB_2)
print CDR3Shared(clonotype_TRB_1,clonotype_TRB_2)

print "mixcr2CDR3 works as expected!"



#-------


if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)


updated_manifest="/u/home/s/serghei/collab/metadata/updated_manifest.csv"


#File_Name
#sample

#example
#G59822.GTEX-11P82-1326.2 GTEX-11P82-1326-SM-5HL62



samplesSet=set()
File_NameSet_RNASeq=set()
File_NameSet_rl_76=set()



dict__sample__File_Name={}
dict__File_Name__sample={}

dict__sample__body_site_s={}
dict__individual_samples={}
dict_tissuePairs_distance={}
dict_tissuePairs_pairs={}





dict_indv_tissues={}
dict_tissuePairs_exactShared={}
dict_tissuePairs2={}

listTissuePairs=[]

individuals=set()
tissues=set()
samples=set()

dict_tissue_sampleSet={}


dict_tissue_exactShared={}
dict_tissue_distance={}




chains=["IGH","IGK","IGL","TRA","TRB","TRG","TRD"]

for c in chains:
    dict_tissuePairs_exactShared[c]={}
    dict_tissuePairs_distance[c]={}
    dict_tissuePairs_pairs[c]={}

    #the same tissues, different individuals
    dict_tissue_exactShared[c]={}
    dict_tissue_distance[c]={}




print "Open updated_manifest.csv take only RNA-Seq and make dictionaries"

# Make a set of eixisting tissue pairs

with open(updated_manifest, 'r') as f:
    csvFile = csv.reader(f)
    next(csvFile)
    for line in csvFile:
        
        

        if line[2]=="RNA-Seq":
            individuals.add(line[24])
            tissues.add(line[16])
            samples.add(line[1])

            dict__sample__File_Name[line[14]]=line[1]
            samplesSet.add(line[14])
            File_NameSet_RNASeq.add(line[1])

            dict__sample__body_site_s[line[1]]=line[16]

#print line[14],line[16]


#make dict_tissue_sampleSet, to store a set of samples per tissue
for t in tissues:
    for chain in chains:
        dict_tissue_exactShared[chain][t]=[]
        dict_tissue_distance[chain][t]=[]

    dict_tissue_sampleSet[t]=set()

for s in samples:
    tissue=dict__sample__body_site_s[s]
    dict_tissue_sampleSet[tissue].add(s)




for i in individuals:
    dict_indv_tissues[i]=set()


with open(updated_manifest, 'r') as f:
    csvFile = csv.reader(f)
    next(csvFile)
    for line in csvFile:
        if line[2]=="RNA-Seq":
            indv=line[24]
            tissue=line[16]
            dict_indv_tissues[indv].add(tissue)

#print dict_indv_tissues



print "Obtain tissue pairs"

for key, value in dict_indv_tissues.iteritems():
    for i in list(itertools.combinations(value, 2)):
        listTissuePairs.append(i)




validTissuePairs=set()






for key, value in Counter(listTissuePairs).iteritems():
    if value>10:
        #print key,value
        validTissuePairs.add(key)



print "Save validTissuePairs into file"
fileOut=open(args.outDir+"/validTissuePairs.csv","w")

for t in validTissuePairs:
    fileOut.write(t[0]+","+t[1])
    fileOut.write("\n")

fileOut.close()



#valid - supported by more then 10 samples
for t in validTissuePairs:
    for chain in chains:
        dict_tissuePairs_distance[chain][t]=[]
        dict_tissuePairs_exactShared[chain][t]=[]


print "Add pairs of the same tissues. TO e"
for t in tissues:
   for chain in chains:
        dict_tissuePairs_distance[chain][(t,t)]=[]
        dict_tissuePairs_exactShared[chain][(t,t)]=[]




print "Number of tissue pairs supported by at least 10 individuals -", len([i for i in Counter(listTissuePairs).values() if i > 10])
print "Tissue pair is supported on average by this many individuals - ", np.mean([i for i in Counter(listTissuePairs).values() if i > 10])
print "Number of tissue pairs supported by 100 individuals ",len([i for i in Counter(listTissuePairs).values() if i > 100])
print "Number of tissues per individual obtained from updated_manifest.csv (RNA-Seq)",np.mean(Counter(listTissuePairs).values())





#------------
print "Scan througth every sample"
print "Use unmapped reads to get the read length"
print "Consider only samples with read length  =76"

print "TO DO : use replicates to study how consistent are the results"

numberSampleOtherReadLength=[]






print "Select samples with read length  =76 and save to the set File_NameSet_rl_76 "


file_rl_not76=open(args.outDir+"/samples_rl_not76.txt","w")

k=0
for i in individuals:
    k+=1
    string1='%s*%s*' %(args.dir,i)
    string_fastq='%s*%s*' %(args.dirFastq,i)
    
    
    file=[]
    
    
    
    
    file = glob(string1)
    file2 = glob(string_fastq)
    
    
    for f in file2:
        rl=0
        rl=fastq2readLength(f)
        if rl==76:
            fileName=f.split("afterQC_lostHuman_Fasta/")[1].split("_")[0].split(".unmapped")[0]
            if fileName in File_NameSet_RNASeq:
                #print f,fileName
                File_NameSet_rl_76.add(fileName)
        else:
            file_rl_not76.write(str(rl)+","+fileName)
            file_rl_not76.write("\n")



file_rl_not76.close()


#fix later
print "Save file with read length 76 to ..."

fileOut_76=open(args.outDir+"/samples_rl_76.csv","w")


for f in File_NameSet_rl_76:
    fileOut_76.write(f)
    fileOut_76.write("\n")
fileOut_76.close()





#=======================================================================
print "2.) The same tissue different individuals"

k=0


tissues=list(tissues)
print "Tissues ",tissues[int(args.tissue)]


for t in [tissues[int(args.tissue)]]:
    files=set()
    files.clear()
    #print t, len(dict_tissue_sampleSet[t])
    for s in dict_tissue_sampleSet[t]:
        string1='%s*%s*' %(args.dir,s)
        file = glob(string1)
        if file:
            files.add(file[0])



    for p in list(itertools.combinations(list(files), 2)): #temp
        fileName1=p[0].split("_")[1].split(".unmapped")[0]
        fileName2=p[1].split("_")[1].split(".unmapped")[0]
        
        if fileName1 in File_NameSet_rl_76      and fileName2 in  File_NameSet_rl_76:
            
            
            
            for chain in chains:
                
                clonotypes1={}
                clonotypes2={}
                clonotypes1=mixcr2CDR3(p[0],chain)
                clonotypes2=mixcr2CDR3(p[1],chain)
                
                
                
                shared_CDR3=len(CDR3Shared(clonotypes1,clonotypes2))
                dict_tissue_exactShared[chain][t].append(shared_CDR3)
                d=distanceCDR3(clonotypes1,clonotypes2)
                if d!=-1:
                    dict_tissue_distance[chain][t].append(d)





print "Save to files ------"



for chain in chains:
    print chain
    
    #number of CDR3 shared between of samples from the same tissues
    
    fileOut=open(args.outDir+"/mean_sharedCDR3_acrossInd_exactShared_"+chain+".csv","w")
    for key,value in dict_tissue_exactShared[chain].iteritems():
        if value:
            if len(value)>1:
                fileOut.write(str(key)+","+str(np.mean(value)))
                fileOut.write("\n")


    fileOut.close()
    
    
    #distance between of samples from the same tissues
    fileOut2=open(args.outDir+"/mean_sharedCDR3_acrossInd_distance_mean_"+chain+".csv","w")
    for key,value in dict_tissue_distance[chain].iteritems():
        if value:
            if len(value)>1:
                fileOut2.write(str(key)+","+str(np.mean(value)))
                fileOut2.write("\n")

    fileOut2.close()
    



    fileOut4=open(args.outDir+"/sharedCDR3_acrossInd_exactShared_"+chain+".csv","w")
    for key,value in dict_tissue_exactShared[chain].iteritems():
        if value:
            for v in value:
                fileOut4.write(str(key)+","+str(v))
                fileOut4.write("\n")


    fileOut.close()
    
    
    #distance between of samples from the same tissues
    fileOut5=open(args.outDir+"/sharedCDR3_acrossInd_distance_"+chain+".csv","w")
    for key,value in dict_tissue_distance[chain].iteritems():
        if value:
            for v in value:
                fileOut5.write(str(key)+","+str(v))
                fileOut5.write("\n")

    fileOut2.close()
    



print "DONE!!!!"
                               
                               
                               
                               
                               
                               




