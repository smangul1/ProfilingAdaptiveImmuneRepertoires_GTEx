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

#===================================

ap = argparse.ArgumentParser("python")

necessary_arguments = ap.add_argument_group("Necessary Inputs")
necessary_arguments.add_argument("dir",
                                 help="directory to save results of the analysis")
necessary_arguments.add_argument("dirFastq",
                                 help="dir with unmapped reads")
necessary_arguments.add_argument("outDir",
                                 help="outDir")





args = ap.parse_args()




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
dict__tissue_sample={}


dict_indv_tissues={}
dict_tissuePairs={}
dict_tissuePairs2={}

listTissuePairs=[]

individuals=set()
tissues=set()


print "Open updated_manifest.csv take only RNA-Seq and make dictionaries"

# Make a set of eixisting tissue pairs

with open(updated_manifest, 'r') as f:
    csvFile = csv.reader(f)
    next(csvFile)
    for line in csvFile:
        
        

        if line[2]=="RNA-Seq":
            individuals.add(line[24])
            tissues.add(line[16])

            dict__sample__File_Name[line[14]]=line[1]
            samplesSet.add(line[14])
            File_NameSet_RNASeq.add(line[1])

            dict__sample__body_site_s[line[1]]=line[16]

#print line[14],line[16]




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


#valid - supported by more then 10 samples
for t in validTissuePairs:
    dict_tissuePairs[t]=[]
    dict_tissuePairs2[t]=[]





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
            #print f,fileName
            File_NameSet_rl_76.add(fileName)
        else:
            file_rl_not76.write(str(rl)+","+f)
            file_rl_not76.write("\n")


    if k>10:
        break

file_rl_not76.close()







k=0 #temporary for TESTING
for i in individuals:
    k+=1
    string1='%s*%s*' %(args.dir,i)
    file = glob(string1)
    for p in list(itertools.combinations(file, 2)):
        #print p[0]
        #print p[1]
        fileName1=p[0].split("_")[1].split(".unmapped")[0]
        fileName2=p[1].split("_")[1].split(".unmapped")[0]

        if fileName1 in File_NameSet_RNASeq and fileName2 in  File_NameSet_RNASeq:
    
    


            clonotypes1={}
            clonotypes2={}
            clonotypes1=mixcr2CDR3(p[0],"TR")
            clonotypes2=mixcr2CDR3(p[1],"TR")
            #print p[0],clonotypes1
            #print p[1],clonotypes2
            
            
            
            
            
            
            
            beta=BrayCurtis(clonotypes1,clonotypes2)
            
            t1=dict__sample__body_site_s[fileName1]
            t2=dict__sample__body_site_s[fileName2]

            
            #print beta, dict__sample__body_site_s[fileName1], dict__sample__body_site_s[fileName2],len(CDR3Shared(clonotypes1,clonotypes2))
            if  (t1,t2) in validTissuePairs:
                shared_CDR3=len(CDR3Shared(clonotypes1,clonotypes2))
                dict_tissuePairs[(t1,t2)].append(shared_CDR3)
                dict_tissuePairs[(t1,t2)].append(beta)






fileOut=open(args.outDir+"/sharedCDR3_acrossTissues.txt","w")

for key,value in dict_tissuePairs.iteritems():
    if value:
        print key,value
        fileOut.write(str(key)+","+str(np.mean(value)))
        fileOut.write("\n")


fileOut.close()

fileOut2=open(args.outDir+"/sharedCDR3_acrossTissues_beta.txt","w")

for key,value in dict_tissuePairs2.iteritems():
    if value:
        print key,value
        fileOut2.write(str(key)+","+str(np.mean(value)))
        fileOut2.write("\n")

fileOut2.close()



