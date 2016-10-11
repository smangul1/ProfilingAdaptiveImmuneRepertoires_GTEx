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


def BrayCurtis(dict1,dict2):
    # Given a hash { 'species': count } for each sample, returns Compositional dissimilarity
    #BrayCurtis({'a': 0.2, 'b': 0.5, 'c': 0.3,},{'a': .4, 'b': 0.1, 'c':0.5 ,})
    #0.4
    #    # confirmed by vegan R package
    #http://www.inside-r.org/packages/cran/vegan/docs/vegdist
    # this formula is equivalent to formula with differnece
    
    s=0
    s1=0.0
    s2=0.0
    
    shared=set(dict1.keys()) & set(dict2.keys())
    for i in shared:
        s+=min(dict1[i],dict2[i])
    return 1-2.0*s/(sum(dict1.values())+sum(dict2.values()))



def fastq2readLength(inFile):
    readString=""
    infile = gzip.open(inFile, 'r')
    firstLine = infile.readline()
    
    for line in infile:
        if ">" in line:
            return len(readString.replace('\n',''))
        readString+=line
        




def mixcr2CDR3(inFile, clonotypes):
    
    clonotypesTemp=set()
    
    with open(inFile, 'r') as f:
        csvFile = csv.reader(f,delimiter='\t')
        next(csvFile)
        for line in csvFile:
            
            if  line[32] in clonotypesTemp:
                clonotypes[line[32]]+=float(line[2])
            else:
                clonotypesTemp.add(line[32])
                clonotypes[line[32]]=float(line[2])

#===================================

ap = argparse.ArgumentParser("python")

necessary_arguments = ap.add_argument_group("Necessary Inputs")
necessary_arguments.add_argument("dir",
                                 help="directory to save results of the analysis")
necessary_arguments.add_argument("dirFastq",
                                 help="dir with unmapped reads")
necessary_arguments.add_argument("out",
                                 help="out")






args = ap.parse_args()


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

k=0



print "Select samples with read length  =76 and save to the set File_NameSet_rl_76 "

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
            print f,fileName
            File_NameSet_rl_76.add(fileName)
        else:
            print rl, f

    
    sys.exit(1)

    for f1 in file:
        for f2 in file:
            
            fileName1=f1.split("_")[1].split(".unmapped")[0]
            fileName2=f2.split("_")[1].split(".unmapped")[0]
            
            
            
            
            if fileName1>fileName2 and os.stat(f1).st_size != 0 and os.stat(f2).st_size != 0 and fileName1 in File_NameSet and fileName2 in File_NameSet:
                fileName1=f1.split("_")[1].split(".unmapped")[0]
                fileName2=f2.split("_")[1].split(".unmapped")[0]

                clonotypes1={}
                clonotypes2={}
                mixcr2CDR3(f1, clonotypes1)
                mixcr2CDR3(f2, clonotypes2)
                
                #print fileName1,fileName2,len(clonotypes1), len(clonotypes2), len(clonotypes1 & clonotypes2)


                tissue1=dict__sample__body_site_s[fileName1]
                tissue2=dict__sample__body_site_s[fileName2]

                if tissue1>=tissue2:
                    tp=tissue1+"---"+tissue2
                else:
                    tp=tissue2+"---"+tissue1
                
                if tissue1==tissue2:
                    print "!!!",f1,f2,tissue1
                beta=BrayCurtis(clonotypes1,clonotypes2)
                tissues_pairs_dict2[tp].append(beta)

    if k>10:
        break




        
target = open(args.out, 'w')

for key, value in tissues_pairs_dict2.iteritems():
    
    if len(value)>0:

        target.write( key+","+str(np.mean(value)) + ","+str(np.std(value)))
        target.write("\n")







