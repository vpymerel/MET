#!/usr/bin/env python

# Importing modules
import argparse

# import the builtin time module
import time
# Grab Currrent Time Before Running the Code
start_time = time.time()


parser = argparse.ArgumentParser(description="""           
Description
-----------
Well ... Well, well, well ...
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Author
-------
    Vincent Merel
""")

#parsing arguments
parser.add_argument("--groups",
                    type=str,
                    required=True,
                    dest="groups",
                    default=None,
                    help="A tabulated file with group info (i.e. sex, pops ...) for each ind. in the vcf \nind. group\none output will be produced per group")
                    
parser.add_argument("--vcf",
                    type=str,
                    required=True,
                    dest="vcf",
                    default=None,
                    help="The vcf to process")

                    
parser.add_argument("--contigs",
                    type=str,
                    required=True,
                    dest="contigs",
                    default=None,
                    help="contigs to keep")
                    
parser.add_argument("--output",
                    type=str,
                    required=False,
                    dest="output",
                    default=None,
                    help="the output track prefix")

#Parsing arguments
args = parser.parse_args()

Inds_Groups={}
Unique_Groups=[]

File = open(args.groups, "r")
output_geno = open(args.output+".geno", "w")
output_pos = open(args.output+".positions", "w")

for line in File:
  
  Ind, Group = line.split()[0:2]
  
  Inds_Groups[Ind] = Group
  
  if Group not in Unique_Groups:
    Unique_Groups.append(Group)

print(Unique_Groups)
File.close()

#Parsing contigs
Contigs=[]

File = open(args.contigs, "r")

for line in File:
  
  Contig = line.split()[0]

  if Contig not in Contigs:
    Contigs.append(Contig)
  
File.close()

N_Pos = 0
Written = 0

VCF = open(args.vcf, "r")

for line in VCF:
  
  if "CHROM" in line:
    
    Inds = line.split()[9:]
    Groups = [Inds_Groups[Ind] for Ind in Inds]

  elif "AF=" in line:
    
    N_Pos = N_Pos + 1
    
    if N_Pos%1000000==0:
      print(str(N_Pos)+" positions have been processed in "+str((time.time()-start_time)/60)+" min")
      
    Contig = line.split()[0]
    
    Pos = int(line.split()[1])
    
    #['0/0:0,51,255', '0/0:0,39,255' ...]
    Individuals = line.split()[9:]
    #['0/0', '0/0' ...]
    Genotypes = [Individual.split(':')[0] for Individual in Individuals]

    Ref = [sum([0 if (Allele=='1' or Allele=='.') else 1 for Allele in Genotype.split("/")]) for Genotype in Genotypes]
    NonRef = [sum([0 if (Allele=='0' or Allele=='.') else 1 for Allele in Genotype.split("/")]) for Genotype in Genotypes]

    ToWrite = []
    N_Pop_with_missing_data = 0
    
    for Unique_Group in Unique_Groups:

        Ref_Values = [Value for group, Value in zip(Groups, Ref) if group == Unique_Group] 
        Ref_Count = sum(Ref_Values)
        NonRef_Values = [Value for group, Value in zip(Groups, NonRef) if group == Unique_Group] 
        NonRef_Count = sum(NonRef_Values)
        
        N_Ind = sum([1 for group in Groups if group == Unique_Group] )
        
        ToWrite.append(str(Ref_Count))
        ToWrite.append(str(NonRef_Count))

    if Contig in Contigs:
      output_geno.write("\t".join(ToWrite)+"\n")
      output_pos.write("Myr\t"+Contig+"\t"+str(Pos)+"\t"+str(Pos)+"\n")
        
        

  '''

    NonRef = [sum(X) for X in NonRef_Alleles]

    return(NonRef)



        if (Ref_Count+NonRef_Count)<2*N_Ind:
          N_Pop_with_missing_data = N_Pop_with_missing_data + 1
          Ref_Count = 0
          NonRef_Count = 0

    if N_Pop_with_missing_data/len(Unique_Group)<0.25:
cd /home/vincent/met/GenomeScans

python3 vcfToBayPass.py \
--groups /home/vincent/MSE/RADseq/Ys.txt \
--vcf /home/vincent/MET/GenomeScans/myae.vcf \
--output /home/vincent/MET/GenomeScans/MonCul
'''     
