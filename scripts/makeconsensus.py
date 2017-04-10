#!/usr/bin/python
import argparse
import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def numpos(var):
    newpos = []
    for i in range(len(var)):
        if len(var.ix[i,'REF'])==1:
            newpos.append([var.ix[i,'POS']-1])
        else:
            a=[]
            for x in range(len(var.ix[i,'REF'])):
                a.append(var.ix[i,'POS']+x-1)
            newpos.append(a)
    return newpos

def consensus_list(ref,var):
    genome = str(ref.seq)
    seq = []
    snps=[]
    index = 0
    for i in range(len(var)):
        entry = var.ix[i]
        r = entry['REF']
        v = entry['ALT']
        positions = entry['newpos']
        if len(positions)==1:
            seq.append(genome[index:int(positions[0])])
            index = int(positions[0]+1)
            seq.append(v.lower())
        else:
            seq.append(genome[index:int(positions[0])])
            index = int(positions[-1]+1)
            seq.append(v.lower())
    seq.append(genome[index:])
    return seq


def joinList(consensus,var):
    seq = ""
    snps = []
    index=0
    varcount = 0
    snps_ref = []
    ref_pos = []
    for i in consensus:
        if i.islower():
            for x in range(len(i)):
                snps_ref.append(varcount)
                snps.append(index+1+x)
            varcount+=1
        seq=seq+i.upper()
        index = len(seq)-1
    for j in snps_ref:
        ref_pos.append(var.ix[j,'POS'])
    return seq,snps,ref_pos

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--reference',help='Input reference fasta',required=True)
    parser.add_argument('-f','--file',help='Input vcf formatted SNPS file',required=True)
    parser.add_argument('-o','--output',help='Output folder',required=True)
    args=parser.parse_args()

    reference = list(SeqIO.parse(args.reference,"fasta"))[0]
    vcfFile = args.file
    outputfolder = args.output
    sample_name = vcfFile.split('.vcf')[0]
    variants = pd.read_csv(vcfFile,skiprows=23,sep='\t')[['POS','ID','REF','ALT']]
    variants = variants.drop_duplicates()
    variants = variants.sort_values('POS')
    variants = variants.reset_index()
    variants['REFLEN']=variants['REF'].apply(len)
    variants['newpos']=numpos(variants)

    if not os.path.exists(outputfolder):
    	os.mkdir(outputfolder)

    consensus = consensus_list(reference,variants)
    cns_string,snp_locations,ref_snp_locs = joinList(consensus,variants)

    print(sample_name)
    samplename = sample_name.split("/")[-1].split('.')[0]
    print(samplename)
    cns_seq = Seq(cns_string)
    record = SeqRecord(cns_seq, id=samplename, name=sample_name, description="consensus genome")

    output_handle = open(outputfolder + "/" + samplename + ".fasta", "w")
    SeqIO.write(record, output_handle, "fasta")
    output_handle.close()

main()
