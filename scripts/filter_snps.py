import argparse
import pandas as pd
import os

def filter_phage(regions,pos,alt):
    found = False
    for i in range(len(regions)):
        j=regions.ix[i]
        if pos >= j.start and pos <= j.stop:
            return 'N'
    return alt

def main():

    #argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--regions',help='Input tab separated regions file (no header)',required=True)
    parser.add_argument('-f','--file',help='Input vcf formatted SNPS file',required=True)
    parser.add_argument('-o','--output',help='Output folder',required=True)
    args=parser.parse_args()

    outputfolder = args.output
    vcfFile = args.file
    regions = args.regions
    name = vcfFile.split('/')[-1].split('.')[0]
    
    #create output folder if it doesnt exist
    if not os.path.exists(outputfolder):
    	os.mkdir(outputfolder)

    regions=pd.read_csv(regions,sep='\t',names=['id','start','stop']).sort_values('start')
    vcf = pd.read_csv(vcfFile,sep='\t',skiprows=23)
    vcf['ALT']=vcf.apply(lambda row: filter_phage(regions,row['POS'], row['ALT']), axis=1)

    with open(vcfFile) as myfile:
        head = [next(myfile) for x in xrange(23)]

    with open(outputfolder+'/'+name+'_filtered.vcf','w') as newfile:
        newfile.writelines(head)
        vcf.to_csv(newfile,mode='a',sep='\t',index=False)

main()




