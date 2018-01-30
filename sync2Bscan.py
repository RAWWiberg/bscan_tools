#!/usr/bin/env python2.7
import sys
import argparse
import os
import csv
import operator
"""
This script takes a sync file and outputs a file for input to bayescan.

It will also output a separate file giving the chromosomal coordinates
of SNPs and their IDs in the BayeScan file.

The population ID will refer to the column in .sync file

Counts can be given as raw read counts, scaled to a fixed value
or scaled to neff
(for details of "neff" see Kolaczkowski et al., 2011 Genetics 187:245-260)

Special Feature/Bug: If each pool has an unequal number of chromosomes (e.g.
sex chromosome) then supply a space separated list of known CHROMOSOME numbers
to -n. The script will then scale counts in each sample to the sample specific
chromosome number (scale argument must be set to "neff").
"""

def neff(totc,n):
    # n is the number of individuals in the pool
    neffc= float((totc*(2*n))-1)/float(totc+(2*n))
    return int(round(neffc))

def neff_c(totc,n):
    # n is the number of individuals in the pool
    neffc= float((totc*(n))-1)/float(totc+n)
    return int(round(neffc))

def GetMajorAlleles(cnts,mincnt,minc):
    # checks if a SNP is truly biallelic across all pops
    # returns the major and minor alleles across all pops
    # also checks if there are too many INDELS
    Acnt = 0
    Tcnt = 0
    Ccnt = 0
    Gcnt = 0
    Dcnt = 0
    allele_i = {"A":0, "T":1, "C":2, "G":3,"INDEL":5}
    all_cnt = 0
    for i in cnts:
        alleles = i.split(":")
        Acnt = Acnt + int(alleles[allele_i['A']])
        Tcnt = Tcnt + int(alleles[allele_i['T']])
        Ccnt = Ccnt + int(alleles[allele_i['C']])
        Gcnt = Gcnt + int(alleles[allele_i['G']])
        Dcnt = Dcnt + int(alleles[allele_i['INDEL']])
    # only count an allele if the number of reads >= mincnt
    if Acnt >= mincnt:
        all_cnt = all_cnt + 1
    if Tcnt >= mincnt:
        all_cnt = all_cnt + 1
    if Ccnt >= mincnt:
        all_cnt = all_cnt + 1
    if Gcnt >= mincnt:
        all_cnt = all_cnt + 1
    # if there are not enough alleles, the position is not a SNP
    #print Acnt, Tcnt, Ccnt, Gcnt, all_cnt
    if all_cnt <= 1:
        return "not SNP" #not enough alleles
    counts = [("A",Acnt), ("T",Tcnt), ("C",Ccnt), ("G",Gcnt)]
    counts.sort(key=operator.itemgetter(1))
    # last two alleles are the major and minor alleles
    major_alleles = counts[-2:]
    #print major_alleles, mincnt
    # if Dcnt has >= mincnt then allele is "tainted" by INDELS
    if Dcnt >= mincnt:
        return "not SNP"# SNP tainted by INDELS
    #print Dcnt
    # if the *third* most common allele has count > mincnt then the SNP
    # is not truly biallelic. There are more than 2 alleles.
    #print counts
    if counts[1][1] > mincnt:
        #print counts[1][1]
        return "not SNP"# SNP is not biallelic
    # if the *second* most common allele has count < mincnt then the
    # allele is fixed
    #print major_alleles
    for i in major_alleles:
        if i[1] < mincnt:
            return "not SNP"# major allele considered fixed
    return major_alleles


def checkSNP(line,maxc,minc,mincnt):
    # checks the SNP line for coverage
    # 1) each population has at least minc
    # 2) no population has > maxc
    # 3) no indels
    cnts = line[3:]
    #print cnts #script tester line
    #print type(maxc), type(minc) #script tester line
    for i in range(0,len(cnts)):
        pop = [int(j) for j in cnts[i].split(":")]
        #print pop, type(sum(pop)), maxc, sum(pop) > maxc, minc, sum(pop) < minc #script tester line
        # If the coverage within a population is > maxc or < minc
        # then the SNP is considered invalid
        # coverage is counted across As,Ts,Cs and Gs,
        # Ns and INDELs *not* counted
 	if sum(pop[:-2]) > maxc or sum(pop[:-2]) < minc:
            return "not SNP"#: coverage too high or too low
        
    return cnts

def Bscandict(n,cnts,major_alleles,rescale,scale,snp_dict,snp):
    # This function adds the SNP data to a dictionary.
    # Also writes SNP chromosome and position information to a file.

    # Set allele indexes. This gives the index of each allele
    # in the sync format file (0:0:0:0)
    allele_i = {"A":0, "T":1, "C":2, "G":3}
    for p in range(0,len(cnts)):
        pop = cnts[p]
        if snp_dict.has_key(str(p+1)):
            # For each population
            maj_al=major_alleles[1]
            min_al=major_alleles[0]
            maj_alc=pop.split(":")[allele_i[maj_al[0]]]
            min_alc=pop.split(":")[allele_i[min_al[0]]]
            totc=int(maj_alc)+int(min_alc)
            #print pop,pop.split(":")[allele_i[maj_al[0]]],major_alleles
            if rescale == "1":
                # Rescale the counts
                if scale == "neff":
                    # Scale counts to the effective sample size (neff)
                    # Compute the major allele frequency
                    maj_af=float(maj_alc)/float(totc)
                    # Compute a new totc as neff
                    if type(n) == list:
                        # n is a list, pool sizes are uneven.
                        #print n[p]
                        totc=neff_c(int(totc),int(n[p]))
                    else:
                        totc=neff(int(totc),int(n))
                    # Make new major and minor allele counts from the new totc
                    # and the major allele frequency
                    maj_alc=maj_af*totc
                    maj_alc=int(round(maj_alc))
                    min_alc=totc-maj_alc
                    # Add to dictionary
                    snp_dict[str(p+1)][str(snp)]=[snp,str(totc),'2',maj_alc,min_alc]
                else:
                    # Scale counts to be out out of a defined value (scale)
                    # Compute the major allele frequency
                    maj_af=float(maj_alc)/float(totc)
                    # Make new totc == scale
                    totc=int(scale)
                    # Make new major and minor allele counts from the new totc
                    # and the major allele frequency
                    maj_alc=maj_af*totc
                    maj_alc=int(round(maj_alc))
                    min_alc=totc-maj_alc
                    # Add to dictionary
                    snp_dict[str(p+1)][str(snp)]=[snp,str(totc),'2',maj_alc,min_alc]

            else:
                # Use the raw counts
                # Add to dictionary
                snp_dict[str(p+1)][str(snp)]=[snp,str(totc),'2',maj_alc,min_alc]
        else:
            snp_dict[str(p+1)] = {}
            maj_al=major_alleles[1]
            min_al=major_alleles[0]
            # For each population
            #print pop,pop.split(":")[allele_i[maj_al[0]]],major_alleles
            maj_alc=pop.split(":")[allele_i[maj_al[0]]]
            min_alc=pop.split(":")[allele_i[min_al[0]]]
            totc=int(maj_alc)+int(min_alc)
            if rescale == "1":
                # rescale the counts
                if scale == "neff":
                    # Compute the major allele frequency
                    maj_af=float(maj_alc)/float(totc)
                    # Compute a new totc as neff
                    #print type(n), n
                    if type(n) == list:
                        totc=neff_c(int(totc),int(n[p]))
                    else:
                        totc=neff(int(totc),int(n))
                    # Make new major and minor allele counts from the new totc
                    # and the major allele frequency
                    maj_alc=maj_af*totc
                    maj_alc=int(round(maj_alc))
                    min_alc=totc-maj_alc
                    # Add to dictionary
                    snp_dict[str(p+1)][str(snp)]=[snp,str(totc),'2',maj_alc,min_alc]
                else:
                    # Compute the major allele frequency
                    maj_af=float(maj_alc)/float(totc)
                    # Make the new totc == scale
                    totc=int(scale)
                    # Make new major and minor allele counts from the new totc
                    # and the major allele frequency
                    maj_alc=maj_af*totc
                    maj_alc=int(round(maj_alc))
                    min_alc=totc-maj_alc
                    # Add to dictionary
                    snp_dict[str(p+1)][str(snp)]=[snp,str(totc),'2',maj_alc,min_alc]

            else:
                # use the raw counts
                # Add to dictionary
                snp_dict[str(p+1)][str(snp)]=[snp,str(totc),'2',maj_alc,min_alc]


def printBscandict(snp_dict,bscanfil_writer):
    # This function cycles through a created dictionary and
    # prints the data in the BayeScan format.
    popkeys=[int(pop) for pop in snp_dict.keys()]
    popkeys.sort()
    for pop in popkeys:
        bscanfil_writer.writerow(["[pop]="+str(pop)])
        snpkeys = [int(snp) for snp in snp_dict[str(pop)].keys()]
        snpkeys.sort()
        for snp in snpkeys:
            bscanfil_writer.writerow(snp_dict[str(pop)][str(snp)])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    #MAIN CONTROLS
    parser.add_argument('-filename',
                       help = 'tabix indexed .vcf file name')

    #OTHER CONTROLS
    parser.add_argument('-n', default = 50, nargs='*',
                        help = 'Pool size (default: 50).')

    parser.add_argument('-N',
                        help = 'Number of populations.')
   
    parser.add_argument('-min_cov', default = 10,
                        help = 'Minimum coverage (default: 10).')

    parser.add_argument('-max_cov', default = 500,
                        help = 'Maximum coverage (default: 500).')

    parser.add_argument('-min_count', default = 5,
                        help = 'Minimum count to call allele (default: 4).')

    parser.add_argument('-rescale', default = '0',
                        help = 'Rescale (1) or not (0) the allele frequencies.')

    parser.add_argument('-scale', default = 'neff',
                        help = 'Scale if rescaling')
  
    parser.add_argument('-prefix', default = 'bscan',
                        help = 'Prefix for snp and pop files.')

    parser.add_argument('-outdir', default = 'none',
                        help = 'Output directory for snp and pop files.')
    
    parser.add_argument('-region',default = "None",
                        help = 'Contig/Chromosome id for the regions'+\
                        ' to consider')
    
    args = vars(parser.parse_args())
    try:
        #print args #script tester line
        snp=1
        npops=args['N']
        lines = open(args['filename'], 'rb')
        # INITIALISE A DICTIONARY
        args['max_cov']=int(args['max_cov'])
        args['min_cov']=int(args['min_cov'])
        args['min_count']=int(args['min_count'])
        #print len(args['n']), type(args['n'])
        if len(args['n']) > 1 and len(args['n']) < int(npops):
            sys.stderr.write("Either supply one pool size or the number of chromosomes in each pool to -n\n")
            sys.exit()
        elif len(args['n']) == 1:
            args['n'] = int(args['n'][0])
        elif len(args['n']) == int(npops):
            pass
        snp_dict={}
        if args['outdir'] == "none":
            args['outdir'] = os.getcwd()
        if args['region'] != "None":
            # A REGION HAS BEEN SPECIFIED.
            # ONLY CONSIDER LINES THAT ARE IN THIS CHROMOSOME.
            #print type(args['region']) #script tester line
            if type(args['region']) == list:
                # REGIONS ARE IN A LIST. MULTIPLE REGIONS HAVE BEEN
                # SPECIFIED.
                sys.stderr.write("THIS FEATURE IS NOT ADDED YET, SORRY...\n")
                sys.exit
                for chrom in args['region']:
                    # INITIALISE A SNP FILE
                    snpfil=open(args['outdir']+'/'+prefix+'_'+chrom+'_snps.txt','w')
                    snpfil_writer = csv.writer(snpfil,delimiter = "\t")
                    for line in lines:
                        line=[i.replace("\n","") for i in line.split("\t")]
                        #if line[0] in chroms:
                            #print line #script tester line
                        #print cont #script tester line
                    # INITIALISE A BSCAN INPUT FILE
                    bscanfil=open(args['outdir']+'/'+args['prefix']+'_'+chrom+'_bscanin.txt','w')
                    bscanfil_writer = csv.writer(bscanfil,delimiter = "\t")
                    # PRINT ALL THE DATA FROM snp_dict
                    printBscandict(snp_dict,bscanfil_writer)

            elif type(args['region']) != str:
                # REGION SPECIFED IS A NUMBER. ASSUME IT SPECIFIES A
                # CHROMOSOME.
                sys.stderr.write("THIS FEATURE IS NOT ADDED YET, SORRY...\n")
                sys.exit
                chrom = str(args['region'])
                for line in lines:
                    line=[i.replace("\n","") for i in line.split("\t")]
                    #if line[0] == chrom:
                        #print line #script tester line
                # PRINT ALL THE DATA FROM snp_dict

            else:
                # REGION SPECIFIED IS A SINGLE STRING. ASSUME IT
                # SPECIFIES A SINGLE CHROMOSOME.
                chrom = args['region']
                # INITIALISE A SNP FILE
                snpfil=open(args['outdir']+'/'+args['prefix']+'_'+chrom+'_snps.txt','w')
                snpfil_writer = csv.writer(snpfil,delimiter = "\t")
                for line in lines:
                    line=[i.replace("\n","") for i in line.split("\t")]
                    if line[0] == chrom:
                        #print line #script tester line
                        SNP = checkSNP(line,int(args['max_cov']),int(args['min_cov']),
                               int(args['min_count']))
                        if SNP != "not SNP":
                            major_alleles = GetMajorAlleles(SNP,args['min_count'],
                                                    args['min_cov'])
                            if major_alleles != "not SNP":
                                #print major_alleles
                                #print SNP
                                Bscandict(args['N'],SNP,major_alleles,args['rescale'],
                                          args['scale'],snp_dict,snp)
                                snpfil_writer.writerow([line[0],line[1],snp])
                                snp=snp+1
                # INITIALISE A BSCAN INPUT FILE
                bscanfil=open(args['outdir']+'/'+args['prefix']+'_'+chrom+'_bscanin.txt','w')
                bscanfil_writer = csv.writer(bscanfil,delimiter = "\t")
                # PRINT ALL THE DATA FROM snp_dict
                bscanfil_writer.writerow(["[loci]="+str(snp)])
                bscanfil_writer.writerow(["[populations]="+str(npops)])
                printBscandict(snp_dict,bscanfil_writer)

                   
        else:
            sys.stderr.write("Warning: No region given, creating file for entire dataset...\n")
            chrom = "all"
            # INITIALISE A SNP FILE
            snpfil=open(args['outdir']+'/'+args['prefix']+'_'+chrom+'_snps.txt','w')
            snpfil_writer = csv.writer(snpfil,delimiter = "\t")
            for line in lines:
                line=[i.replace("\n","") for i in line.split("\t")]
                #print line #script tester line
                #print line #script tester line
                SNP = checkSNP(line,args['max_cov'],args['min_cov'],
                               args['min_count'])
                #print SNP #script tester line
                if SNP != "not SNP":
                    major_alleles = GetMajorAlleles(SNP,args['min_count'],
                                                    args['min_cov'])
                    if major_alleles != "not SNP":
                        #print major_alleles
                        #print SNP
                        Bscandict(args['n'],SNP,major_alleles,args['rescale'],
                                  args['scale'],snp_dict,snp)
                        snpfil_writer.writerow([line[0],line[1],snp])
                        snp=snp+1
                        
            # INITIALISE A BSCAN INPUT FILE
            bscanfil=open(args['outdir']+'/'+args['prefix']+'_'+chrom+'_bscanin.txt','w')
            bscanfil_writer = csv.writer(bscanfil,delimiter = "\t")
            # PRINT ALL THE DATA FROM snp_dict
            bscanfil_writer.writerow(["[loci]="+str(snp-1)])
            bscanfil_writer.writerow(["[populations]="+str(npops)])
            printBscandict(snp_dict,bscanfil_writer)

        

            
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
