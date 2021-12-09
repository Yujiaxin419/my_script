#!/bin/usr/env python
#coding=utf-8

## Version:     python3.6
## Author:      yujiaxin54@126.com
## Description: This script statistics basic infomation (Gene number/fpkm/GC content) from .bed file such as A/B compartment file or TAD file
## Dependencise: Bedtools
## Usage:       python stat_ABcompatrment.py -i <in.fasta> -b <in.bed> -f <fpkm> -g <gff> -o <out.file>


import argparse
import re
import os
import sys

def getOverlapGene(bed1,ingff):
    outBed = '.'.join(ingff.split('.')[:-1]) + '.bed'
    overlapBed = '.'.join(ingff.split('.')[:-1]) + '.overlap.bed'
    print('Converting gff to bed file...')
    os.system('less ' + ingff + ' |grep \'gene\' |awk \'BEGIN{FS="\t|=|;";OFS="\t"}{print $1,$4,$5,$10}\' > ' + outBed)
    print('Finding overlap...')
    try:
        os.system('bedtools intersect -a ' + bed1 + ' -b ' + outBed + ' -loj > ' + overlapBed)
    except:
        print('bedtools no found')
        sys.exit
    overlapGeneDic = {}
    bed1Dic = {}
    ###
    print('Reading '+ bed1)
    with open(bed1, 'r') as IN:
        for line in IN:
            tmpLst = line.rstrip().split('\t')
            chr, RegionStart, RegionEnd = tmpLst[:3]
            bed1Dic[(chr, RegionStart, RegionEnd)] = line.rstrip()
    ###
    print('Reading '+ overlapBed)
    with open(overlapBed,'r') as IN:
        for line in IN:
            tmpLst = line.rstrip().split('\t')
            chr, RegionStart, RegionEnd = tmpLst[:3]
            overlapGene = tmpLst[-1]
            geneLst = overlapGeneDic.setdefault((chr, RegionStart, RegionEnd), [])
            if overlapGene != '.':
                geneLst.append(overlapGene)
                overlapGeneDic[(chr, RegionStart, RegionEnd)] = geneLst
            else: pass
    return bed1Dic,overlapGeneDic

def readfpkm(infpkm):
### Gene    fpkm
### ZS110001    0.002
    fpkmDic = {}
    with open(infpkm, 'r') as IN:
        IN.readline()
        for line in IN:
            tmpLst = line.strip().split('\t')
            fpkmDic[tmpLst[0]] = tmpLst[1]
    return fpkmDic

def readFasta(inFasta):
    print("Reading genome...")
    fastaDic = {}
    with open(inFasta,'r') as IN:
        fastaName = IN.readline().strip()[1:]
        fa = ''
        for line in IN:
            if line.startswith('>'):
                fastaDic[fastaName] = fa
                fastaName = line.strip()[1:]
                fa = ''
            else:
                fa += line.rstrip()
        fastaDic[fastaName] = fa
    return fastaDic

def statInfomation(bed1Dic, overlapGeneDic, fpkmDic, fastaDic, output):
    print("Statisticsing infomation...")
    with open(output, 'w') as OUT:
        OUT.write('#oriInfo\tgeneNumber\tfpkm\tGCcontent\n')   
        for region in overlapGeneDic:
            ### calc gene number
            geneNumber = len(overlapGeneDic[region])
            ### calc fpkm
            fpkmSum = float(0)
            for gene in overlapGeneDic[region]:
                fpkmSum += float(fpkmDic[gene])
            if geneNumber != 0:
                fpkmAva = fpkmSum/float(geneNumber)
            else :pass
            ### calc GC content
            Chr, start, end = region
            subFasta = fastaDic[Chr][int(start):int(end)]
            gcSum = len(re.findall('[gcGC]', subFasta))/(int(end)-int(start))
            ### output
            oriInfo = bed1Dic[region]
            OUT.write("{}\t{}\t{}\t{}\n".format(oriInfo, geneNumber, fpkmAva, gcSum))

def main(bed1, inGff, inFasta, inFpkm, output):
    bed1Dic, overlapGeneDic = getOverlapGene(bed1, inGff)
    fastDic = readFasta(inFasta)
    fpkmDic = readfpkm(inFpkm)
    statInfomation(bed1Dic, overlapGeneDic, fpkmDic, fastDic, output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script statistics basic infomation (Gene number/fpkm/GC content) from .bed file \
        such as A/B compartment file or TAD file")
    parser.add_argument('-i', '--fasta', required=True,
                        help='<filepath> The fasta file of genome ')
    parser.add_argument('-b', '--bed', required=True,
                        help='<filepath> The bed file ')
    parser.add_argument('-g', '--gff', required=True,
                        help='<filepath> The gff file ')
    parser.add_argument('-f', '--fpkm', required=True,
                        help='<filepath> The fpkm file ')
    parser.add_argument('-o', '--output', required=True,
                        help='<output>  The output file')
    args = parser.parse_args()
    main(args.bed, args.gff,args.fasta, args.fpkm, args.output)
