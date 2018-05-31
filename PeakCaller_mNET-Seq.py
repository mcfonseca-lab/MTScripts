#!/usr/bin/python
#############################################################################
###                                                                       ###
###   Antisense Project - 2018                                            ###
###   PeakCalling mNET-Seq                                                 ###
###   Author: Rui Luis (MCFonsecaLab)                                     ###
###                                                                       ###
#############################################################################

#
#Input: This script expect as input a BED file and a BAM file
#Function: From the regions defined in bed file, this script will search
#the peaks following the metric:
#
#   a significant peak is one that is 3 standard deviation bigger than the mean peaks in
#   the around 200 bp (100 upstream and 100 downstream).
#Output: A bed file with all significant peaks found in the regions.

#Import Packages
import subprocess
import argparse
import os
import pysam
import random
import string
import numpy as np


#Input Arguments
parser = argparse.ArgumentParser(description='PeakCalling mNET-Seq')
parser.add_argument('bedFile',  help='Base bedFile where are defined the regions to search peaks')
parser.add_argument('bamFile', help='Base bamFile where peaks are investigated')
parser.add_argument('--outNameFile', dest='outNameFile', default= None, help='Base bamFile where peaks are investigated')

args = parser.parse_args()
if args.outNameFile == None:
    args.outNameFile = "PeakCallingOutput.bed"

def generateRandomLetterNumber(N):
    return ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in range(N))


#From bamFile Header, return chr and length of them
def getBamChrOrder(bamFile):
    return [list((x.split("\t")[1].split(":")[1], x.split("\t")[2].split(":")[1])) for x in pysam.view('-H', bamFile).split("\n") if x.startswith("@SQ")]

#Create Coord by Coord for each Feature
#Function's Author: Kenny Rebelo (MCFonseca Lab - iMM), adapted in this script to return list and not file
#Input: Bed File
#Output: Bed File divided Coordinate by Coordinate
def getCoords(features):
    splice_events = []
    #featureCoords = []

    with open(features) as featurelist:
        for feature in featurelist:
            splice_events.append(feature.strip("\n\r").split("\t"))

    output = []
    counter = 1

    for line in splice_events:
        d = 1
        dneg = int(line[2]) - int(line[1])

        if line[5] == '+':
            for i in range(int(line[1]), int(line[2])):
                output.append(line[0] + "\t" + str(i) + "\t" + str(i + 1) + "\t" + line[3] + "\t" + str(counter) + "\t" + \
                             line[5] + "\t" + "DistFromTSS=" + str(d) + "\t" + "FeatureSize=" + str(
                    int(line[2]) - int(line[1])) + "\t" + line[6] +"\n")
                d = d + 1
                counter = counter + 1
        else:
            for i in range(int(line[1]), int(line[2])):
                output.append(line[0] + "\t" + str(i) + "\t" + str(i + 1) + "\t" + line[3] + "\t" + str(counter) + "\t" + \
                             line[5] + "\t" + "DistFromTSS=" + str(dneg) + "\t" + "FeatureSize=" + str(
                    int(line[2]) - int(line[1])) + "\t" + line[6]+"\n")
                dneg = dneg - 1
                counter = counter + 1

    return [i.split() for i in output]

def sortbedToolsPerFile(bedFile,chrNames):
    args = "bedtools sort -faidx {} -i {} ".format(chrNames,bedFile,bedFile).split()
    return subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]

def coverageBed(BedFile,BamFile,chrFileSortedNames):
    args = "bedtools coverage -a {} -b {} -g {} -sorted -s -counts".format(BedFile,BamFile,chrFileSortedNames).split()
    return subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]

#Returns a list of Final Features with clips of 200 nt in each side
def addSideClips(features, N=100):
    output=[]
    for feature in features:
        output.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(feature[0],str(int(feature[1])-N),feature[1],feature[3],feature[4],feature[5],"CLIP_TO_REMOVE"))
        output.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(feature[0],feature[2],str(int(feature[2])+N),feature[3],feature[4],feature[5],"CLIP_TO_REMOVE"))
        output.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(feature[0],feature[1],feature[2],feature[3],feature[4],feature[5],"FromOriginal"))
    return output


def searchNETpeaksExon(coords):  # Input single nucleotide coverage
    # results for the exons of interest
    """Returns the position of NET-seq pauses"""

    nucleotides = []

    with open(coords, 'r') as bedcoords:
        for coord in bedcoords:
            nucleotides.append(coord.strip("\n").split("\t"))

    first = True
    coverage_count = []
    coverage_line = []
    output = []

    for enum,nucleotide in enumerate(nucleotides):

        if first:  # deals with first line
            first = False
            prevID = nucleotide[3]
            coverage_count.append(int(nucleotide[-1]))
            coverage_line.append(nucleotide)

        elif nucleotide[3] == prevID:
            prevID = nucleotide[3]
            coverage_count.append(int(nucleotide[-1]))
            coverage_line.append(nucleotide)

        if nucleotide[3] != prevID or len(nucleotides) == enum + 1: # calculates the means and standard deviation for each possible pause
            if len(coverage_line) >= 201:

                for enum,count_line in enumerate(zip(coverage_count, coverage_line)):

                    if enum >= 100 and enum < len(coverage_line)-100:

                        mean = float(np.mean(coverage_count[enum-100:enum+101]))
                        std = float(np.std(coverage_count[enum-100:enum+101]))

                        if float(count_line[0]) > (mean + (3 * std)) and float(count_line[0]) >= 4:
                            #print "\t".join(count_line[1]) + "\n"
                            output.append("\t".join(count_line[1]) + "\n")

                # deals with first line of the new exon
                prevID = nucleotide[3]
                coverage_count = []
                coverage_line = []
                coverage_count.append(int(nucleotide[-1]))
                coverage_line.append(nucleotide)

    return output


if __name__ == "__main__":
    # Generate a String of ID
    RUN_ID = generateRandomLetterNumber(8)

    #Create file of Chrmosome bedNames
    with open("chr.names_"+RUN_ID+".txt", "w") as w:
        for line in getBamChrOrder(args.bamFile):
            w.write("\t".join(line)+"\n")
        w.close()

    #Introduce 200nt Clips in each feature side, to correct calculate peak calling values
    with open(args.bedFile, "r") as r:
        with open("OriginalBedFileWITHClips_"+RUN_ID+".bed", "w") as w:
            for line in addSideClips([ x.split("\t") for x in r.readlines()]):
                w.write(line)
            w.close()
        r.close()


    #Sort features from BedFile with Clips
    with open("OriginalBedFileWITHClips_sorted_"+RUN_ID+".bed", "w") as w:
        w.write(sortbedToolsPerFile("OriginalBedFileWITHClips_"+RUN_ID+".bed", "chr.names_"+RUN_ID+".txt"))
        w.close()


    #Create a file with coordinate by coordinate position for bed file
    with open("CoorByCoordOriginalBedFileWITHClips_"+RUN_ID+".bed", "w") as w:
        for line in getCoords("OriginalBedFileWITHClips_sorted_"+RUN_ID+".bed"):
            w.write("\t".join(line)+"\n")
        w.close()

    #Sort features from CoorToCoord BedFile
    with open("CoorByCoordOriginalBedFileWITHClips_sorted_"+RUN_ID+".bed", "w") as w:
        for line in sortbedToolsPerFile("CoorByCoordOriginalBedFileWITHClips_"+RUN_ID+".bed", "chr.names_"+RUN_ID+".txt"):
            w.write(line)
        w.close()

    with open("Coverage_File_"+RUN_ID+".bed","w") as w:
        w.write(coverageBed("CoorByCoordOriginalBedFileWITHClips_sorted_"+RUN_ID+".bed", args.bamFile ,"chr.names_"+RUN_ID+".txt"))
        w.close()

    with open(args.outNameFile,"w") as w:
        for line in searchNETpeaksExon("Coverage_File_"+RUN_ID+".bed"):
            w.write(line)

    files = [str(f) for f in os.listdir('.') if os.path.isfile(f)]
    for f in files:
        if string.find(f,RUN_ID) >=1:
            os.remove(f)