#!/usr/bin/python
###################################################################################
###                                                                             ###
###   Antisense Project - 2018                                                  ###
###   Generate Stranded MetaProfile with two TSS (ProteinCoding and Antisense). ###
###   Author: Rui Luis (MCFonsecaLab)                                           ###
###                                                                             ###
###################################################################################



###
# Imported Packages
import os
import subprocess
from deeptools import bamCoverage
from deeptools import computeMatrix
from deeptools import plotProfile
import pysam
import argparse

parser = argparse.ArgumentParser(description='Generate Stranded MetaProfile.')
parser.add_argument('bamFile', help='Base bamFIle to be reported in the MetaProfile')
parser.add_argument('--cpu', dest='NumberCpus', default=8, help='Number of Cpus to use')
parser.add_argument('--yMax', dest='yMax', default="", help='The Maximum y value to plotProfile deeptools function')

args = parser.parse_args()

###
# General Input
bamFile = args.bamFile
basenameNoExtension = os.path.splitext(os.path.basename(bamFile))[0]
NumberCpus = args.NumberCpus

if args.yMax != "":
    yMax = "--yMax " + args.yMax
else:
    yMax = ""

######
# DeepTools

###
# Input Arguments

# bamCoverage
binSize = 10
blackListFileName = os.path.dirname(os.path.realpath(__file__)) + "/Blacklist_hg38.bed"
effectiveGenomeSize = 2913022398
normalizeUsing = "BPM"

# ComputeMatrix
bpbefore = 500
bpafter = 500

TSS1_TSS2_PLUS = '/home/rluis/Rui-testing/0nly_Annoted_APPROACH/Human/Hela/Metagenes/TSS1_TSS2_region_PLUS.bed'
TSS1_TSS2_MINUS = '/home/rluis/Rui-testing/0nly_Annoted_APPROACH/Human/Hela/Metagenes/TSS1_TSS2_region_MINUS.bed'


# Define Arguments
def produce_bamCoverage(Type, bamFile, basenameNoExtension, NumberCpus, binSize, blackListFileName, effectiveGenomeSize,
                        normalizeUsing):
    outFileBW = basenameNoExtension + "_" + Type + "_.bw"
    args_bamCoverage = "-b {} -o {} --numberOfProcessors {} --binSize  {} --blackListFileName {} --effectiveGenomeSize {} --normalizeUsing {} " \
                       "--outFileFormat bigwig".format(bamFile, outFileBW, NumberCpus, binSize, blackListFileName,
                                                       effectiveGenomeSize, normalizeUsing).split()
    bamCoverage.main(args_bamCoverage)
    return outFileBW


###
# bamCoverage
def produce_computeMatrix(Type, bpbefore, bpafter, outFileBW, binSize, region):
    outFile_computeMatrix = "matrix" + "_" + Type + "_.gz"
    args_computeMatrix = 'scale-regions -b {} -a {} -S {}  --binSize {} --outFileName {} -R {}'.format(
        bpbefore, bpafter, outFileBW, binSize, outFile_computeMatrix, region).split()
    computeMatrix.main(args_computeMatrix)
    return outFile_computeMatrix


###
# plotMetaProfile
def produce_MetaProfile(Type, basenameNoExtension, outFile_computeMatrix, yMax=yMax):
    outFile_plotProfile = "MetaProfile_" + basenameNoExtension + "_" + Type + "_.eps"
    args_plotMetaProfile = "-m {} -out {} --plotType se --plotFileFormat eps {} --yMin 0 --samplesLabel '' --plotHeight 15 --plotWidth 25 --startLabel TSS1 --endLabel TSS2 ".format( outFile_computeMatrix, outFile_plotProfile, yMax).split()
    plotProfile.main(args_plotMetaProfile)


######
# Samtools

# Samtools Inputs
StrandnessType = "fr-firstrand"  # fr-secondstrand or fr-firstrand


def filter_by_flag(bamFile, Flag):
    NameToReturn = basenameNoExtension + "_" + str(Flag) + ".bam"
    reads = pysam.view('-@', str(NumberCpus), '-hbf', Flag, bamFile)
    with open(NameToReturn, 'w') as f:
        for x in reads:
            f.write(x)


def merge_bamFiles(basenameNoExtension):
    if StrandnessType == "fr-firstrand":
        pysam.merge('-f', basenameNoExtension + "_Reverse.bam", basenameNoExtension + "_99.bam",
                    basenameNoExtension + "_147.bam")
        pysam.merge('-f', basenameNoExtension + "_Forward.bam", basenameNoExtension + "_163.bam",
                    basenameNoExtension + "_83.bam")
    elif StrandnessType == "fr-secondstrand":
        pysam.merge(basenameNoExtension + "_Forward.bam", basenameNoExtension + "_99.bam",
                    basenameNoExtension + "_147.bam")
        pysam.merge(basenameNoExtension + "_Reverse.bam", basenameNoExtension + "_163.bam",
                    basenameNoExtension + "_83.bam")


def createBAMIndex(basenameNoExtension):
    args_index_Reverse = "samtools index {}_Reverse.bam".format(basenameNoExtension).split()
    subprocess.call(args_index_Reverse)
    args_index_Forward = "samtools index {}_Forward.bam".format(basenameNoExtension).split()
    subprocess.call(args_index_Forward)


def rbindMatrix(matrix1, matrix2, outputName):
    args_rbind = "computeMatrixOperations rbind -m {} {} -o {}".format(matrix1, matrix2, outputName).split()
    subprocess.call(args_rbind)


if __name__ == '__main__':
    filter_by_flag(bamFile, '147')
    filter_by_flag(bamFile, '99')
    filter_by_flag(bamFile, '163')
    filter_by_flag(bamFile, '83')
    merge_bamFiles(basenameNoExtension)
    createBAMIndex(basenameNoExtension)

    # Create Forward and Reverse Bigwig files
    Types = ["Reverse", "Forward"]
    outFileBW = []
    for enum, directionBamFile in enumerate([basenameNoExtension + "_Reverse", basenameNoExtension + "_Forward"]):
        outFileBW.append(
            produce_bamCoverage(Types[enum], directionBamFile + ".bam", basenameNoExtension, NumberCpus, binSize,
                                blackListFileName, effectiveGenomeSize, normalizeUsing))

    # Reads Forward
    for num, AnotationType in enumerate([TSS1_TSS2_MINUS, TSS1_TSS2_PLUS]):
        outFile_computeMatrix = produce_computeMatrix(Types[num], bpbefore, bpafter, outFileBW[num],
                                                      binSize, AnotationType)

    rbindMatrix("matrix" + "_Forward_.gz", "matrix" + "_Reverse_.gz", "Final_Matrix_ProteinCoding.gz")
    produce_MetaProfile("TSS2_TSS1_FORWARD", basenameNoExtension, "Final_Matrix_ProteinCoding.gz")

    # Reads Reverse
    for num, AnotationType in enumerate([TSS1_TSS2_PLUS, TSS1_TSS2_MINUS]):
        outFile_computeMatrix = produce_computeMatrix(Types[num], bpbefore, bpafter, outFileBW[num],
                                                      binSize, AnotationType)

    rbindMatrix("matrix" + "_Forward_.gz", "matrix" + "_Reverse_.gz", "Final_Matrix_ProteinCoding.gz")
    produce_MetaProfile("TSS2_TSS1_REVERSE", basenameNoExtension, "Final_Matrix_ProteinCoding.gz")

    # Remove Files
    os.remove(basenameNoExtension + "_Reverse.bam")
    os.remove(basenameNoExtension + "_Reverse.bam.bai")
    os.remove(basenameNoExtension + "_Forward.bam")
    os.remove(basenameNoExtension + "_Forward.bam.bai")
    os.remove("matrix" + "_Forward_.gz")
    os.remove("matrix" + "_Reverse_.gz")
    os.remove("Final_Matrix_ProteinCoding.gz")
    os.remove(basenameNoExtension + "_Forward_.bw")
    os.remove(basenameNoExtension + "_Reverse_.bw")
    os.remove(basenameNoExtension + "_99.bam")
    os.remove(basenameNoExtension + "_147.bam")
    os.remove(basenameNoExtension + "_163.bam")
    os.remove(basenameNoExtension + "_83.bam")
