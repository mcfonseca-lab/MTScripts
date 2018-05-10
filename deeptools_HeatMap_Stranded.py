#!/usr/bin/python
#############################################################################
###                                                                       ###
###   Antisense Project - 2018                                            ###
###   Generate Stranded HeatMap.                                          ###
###   Author: Rui Luis (MCFonsecaLab)                                     ###
###                                                                       ###
#############################################################################

###
#Imported Packages
import os
import subprocess
from deeptools import bamCoverage
from deeptools import computeMatrix
from deeptools import plotHeatmap
import pysam
import argparse



parser = argparse.ArgumentParser(description='Generate Stranded HeatMap.')
parser.add_argument('bamFile', help='Base bamFIle to be reported in the HeatMap')
parser.add_argument('--cpu', dest='NumberCpus', default=8 , help='Number of Cpus to use')
parser.add_argument('--zMax', dest='zMax', default="", help='The Maximum z value to plotHeatMap deeptools function')
args = parser.parse_args()

###
#General Input
bamFile=args.bamFile
basenameNoExtension = os.path.splitext(os.path.basename(bamFile))[0]
NumberCpus=args.NumberCpus
if args.zMax != "":
    zMax="--zMax " + args.zMax
else:
    zMax=""

######
#DeepTools

###
#Input Arguments

#bamCoverage
binSize=10
blackListFileName=os.path.dirname(os.path.realpath(__file__)) +  "/Blacklist_hg38.bed"
effectiveGenomeSize=2913022398
normalizeUsing="BPM"

#ComputeMatrix
referencePoint="TSS"
bpbefore=500
bpafter=500
AntisensePLUS = '/home/rluis/Rui-testing/0nly_Annoted_APPROACH/Human/Hela/Metagenes/TSS2_Isoforms_PLUS.bed'
AntisenseMINUS = '/home/rluis/Rui-testing/0nly_Annoted_APPROACH/Human/Hela/Metagenes/TSS2_Isoforms_MINUS.bed'
ProteinCodingMINUS = '/home/rluis/Rui-testing/0nly_Annoted_APPROACH/Human/Hela/Metagenes/ProteinCodingIsoforms_MINUS.bed'
ProteinCodingPLUS = '/home/rluis/Rui-testing/0nly_Annoted_APPROACH/Human/Hela/Metagenes/ProteinCodingIsoforms_PLUS.bed'

#PlotHeatMap
heatMapColor="Greys"

###
#bamCoverage

#Define Arguments
def produce_bamCoverage(Type,bamFile,basenameNoExtension, NumberCpus, binSize, blackListFileName,effectiveGenomeSize,normalizeUsing):
    outFileBW=basenameNoExtension+"_"+Type+"_.bw"
    args_bamCoverage="-b {} -o {} --numberOfProcessors {} --binSize  {} --blackListFileName {} --effectiveGenomeSize {} --normalizeUsing {} " \
         "--outFileFormat bigwig".format(bamFile, outFileBW, NumberCpus, binSize, blackListFileName,effectiveGenomeSize,normalizeUsing).split()
    bamCoverage.main(args_bamCoverage)
    return outFileBW


###
#bamCoverage
def produce_computeMatrix(Type,referencePoint,  bpbefore, bpafter, outFileBW, binSize,region):
    outFile_computeMatrix = "matrix"+"_"+Type+"_.gz"
    args_computeMatrix = 'reference-point --referencePoint {} -b {} -a {} -S {}  --binSize {} --outFileName {} -R {}'.format(referencePoint,  bpbefore, bpafter, outFileBW, binSize,outFile_computeMatrix,region).split()
    computeMatrix.main(args_computeMatrix)
    return outFile_computeMatrix

###
#plotHeatmap
def produce_plotHeadmap(Type,basenameNoExtension,outFile_computeMatrix,heatMapColor,zMax=zMax):
    outFile_plotHeatmap = "HeatMap_"+basenameNoExtension+"_"+Type+"_.eps"
    args_plotHeatMap="-m {} -out {} --colorMap {} --sortRegions keep --plotFileFormat eps {} --zMin 0 ".format(outFile_computeMatrix, outFile_plotHeatmap,heatMapColor,zMax).split()
    plotHeatmap.main(args_plotHeatMap)


######
#Samtools

#Samtools Inputs
StrandnessType="fr-secondstrand" #fr-secondstrand or fr-firstrand

def filter_by_flag(bamFile, Flag):
    NameToReturn = basenameNoExtension+"_"+str(Flag)+".bam"
    reads = pysam.view('-@',str(NumberCpus), '-hbf', Flag ,  bamFile)
    with open (NameToReturn, 'w') as f:
        for x in reads:
            f.write(x)

def merge_bamFiles(basenameNoExtension):
    if StrandnessType == "fr-firstrand":
        pysam.merge('-f', basenameNoExtension + "_Reverse.bam", basenameNoExtension+ "_99.bam", basenameNoExtension + "_147.bam")
        pysam.merge('-f', basenameNoExtension + "_Forward.bam", basenameNoExtension+ "_163.bam", basenameNoExtension + "_83.bam")
    elif StrandnessType == "fr-secondstrand":
        pysam.merge(basenameNoExtension + "_Forward.bam", basenameNoExtension+ "_99.bam", basenameNoExtension + "_147.bam")
        pysam.merge(basenameNoExtension + "_Reverse.bam", basenameNoExtension+ "_163.bam", basenameNoExtension + "_83.bam")

def createBAMIndex(basenameNoExtension):
    pysam.index(basenameNoExtension + "_Reverse.bam")
    pysam.index(basenameNoExtension + "_Forward.bam")

def rbindMatrix(matrix1,matrix2,outputName):
    args_rbind="computeMatrixOperations rbind -m {} {} -o {}".format(matrix1, matrix2,outputName).split()
    subprocess.call(args_rbind)

if __name__ == '__main__':
    filter_by_flag(bamFile,'147')
    filter_by_flag(bamFile, '99')
    filter_by_flag(bamFile, '163')
    filter_by_flag(bamFile, '83')
    merge_bamFiles(basenameNoExtension)
    createBAMIndex(basenameNoExtension)

    #Create Forward and Reverse Bigwig files
    Types=["Reverse","Forward"]
    outFileBW = []
    for enum,directionBamFile in enumerate([basenameNoExtension + "_Reverse",basenameNoExtension + "_Forward"]):
        outFileBW.append(produce_bamCoverage(Types[enum],directionBamFile+".bam",basenameNoExtension, NumberCpus, binSize, blackListFileName,effectiveGenomeSize,normalizeUsing))
    
    #ProteinCoding Forward
    for num,AnotationType in enumerate([ProteinCodingMINUS,ProteinCodingPLUS]):
        outFile_computeMatrix = produce_computeMatrix(Types[num],referencePoint, bpbefore, bpafter, outFileBW[num], binSize,AnotationType)

    rbindMatrix("matrix"+"_Forward_.gz","matrix"+"_Reverse_.gz","Final_Matrix_ProteinCoding.gz")
    produce_plotHeadmap("ProteinCoding_FORWARD",basenameNoExtension, "Final_Matrix_ProteinCoding.gz", heatMapColor)


    #ProteinCoding Reverse
    for num,AnotationType in enumerate([ProteinCodingPLUS,ProteinCodingMINUS]):
        outFile_computeMatrix = produce_computeMatrix(Types[num],referencePoint, bpbefore, bpafter, outFileBW[num], binSize, AnotationType)

    rbindMatrix("matrix"+"_Forward_.gz","matrix"+"_Reverse_.gz","Final_Matrix_ProteinCoding.gz")
    produce_plotHeadmap("ProteinCoding_REVERSE",basenameNoExtension, "Final_Matrix_ProteinCoding.gz", heatMapColor)


    #Antisense Forward
    for num,AnotationType in enumerate([AntisenseMINUS,AntisensePLUS]):
        outFile_computeMatrix = produce_computeMatrix(Types[num],referencePoint, bpbefore, bpafter, outFileBW[num], binSize, AnotationType)

    rbindMatrix("matrix"+"_Forward_.gz","matrix"+"_Reverse_.gz","Final_Matrix_ProteinCoding.gz")
    produce_plotHeadmap("Antisense_FORWARD",basenameNoExtension, "Final_Matrix_ProteinCoding.gz", heatMapColor)

    #Antisense Reverse
    for num,AnotationType in enumerate([AntisensePLUS,AntisenseMINUS]):
        outFile_computeMatrix = produce_computeMatrix(Types[num],referencePoint, bpbefore, bpafter, outFileBW[num], binSize, AnotationType)

    rbindMatrix("matrix"+"_Forward_.gz","matrix"+"_Reverse_.gz","Final_Matrix_ProteinCoding.gz")
    produce_plotHeadmap("Antisense_REVERSE",basenameNoExtension, "Final_Matrix_ProteinCoding.gz", heatMapColor)

    #Remove Files
    os.remove(basenameNoExtension + "_Reverse.bam")
    os.remove(basenameNoExtension + "_Reverse.bam.bai")
    os.remove(basenameNoExtension + "_Forward.bam")
    os.remove(basenameNoExtension + "_Forward.bam.bai")
    os.remove("matrix"+"_Forward_.gz")
    os.remove("matrix" + "_Reverse_.gz")
    os.remove("Final_Matrix_ProteinCoding.gz")
    os.remove(basenameNoExtension + "_Forward_.bw")
    os.remove(basenameNoExtension + "_Reverse_.bw")
    os.remove(basenameNoExtension + "_99.bam")
    os.remove(basenameNoExtension+ "_147.bam")
    os.remove(basenameNoExtension+ "_163.bam")
    os.remove(basenameNoExtension+ "_83.bam")
