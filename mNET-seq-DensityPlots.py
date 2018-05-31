#!/usr/bin/python
#############################################################################
###                                                                       ###
###   Antisense Project - 2018                                            ###
###   Density Plots From Output of PeakCaller_mNET-Seq (Rui Luis)         ###
###   Author: Rui Luis (MCFonsecaLab)                                     ###
###                                                                       ###
#############################################################################

#Import Packages
import argparse
import pandas as pd
import matplotlib.pyplot as plt

#Input Arguments
parser = argparse.ArgumentParser(description='DensityPlot mNET-Seq')
parser.add_argument('PeakCallingFiles', nargs='+', help='Output file from PeakCalling Script for mNEt-seq (Rui Luis, 2018). For more than one replicate '
                                                                 'It will merge the results, according to the given function.')

parser.add_argument('--outNameFile', dest='outNameFile', default= None, help='Base bamFile where peaks are investigated')

args = parser.parse_args()

globalDT = pd.DataFrame()

for replicate in args.PeakCallingFiles:
    print replicate

    x = pd.read_csv(replicate, header=None, delimiter="\t", names =["chr","start","end","FeatureID",
                                                                        "Nline","Strand","DistFromTSS",
                                                                        "FeatureSize","Extrainfo","Coverage"])
    x.sort_values(["Nline"], inplace=True,ascending=True)
    x.reset_index(drop=True, inplace=True)


    globalDT = globalDT.append(x)


globalDT.sort_values(["Nline"],inplace=True)

plt.figure(1)
plt.title('Conserved peaks between replicates per Feature')
plt.xlabel("Distance From TSS")
counts = globalDT.groupby("Nline").filter(lambda x: len(x) == len(args.PeakCallingFiles)).drop_duplicates(subset=["chr","start","end","FeatureID",
                                                                        "Nline","Strand","DistFromTSS",
                                                                        "FeatureSize","Extrainfo"])
counts.sort_values(["Nline"],inplace=True)
try:
    counts["DistFromTSS"].str.strip("DistFromTSS=").astype('int64').where(counts["DistFromTSS"].str.strip("DistFromTSS=").astype('int64') <= 200)\
        .plot.hist(stacked=True, alpha=0.5, bins=200, xlim=[0,200])

    plt.show(1)
except:
    pass

plt.figure(2)
plt.title('All peaks per Feature')
plt.xlabel("Distance From TSS")

globalDT["DistFromTSS"].str.strip("DistFromTSS=").astype('int64').where(globalDT["DistFromTSS"].str.strip("DistFromTSS=").astype('int64') <= 1000)\
    .plot.hist(stacked=True, alpha=0.5, bins=200, xlim=[0,1000])
plt.show(2)

plt.figure(3)
plt.title('All peaks per Feature')
plt.xlabel("Distance From TSS")

globalDT["DistFromTSS"].str.strip("DistFromTSS=").astype('int64').where(globalDT["DistFromTSS"].str.strip("DistFromTSS=").astype('int64') <= 200)\
    .plot.hist(stacked=True, alpha=0.5, bins=75, xlim=[0,200])
plt.show(3)

########################################################
## Highest peak per Gene
plt.figure(4)
plt.title('Highest peak per Feature')
plt.xlabel("Distance From TSS")

globalDT["DistFromTSS"] = globalDT["DistFromTSS"].str.strip("DistFromTSS=")
globalDT["DistFromTSS"] = globalDT["DistFromTSS"].astype('int64')
globalDT =  globalDT.reset_index(drop=True)

CoverageDataMAX = globalDT.loc[globalDT.reset_index().groupby(['FeatureID'])['Coverage'].idxmax()]
CoverageDataMAX = CoverageDataMAX[CoverageDataMAX['DistFromTSS'] < 1200]
CoverageDataMAX['DistFromTSS'].plot.hist(stacked=True, alpha=0.5, bins=200, xlim=[0,1200])


print "0-200: {} %".format(str(sum((CoverageDataMAX['DistFromTSS'] < 200) & (CoverageDataMAX['DistFromTSS'] > 0))/65.00))
print "200-850: {} %".format(str(sum((CoverageDataMAX['DistFromTSS'] < 850) & (CoverageDataMAX['DistFromTSS'] > 200))/65.00))
print "850-1200: {} %".format(str(sum((CoverageDataMAX['DistFromTSS'] < 1200) & (CoverageDataMAX['DistFromTSS'] > 850))/65.00))
print "1200-inf: {} %".format(str(sum(CoverageDataMAX['DistFromTSS'] > 1200)/65.00))

plt.show(4)

plt.figure(5)
plt.title('Highest peak per Feature')
plt.xlabel("Distance From TSS")
CoverageDataMAX = CoverageDataMAX[CoverageDataMAX['DistFromTSS'] < 200]
CoverageDataMAX['DistFromTSS'].plot.hist(stacked=True, alpha=0.5, bins=75, xlim=[0,200])
plt.show(5)
