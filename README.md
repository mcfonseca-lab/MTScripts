# MTScripts

###Scripts developed during Master Thesis 

####[Deeptools_HeatMap_Stranded.py](https://github.com/rluis/MTScripts/blob/Dev/deeptools_HeatMap_Stranded.py) 
\- A Script to
produce 4 HeatMaps, each of them to ProteinCoding (Forward and Reverse Reads) and Antisense (Forward and Reverse reads).
Later on, If It would be necessary for the project, you can merge the 4 plots in only one (Using p.e. Inkscape or Illustrator).
In this case please give a --zMax argument. Only in this scenario them are comparable.
        
        How to Run It:
        
        python /home/rluis/MTScripts_DEV/deeptools_HeatMap_Stranded.py ../../Y1P_mNET/rep1/hisat2/Y1P_mNET_rep1_HeLa_sorted.bam --cpu 40 --zMax 2.2



###SBATCHs

####deeptools_HeatMap_Stranded.sbatch

Sbatch to run Deeptools_HeatMap_Stranded.py script:

    How to run it:
    
    sbatch deeptools_HeatMap_Stranded.sbatch  ../../Y1P_mNET/rep1/hisat2/Y1P_mNET_rep1_HeLa_sorted.bam 
