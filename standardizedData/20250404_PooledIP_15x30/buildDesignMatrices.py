#!/usr/bin/env python3
import numpy as np
import re
import pandas as pd 
import sys
import os
def main():

	#Identifies the relevant files
	ipmsBase=os.environ["IPMS_DIR"]
	batchName="20250404_PooledIP_15x30"
	inFile="20250407_IP_Ab_mixing_scheme.txt"
 
	annotations = {
        "ANAPC5":"Q9UJX4",
        "CCT7":"Q99832",
        "DYNC1I2":"Q13409",
        "IgG":"",
        "PSMA3":"P25788",
        "DLST":"P36957",
        "CCT2":"P78371",
        "SMG1":"Q96Q15",
        "V5":"",
        "CPSF4":"O95639",
        "MCM2":"P49736",
        "EXOSC10":"Q01780",
        "SMARCC1":"Q92922",
        "RPL23":"P62829",
        "DCTN1":"Q14203",
        "DLD":"P09622",
        "EIF3A":"Q14152",
        "RANBP9":"Q96S59",
        "SEPT2":"Q15019",
        "CSTF2":"P33240",
        "COPB1":"P53618",
        "EXOSC2":"Q13868",
        "FUS":"P35637",
        "RPL23A":"P62750",
        "RUVBL2":"Q9Y230",
        "GINS1":"Q14691",
        "PPP2CA;PPP2CB":"P67775",
        "PFDN1":"O60925",
        "DYNC1H1":"Q14204",
        "ANAPC2":"Q9UJX6"
    }
	mixingMatrix = pd.read_csv("%s/rawData/%s/%s"%(ipmsBase, batchName, inFile), sep='\t')
	abNames = mixingMatrix["Ab"]
	#Removes the antibody column, transposes
	mixingMatrix = mixingMatrix.drop("Ab", axis=1).T.fillna(0.0)
	mixingMatrix.columns = abNames
	
	with open("%s/standardizedData/%s/%s.mixingMatrix.tsv"%(ipmsBase, batchName, batchName), 'w') as f:
		f.write("\t%s\n"%("\t".join(mixingMatrix.columns)))
		f.write("\t%s\n"%("\t".join([annotations[f] for f in mixingMatrix.columns])))
		mixingMatrix.to_csv(f, sep='\t', header=False )
	
	
main()

