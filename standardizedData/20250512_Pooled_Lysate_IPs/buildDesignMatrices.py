#!/usr/bin/env python3
import numpy as np
import re
import pandas as pd 
import sys
import os
def main():

	#Identifies the relevant files
	ipmsBase=os.environ["IPMS_DIR"]
	batchName="20250512_Pooled_Lysate_IPs"
	inFile="20250507_IP_lysate_mixing_scheme.txt"
 
	annotations = {
        "Fus3":"P16892",
        "Ndc80":"P40460",
        "Rrp6":"Q12149",
        "Sod1":"P00445",
        "Sey1":"Q99287",
        "Stu2":"P46675",
        "Nuf2":"P33895",
        "Cox20":"Q04935",
        "Set2":"P46995",
        "Upf2":"P38798",
        "Num1":"Q00402",
        "Ded1":"P06634",
        "Tma46":"Q12000",
        "Def1":"P35732",
        "Rbg":"P39729",
        "Bim1":"P40013",
        "Sub1":"P54000",
        "Rpo41":"P13433",
        "Cdc48":"P25694",
        "Snf5":"P18480",
        "Ndt80":"P38830",
        "Dbp1":"P24784",
        "Tub4":"P53378",
        "Stm1":"P39015",
        "Nup1":"P20676",
        "Swi14":"P53965",
        "Prp4":"P20053",
        "Mrp51":"Q02950",
        "Tom70":"P07213",
        "Nmd4":"Q12129"
    }

	mixingMatrix = pd.read_csv("%s/rawData/%s/%s"%(ipmsBase, batchName, inFile), sep='\t')
	abNames = mixingMatrix["Ab"].tolist()
	#fixing mixup
	(ab0, ab2) = (abNames[0], abNames[2])
	abNames[0] = ab2
	abNames[2] = ab0
    
	#Removes the antibody column, transposes
	mixingMatrix = mixingMatrix.drop("Ab", axis=1).T.fillna(0.0)
    
	fixedABNames = []
	for i in range(len(abNames)):
		tempName = abNames[i]
		if tempName=="Upf2":
			tempName="NMD2" 
		elif tempName=="Rbg":
			tempName="Rbg1"
		elif tempName=="Swi14":
			tempName="Siw14"
		fixedABNames += [ tempName.upper() ]
			
	mixingMatrix.columns = fixedABNames 
	
	with open("%s/standardizedData/%s/%s.mixingMatrix.tsv"%(ipmsBase, batchName, batchName), 'w') as f:
		f.write("\t%s\n"%("\t".join(mixingMatrix.columns)))
		f.write("\t%s\n"%("\t".join([annotations[f] for f in abNames])))
		mixingMatrix.to_csv(f, sep='\t', header=False )
	
	
main()

