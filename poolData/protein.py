import os
import sys
class protein:
    
    gene_symbol = None
    uniprot_set = None 

    def __init__(self,gene_symbol, uniprot_list):
        self.gene_symbol = gene_symbol
        self.uniprot_set = set(uniprot_list)


    def checkMatch(self, queryProtein):
        #Check if the proteins are equivalent based on the current standard.
        #if intersection of two sets is non-empty, then proteins are the same.
        return (self.uniprot_set & queryProtein.uniprot_set)

    def checkMatchQ(self, queryProtein):
        """Checks if there is a match with the query protein, reeturns true if so."""
        return len(self.checkMatch(queryProtein))>0

    def findProteinMatches(self, proteinList):
        #Finds all matching proteins in a list, returns the indices
        return [ i for i in range(len(proteinList)) if self.checkMatch(proteinList[i])]
    
    def updateUniprot(self,otherProtein):
        if(self.checkMatch(otherProtein)):
            self.uniprot_set = self.uniprot_set | otherProtein.uniprot_set

    
    def readProteinProteinAnnotationsUsingPipeline(annotationName, dropHeader=True):
        
        try:
            ipmsDir = os.environ['IPMS_DIR']
        except KeyError:
            sys.stderr.write("ERROR: Environmental variable 'IPMS_DIR' not specified.\n")
            sys.stderr.write("ERROR: To fix this, add the following lines to ~/.profile\n")
            sys.stderr.write("ERROR: export IPMS_DIR=\"/path/to/PooledIPMS\"\n")
            raise Exception("Invalid environment")  
        
        return protein.readProteinProteinAnnotations("%s/referenceData/knownPPIs/%s.tsv"%(ipmsDir, annotationName), dropHeader)
        
    def readProteinProteinAnnotations(tsvFile, dropHeader=True):
        annotationList = []
        firstLine = True
        with open(tsvFile) as f:
            #Skips the header line.
            
            for l in f:
                if firstLine:
                    if dropHeader:
                        f.readline()
                    firstLine=False

                (gene1, gene2, annotation1, annotation2) = l.rstrip().split('\t')
                annotationList.append([protein(gene1, annotation1.split(";")), protein(gene2, annotation2.split(";"))])
                
        return annotationList
    
    def findMatchingProteinPairs(proteinList1, proteinList2):
        #if two different entries in the list have different gene symbols but have some overlap in uniprot
        # loop over all pairs and report all pairs which have matching uniprots
        # return list of tuples which are indices of matching proteins
        pairList = []
        for i in range(len(proteinList1)):
            for j in range(len(proteinList2)):
                if(proteinList1[i].checkMatch(proteinList2[j])):
                    match = (i,j)
                    pairList.append(match)
                else:
                    continue
        return pairList