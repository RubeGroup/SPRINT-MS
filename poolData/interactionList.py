from .protein import protein
#import protein

class interactionList:

    #constructor takes standardized tsv of interactors:
    proteinPairList = None
    
    def __init__(self, tsvfile):
         self.proteinPairList = self.createList(tsvfile)

    def createList(self,tsvfile):
        #take tsvfile and return list of protein object tuples which are interactors
        proteinPairList = []
        with open(tsvfile) as f:
            next(f)
            for l in f:
                d = l.rstrip().split("\t")
                tempGeneSymbol1 = d[0]
                tempGeneSymbol2 = d[1]
                if(len(d)>2):
                    tempAnnotation1 = d[2].split(";")
                    tempAnnotation2 = d[3].split(";")
                else:
                    tempAnnotation1 = ''
                    tempAnnotation2 = ''
                proteinPairList += [(protein(tempGeneSymbol1, tempAnnotation1),protein(tempGeneSymbol2, tempAnnotation2))]
        self.proteinPairList = proteinPairList
        return proteinPairList

    def findMatchesProteinVector(self, proteinOneVector, proteinTwo):
        #Given row/column protein vector, find indices with matching interactions of given protein.
        interactorIndices = []
        for i in range(len(proteinOneVector)):
            for j in range(len(self.proteinPairList)):
                if (self.proteinPairList[j][0].checkMatch(proteinOneVector[i]) and 
                    self.proteinPairList[i][1].checkMatch(proteinTwo)):
                    interactorIndices.append(i)
                elif (self.proteinPairList[j][0].checkMatch(proteinTwo) and 
                      self.proteinPairList[i][1].checkMatch(proteinOneVector[i])):
                    interactorIndices.append(i)
                else:
                    continue
        return interactorIndices
    