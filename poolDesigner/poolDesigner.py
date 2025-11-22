from abc import ABC, abstractmethod
import numpy as np
import itertools
from itertools import takewhile, product
import matplotlib.pyplot as plt
import time




class matrixLoss(ABC):

    def __init__(self, matrix):
        self.matrix = matrix

    @abstractmethod
    def loss(self):
        pass



#############################################################################################
####            Loss functions                                                   ########
#############################################################################################




class basicLoss(matrixLoss):

    def __init__(matrix):
        super().__init__(matrix)
    
    def loss(self):
        mixingMatrix = self.matrix
        cols = mixingMatrix.shape[1]
        nPools = mixingMatrix.shape[0]
        distMatrix = np.zeros((cols,cols))
        total = 0.0
        # all unique unordered column pairs
        for i, j in itertools.combinations(range(cols), 2):
            diff =  mixingMatrix[:, i] - mixingMatrix[:, j]
            dist_sq = np.linalg.norm(diff, 2)**2  # squared L1 norm
            distMatrix[i][j] = 1.0 / (0.001 + dist_sq)
         
            total += 1.0 / (0.001 + dist_sq)
        return round(total,3)
    

class fastLoss(matrixLoss):
    """
        Compute loss = sum(1/distance) for selected columns,
        where distance = squared L2 norm + 0.001,
        excluding self-distances (no diagonal terms).
    """
    def __init__(matrix):
        super().__init__(matrix)
 


    def loss(self,matrix):
        M = matrix
        n_row, n_col = M.shape
        
        # Indices of selected columns (every n_row-th column)
        selected_indices = np.arange(0, n_col, n_row)
        
        # Precompute norms and dot products
        norms = np.sum(M**2, axis=0, keepdims=True)   
        dot_products = M.T @ M                        
        
        # Squared distances
        dists_squared = norms.T + norms - 2 * dot_products
        
        # Add 0.001
        dists_adjusted = dists_squared + 0.001
        
        # Select relevant rows (only those columns we care about)
        selected_dists = dists_adjusted[selected_indices, :]
        
        # Vectorized exclusion of self-distances:
        # set diagonal entries (self-pairs) to inf in one shot
        selected_dists[np.arange(len(selected_indices)), selected_indices] = np.inf
        
        # Reciprocal
        reciprocal = 1.0 / selected_dists
        
        # Final scalar loss
        return round(np.sum(reciprocal),3)
    






    

##############################################################################
######## Matrix structure: currently only for cyclic design ########
##############################################################################


class poolDesigner:
    """
    Class containing all data specifying a pool designer algorithm
    concentrationList = list of concentrations that specific antibody is added to per pool
    concentrationMultiplicities = how many times is each concentration repeated per antibody
    lossFunction = which loss function is used to evaluate matrices
    """


    nPools = None
    concentrationList = None
    concentrationMultiplicities = None
    lossFunction = None

    def __init__(self, nPools, concentrationList, concentrationMultiplicites, lossFunction):
        self.nPools = nPools
        self.concentrationList = concentrationList
        self.concentrationMultiplicities = concentrationMultiplicites
        self.lossFunction = lossFunction


    #############################################################################################
    # Functions to generate and filter vector lists from which we choose a reduced optimal set #
    #############################################################################################

    def generateCombinations(self, nArray,p):
        """Take the elements of nArray which is a list and obtains all combinations of p elements, without replacement"""
        """essentially returns an nCp-length list, where n=len(nArray)"""
        return list(itertools.combinations(nArray, p))

    def generatePermutations(self, nArray,p):
        """Take elements of the given list called nArray and generate all permutations of length p"""
        """eg: from nArray=[1,2,3] we get [(1,2),(2,1),(1,3),(3,1),(2,3),(3,2)]"""
        return list(itertools.permutations(nArray, p))

    def generateVectors(self, baseArray, element, p):
        """Given a nonempty array, place element in p positions in that array and return all possible combinations of the placements """
        """First check if nonEmptyArray has p available positions or not"""
        nonEmptyArrayList = []
        if ((np.size(baseArray) - np.count_nonzero(baseArray))<p):
            print("Not enough empty positions!")
            return nonEmptyArrayList
        filledIndices = np.nonzero(baseArray)
        allIndices=np.arange(len(baseArray))
        remainingIndices = [i for i in allIndices if i not in filledIndices[0]]
        allCombinationslist = self.generateCombinations(remainingIndices,p) #list(itertools.combinations(remainingIndices,p))
        for i in range(len(allCombinationslist)):
            newArray = np.copy(baseArray)
            newArray[[allCombinationslist[i]]] = element
            nonEmptyArrayList.append(newArray)
        return nonEmptyArrayList


    def filterVectorsForFirstStep(self,arrayList):
        """Taking advantage of rotational symmetry we discard a lot of vectors which are in the same equivalence class"""
        """"""
        arraysMod = []
        for i in range(len(arrayList)):
            gaps = {}
            nonZeroIndices = np.nonzero(arrayList[i])[0]
            for j in range(len(nonZeroIndices)):
                if(j==(len(nonZeroIndices)-1)):
                    gaps[nonZeroIndices[j]] = len(arrayList[i]) - nonZeroIndices[j]
                else:
                    gaps[nonZeroIndices[j]] = nonZeroIndices[j+1] - nonZeroIndices[j]
            maxGap = max(gaps.values())
            maxKeys = [k for k,v in gaps.items() if v==maxGap]
            if(nonZeroIndices[-1] in maxKeys and arrayList[i][0]!=0):
                arraysMod.append(arrayList[i])
        return arraysMod

       



    ######################################################################
    # Functions which do the value placements and evaluation in the loop #
    #               after the first (highest) value is placed            # 
    ###################################################################### 

    def evaluatePlacement(self,arrayList, vecCombs):
        #out of all combinations, just choose the ones which minimizes the loss function

        #vecCombs = getPlacements(np.arange(len(arrayList)), nBlocks) #first create a list of n-tuples, where n=number of blocks and each tuple contains combinations of arrays
        rankingDict = {}   #this will contain the vector indices as keys, and the evaluated loss function as value -- sort this to obtain best vector combination
        lossStart = time.process_time()
        for i in range(len(vecCombs)):
            vecsToTest = [arrayList[k] for k in list(vecCombs[i])]
            rankingDict[vecCombs[i]] = self.lossFunction.loss(self.vectorsToMatrix(vecsToTest)) #lossFunction(vectorsToMatrix(vecsToTest))
        rankingDictSorted = dict(sorted(rankingDict.items(), key = lambda item: item[1]))
        minVal = next(iter(rankingDictSorted.values()))
        leastLoss = dict(takewhile(lambda kv: kv[1] == minVal, rankingDictSorted.items()))
       
        #return only those vectors which are tied for the lowest loss function value
        return leastLoss, list(leastLoss.keys())



    def placeIntermediateVals(self,bestKeys, arrayList, element, p):
        """function for placing intermediate values (i.e. all values after highest)"""
        """bestKeys: vector pairs from arrayList with smallest loss"""
        """element: concentration that occurs p times and needs to be tested"""
        """returns all vectors generated by placing values in the restricted vector"""
        """retains "inheritance", i.e. if a set of vectors {V}_i is contained in bestKeys, only their children are used as sets/groups"""
        allArrays = []
        globalListIndex = 0
        allowedCombs = []
        for i in range(len(bestKeys)):
            vecSet = list(bestKeys[i])
            #we need the list of "allowed" combinations to correctly choose vectors for testing later
            listOfIndices = []
            for j in range(len(vecSet)):
                localListIndex = 0
                listOfRestricted = self.generateVectors(arrayList[vecSet[j]], element, p) 
                localListIndex += len(listOfRestricted)
                listOfIndices.append(list(np.arange(globalListIndex, globalListIndex+localListIndex)))
                allArrays += listOfRestricted
                globalListIndex+=localListIndex
            allowedCombs += list(product(*listOfIndices)) #[c for tup in list(product(*listOfIndices)) for c in set(itertools.permutations(tup))]
        """thus allArrays contains all arrays obtained from placing values in restricted arrays, but allowed combs only picks the children of grouped vectors"""
        return allArrays, allowedCombs






    ####################################################################
    # Methods for assembling the mixing matrix from generating vectors #
    ####################################################################

    def blockGenerator(self,array):
        """Given an array, generate all right permutations of the array and stack, to generate block matrix"""
        A = np.copy(array)
        B = A
        for i in range(1,len(array)):
            B = np.vstack((B, np.roll(A,i)))
        blockMatrix = B
        return blockMatrix

    def vectorsToMatrix(self,vectorList):
        """given list of generating vectors, obtain the full mixing matrix"""
        """NOTE: all vectors in vectors to matrix should be of the same size"""

        firstBlock = self.blockGenerator(vectorList[0])
        for i in range(len(vectorList)-1):
            nextBlock =self.blockGenerator(vectorList[i+1])
            firstBlock = np.vstack((firstBlock,nextBlock))

        return np.transpose(firstBlock)




    ###############################################################################################
    # Final matrix search functions #
    ###############################################################################################
        

    def generateCyclicMixingMatrix(self, nBlocks):
        '''
        function that puts everything together and returns the optimal set of mixing matrices
        nPools is the size of the mixing vector to be tested, concList contains list of all concentrations from highest to lowest
        concMultiplicities indicates how many times each concentration occurs from highest to lowest in the Ab's mixing vector
        nBlocks tells how many unique generating vectors we need
        '''

        #for placing the first "HIGHEST" value from concentration list, initialize an empty array then place "HIGHEST" values in it
        nullVec = np.zeros(self.nPools)
        arrayList = self.generateVectors(baseArray = nullVec, element=self.concentrationList[0], p=self.concentrationMultiplicities[0]) 
        
        #add the restriction that we choose only those vectors where the first entry is "high" and the distance between the last entry and first entry is maximized
        #this narrows down the search space and removes (cyclically) degenerate vectors

        arrayListMod = self.filterVectorsForFirstStep(arrayList)


        #now evaluate all nBlock-sets of these vectors using the loss function and return those with smallest loss
        vecCombs = self.generateCombinations(np.arange(len(arrayListMod)), nBlocks)
      
        bestPlacements, bestKeys = self.evaluatePlacement(arrayListMod, vecCombs)

        #now do the placement/evaluation cycle iteratively len(multiplicities)-1 times
        # bestFirstplacements, bestKeys and arrayList keep getting updated as we add all the Abs
        for i in range(len(self.concentrationMultiplicities)-1):
            
            #in every loop, store the list of vector generated
            #also store optimal set indices in bestKeys
            allArrays, allowedCombs = self.placeIntermediateVals(bestKeys, arrayListMod, self.concentrationList[i+1],self.concentrationMultiplicities[i+1])
         
            bestPlacements, bestKeys = self.evaluatePlacement(allArrays, allowedCombs)
            arrayListMod = allArrays

        result = bestKeys  

        #using only the 0'th result here is an arbitrary choice, the results are supposed to be equally good.
        finalVectorList = []
        for i in range(len(result[0])):
            finalVectorList.append(allArrays[result[0][i]]) #[i for i in arrayList if i in list(result[0])]
        bestMixingMatrix = self.vectorsToMatrix(finalVectorList) 
        return bestMixingMatrix 