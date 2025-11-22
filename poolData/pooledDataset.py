import numpy as np
from .protein import protein
import os
import sys
import matplotlib.pyplot as plt
import sklearn.decomposition
            
class pooledDataset:
    """Class containing all data specifcying a pooled IP-MS experiment."""

    #Numpy array containing mixing matrix. Indices: 
    mixingMatrix = None
    #Numpy array containing the signal matrix. Indices:
    signalMatrix = None
    
    #Lists of protein objects
    preyProteins     = None
    abProteins       = None
    
    #Information about intercepts
    intercept        = None
    interceptOption  = None
    customIntercep   = None
    interceptNames   = None
    
    #Dimensions of datasets
    nPreys       = None
    nBaits       = None
    nPools       = None
    nIntercepts  = None

    #constructor takes in standardised tsv files (with uniprot IDs)
    # The files can be specified in three ways:
    # 1. Specify both mixingMatrixPath and signalMatrixPath
    # 2. Specify batchName, meaning
    #       signalMatrixPath=standardizedData/<batchName>/<batchName>.<signalChoice>.tsv 
    #       signalMatrixPath=standardizedData/<batchName>/<batchName>.mixingMatrix.tsv 
    # 3. Specify both batchName and experimentName, meaning 
    #       signalMatrixPath=standardizedData/<batchName>/<experimentName>.<signalChoice>.tsv 
    #       signalMatrixPath=standardizedData/<batchName>/<experimentName>.mixingMatrix.tsv 
    

    def __init__(self, mixingMatrix, signalMatrix, abProteins, preyProteins, intercept, interceptNames, interceptOption, customIntercept=None):
            
        self.mixingMatrix    = mixingMatrix
        self.signalMatrix    = signalMatrix
        self.abProteins      = abProteins
        self.preyProteins    = preyProteins
        
        (self.nPools, self.nBaits) = mixingMatrix.shape
        self.nPreys = len(signalMatrix)

        #Stores the values so that the intercept can be re-computed
        self.interceptOption = interceptOption
        self.customIntercep  = customIntercept


        #Computes the intercept matrix if (if not allready computed)
        if intercept is None:
            (self.intercept, self.interceptNames) = self.buildInterceptMatrix()
        else:
            (self.intercept, self.interceptNames) = (intercept, interceptNames)
            
        
        #Number  of intercepts, and predictors
        self.nIntercepts = self.intercept.shape[1]
        self.interceptColumns = tuple(range(self.nIntercepts))
        
    
    
    def buildInterceptMatrix(self, newSignalMatrix=None):
        if newSignalMatrix is None:
            newSignalMatrix = self.signalMatrix
            
        #Builds thed intercept matrix
        if (self.interceptOption == 'ones'):
            #Classic intercept
            intercept = np.ones(self.nPools)[:,np.newaxis]
            interceptNames = [ protein("BG: Intercept",[])]
            
        elif(self.interceptOption == 'means'):
            #Mean
            colMeans = np.mean(newSignalMatrix, axis=0)
            colMeans /= np.mean(colMeans)
            intercept = colMeans[:,np.newaxis]
            interceptNames = [ protein("BG: Mean",[])]
        elif (self.interceptOption == 'means_ones'):
            #Classic intercept, and mean (subtracting smallest value)
            colMeans = np.mean(newSignalMatrix, axis=0)
            #Subtracting the minimum value
            colMeans -= np.min(colMeans)
            colMeans /= np.mean(colMeans)
            intercept = \
                np.concatenate((colMeans[:,np.newaxis],
                                np.ones(self.nPools)[:,np.newaxis]),axis=1)
            interceptNames = [ protein("BG: Intercept",[]), protein("BG: Mean",[])]
            
        elif (self.interceptOption == "nnmf2"):
            #Non-negative matrix factorization using two components;
            #Runs the factorization
            (w,h,nIt)=sklearn.decomposition.non_negative_factorization(newSignalMatrix, n_components=2)
            intercept = np.transpose(h)/np.mean(h,axis=1)
            interceptNames = [ protein("BG: NNC1",[]), protein("BG: NNC2",[])]
            
        elif(self.interceptOption == 'custom'):
            #Custom intercfept
            intercept = self.customIntercept
            interceptNames = [ protein("BG: Custom %d"%i,[]) for i in range(self.customIntercept.shape[1])]
        else:
            raise Exception("Invalid intercept option: %s"%self.interceptOption)  
        return (intercept, interceptNames)
        
    def resetIntercept(self, interceptOption, customIntercept=None, newSignalMatrix=None):
        self.interceptOption=interceptOption
        self.customIntercep = customIntercept
        (intercept, interceptNames) = self.buildInterceptMatrix(newSignalMatrix=newSignalMatrix)
        return pooledDataset(self.mixingMatrix, self.signalMatrix, self.abProteins, self.preyProteins, intercept, interceptNames, interceptOption, customIntercept = customIntercept)
    
    #Reads dataset from files
    def readDatsetFromFile(mixingMatrixPath, signalMatrixPath, interceptOption, customIntercept=None):

        (abProteins, mixingMatrix)   = pooledDataset.readMixingMatrix(mixingMatrixPath)
        (preyProteins, signalMatrix) = pooledDataset.readSignalMatrix(signalMatrixPath)
        
                    
        return pooledDataset(mixingMatrix, signalMatrix, abProteins, preyProteins, None, None, interceptOption, customIntercept = customIntercept)


    #Reads standardized dataset
    def readStandardizeDataset(batchName, interceptOption, experimentName=None, customIntercept=None, signalChoice='intensity'):
      
        try:
            ipmsDir = os.environ['IPMS_DIR']
        except KeyError:
            sys.stderr.write("ERROR: Environmental variable 'IPMS_DIR' not specified.\n")
            sys.stderr.write("ERROR: To fix this, add the following lines to ~/.profile\n")
            sys.stderr.write("ERROR: export IPMS_DIR=\"/path/to/PooledIPMS\"\n")
            raise Exception("Invalid environment")  
        
        
        if experimentName != None:
            batchExperimentName = "%s/%s"%(batchName, experimentName)
        else:
            batchExperimentName = "%s/%s"%(batchName, batchName)
            
            mixingMatrixPath  = "%s/standardizedData/%s.mixingMatrix.tsv"%(ipmsDir, batchExperimentName)
            signalMatrixPath = "%s/standardizedData/%s.%s.tsv"%(ipmsDir, batchExperimentName, signalChoice)
        
        return pooledDataset.readDatsetFromFile(mixingMatrixPath, signalMatrixPath, interceptOption, customIntercept=customIntercept)

    #Reads the signal matrix file and prey protein information
    def readSignalMatrix(signalMatrixPath):
        preyProteins = []
        tempSignal   = []
        with open(signalMatrixPath) as f:
            next(f)
            for l in f:
                d = l.rstrip().split("\t")
                tempGeneSymbol = d[1]
                tempAnnotations = d[2].split(";")
                preyProteins += [ protein(tempGeneSymbol, tempAnnotations)] 
                tempSignal += [ [float(di) for di in d[3:]]]
        
        return (preyProteins, np.array(tempSignal))

    #Reads the mixing matrix file and ab protein information
    def readMixingMatrix(mixingMatrixPath):
        abProteins = []
        tempMixing = []
        with open(mixingMatrixPath) as f:
            lines = f.read().splitlines()
            tempGeneSymbol = lines[0].split('\t')[1:]
            tempAnnotations = lines[1].split('\t')[1:]
            tempMixing = [list(map(float,line.split('\t')[1:])) for line in lines[2:] if line.strip()]
            abProteins += [protein(gene,[uniprot]) for gene, uniprot in zip(tempGeneSymbol,tempAnnotations)]

        return (abProteins, np.array(tempMixing))
    
    
    def getDesignMatrix(self):
        return np.concatenate((self.intercept,self.mixingMatrix/np.max(self.mixingMatrix)), axis=1)

    
    
    def computeRankedValue(self, n, type="mean"): # [“mean” or “value”]):
        #for each protein, this function computes either “the mean of the top n values” 
        # or the “top n’th value”. Returns a list (over proteins)
        preyList = []
        for i in range(len(self.preyProteins)):
            signalArray = self.signalMatrix[i,:]
            if(type == 'mean'):
                preyList.append(np.mean(np.sort(signalArray)[-n:]))
            else:
                preyList.append(np.sort(signalArray)[-n])
                
        return preyList
    
    def normalizeByRankedValue(self,n, type="mean"): #[“mean” or “value”])
        #calls computeRankedValue and uses the value to normalize the rows, returns a new pooledDataset
        
        inputScale = self.computeRankedValue(n, type)
        for i in range(len(inputScale)):
            if(inputScale[i] ==0):
                inputScale[i] = 1

        return self.rescalePreys(1.0/np.array(inputScale))
    
    def writeSignalMatrix(self, path):
        with open(path, "w+") as f:
            #Gene symbol
            f.write("\t\t%s\n"%"\t".join([p.gene_symbol for p in self.abProteins]))
            #UniProt ID
            f.write("\t\t%s\n"%"\t".join([";".join(p.uniprot_set) for p in self.abProteins]))
            #Writes the  data
            for i in range(self.nPreys):
                prey = self.preyProteins[i]
                f.write("%s\t%s\t%s\n"%(prey.gene_symbol, ";".join(prey.uniprot_set), "\t".join(["%f"%di for di in self.signalMatrix[i]])))
                
        
    ###############################################
    # FUNCTIONS FOR NORMALIZING THE SIGNAL MATRIX # 
    ###############################################
    
    #Rescales the columns in the signal matrix, returns a new matrix
    def rescalePools(self, rescaleVector):
        temp = np.array(rescaleVector).reshape(1,-1)
        matrix_rescaled = self.signalMatrix * temp
        
        (self.intercept, self.interceptNames) = self.buildInterceptMatrix(newSignalMatrix=matrix_rescaled)

        return pooledDataset(self.mixingMatrix, matrix_rescaled, self.abProteins, self.preyProteins, self.intercept, self.interceptNames, self.interceptOption, customIntercept=self.customIntercep)

    #Rescales the rows in the signal matrix, returns a new.
    def rescalePreys(self, rescaleVector):
        
        normalized_matrix = np.copy(self.signalMatrix) * rescaleVector[:,None]
        
        (self.intercept, self.interceptNames) = self.buildInterceptMatrix(newSignalMatrix=normalized_matrix)
   
        return pooledDataset(self.mixingMatrix, normalized_matrix, self.abProteins, self.preyProteins, self.intercept, self.interceptNames, self.interceptOption, customIntercept=self.customIntercep)
    
    def medianNormalizeMatrix(self,normalizeBy):
        matrix = self.signalMatrix
        matrix_normalized = np.zeros((matrix.shape[0], matrix.shape[1]))
        
        if(normalizeBy == 'row'):
            rowMedians = np.median(matrix, axis=1)
            for i in range.len(rowMedians):
                if rowMedians[i] == 0.:
                    rowMedians[i] = 1
                    
                return self.rescalePreys(1/rowMedians)
            
        elif(normalizeBy == 'column'):
            columnMedians = np.median(matrix, axis=0)
            for i in range.len(columnMedians):
                if columnMedians[i] == 0.:
                    columnMedians[i] = 1

                return self.rescalePools(1/columnMedians)
    
    def meanNormalizeMatrix(self, normalizeBy):
        matrix = self.signalMatrix
        matrix_normalized = np.zeros((matrix.shape[0], matrix.shape[1]))
        
        if(normalizeBy == 'row'):
            rowMeans = np.mean(matrix, axis=1)
            for i in range.len(rowMeans):
                if rowMeans[i] == 0.:
                    rowMeans[i] = 1
                    
                return self.rescalePreys(1/rowMeans)
            
        elif(normalizeBy == 'column'):
            columnMeans = np.mean(matrix, axis=0)
            for i in range.len(columnMeans):
                if columnMeans[i] == 0.:
                    columnMeans[i] = 1

                return self.rescalePools(1/columnMeans)
  
    def computePoolNomalizationVector_L2(self, r2cutoff=0.8, maxPreys=20, verbose=False):
        
        #STEP 1 - For each AB/bait, identify the preys with correlation with r^2 above cuttoff (ordered by correlation), giving a dictionary
       
        highCorrelators = {}
        for i in range(len(self.abProteins)):
            for j in range(len(self.preyProteins)):
                if(np.corrcoef(self.mixingMatrix[:,i], self.signalMatrix[j,:])[0,1] > r2cutoff):
                    if(i not in highCorrelators):
                        highCorrelators[i] = []
                    highCorrelators[i].append(j)

        #STEP 2 - Truncate the lists to at most have maxPreys for each AB/bait.
        
        for key,value in highCorrelators.items():
            highCorrelators[key] = highCorrelators[key][:maxPreys]
            
        if verbose:
            print("    Number of correlated proteins: %s"%(",".join(["%d"%len(highCorrelators[key]) for key in highCorrelators])))
        
        #STEP 3 - Create a dictionary indicating what columns of the design matrix should be used for each included prey
        
        inputPredictors = {}
        for baitIndex, preyIndices in highCorrelators.items():
            for item in preyIndices:
                if item not in inputPredictors:
                    inputPredictors[item] = list(self.interceptColumns)  # Starts with the intercept columns.
                inputPredictors[item].append(baitIndex+self.nIntercepts)

        #STEP 4 - Compute the 'Q' matrix: Q_{ab} = Sum_i y_{ia} *  (I - X_i(X_i^TX_i)^-1X_i^T)_{ab} * y_{ib}
        dm        = self.getDesignMatrix()
        nPools = dm.shape[0]
        final_mat = np.zeros((nPools, nPools))

        #for key,value in inputPredictors.items():
        for key,cols in inputPredictors.items():
            Xtemp = dm[:,cols]
            y = self.signalMatrix[key,:]
            P_par = np.dot(Xtemp,np.dot(np.linalg.inv(np.dot(np.transpose(Xtemp),Xtemp)), np.transpose(Xtemp)))
            P_perp = np.identity(nPools) - P_par
            final_mat += np.multiply(P_perp,np.outer(y,y))

        #STEP 5 - Find and return the eigenvector with the smallest eigenvalue
            #Check: All entries should have the same sign! If negative, flip.

        eigenvalues, eigenvectors = np.linalg.eig(final_mat)
        gamma_vec = eigenvectors[:,np.argmin(eigenvalues)]/eigenvectors[:,np.argmin(eigenvalues)].mean()
        if(np.all(gamma_vec)>0):
            return gamma_vec
        else:
            raise Exception("Negative entries in eigenvector")
        
      
    def normalizePools_L2(self, stoppingCriteria, r2cutoff=0.8, maxPreys=20, verbose=False, rankType="mean", rankN=3):
        dm = self.getDesignMatrix()
        nPools = dm.shape[0]
        #Step 1: Initialize trivial normalization vector
    
        curr_norm = np.ones(nPools)
       
        #Step 2: Loop:
            #Step 2a - Create a new datset (temporary) using curr_norm
            #Step 2b - Run computePoolNomalizationVector_L2 and multiply normalizationVector by the output.
            #Step 3c - Terminate after a set number of iterations, or when teh change is small
        if verbose:
           print("> normalizePools_L2:")
        for i in range(20):
            
            #apply rescaling to original raw signal matrix 
            gamma = self.rescalePools(curr_norm).normalizeByRankedValue(rankN, type=rankType).computePoolNomalizationVector_L2(r2cutoff, maxPreys, verbose=verbose)
            curr_norm = curr_norm * gamma 
            diffSquareNorm = np.sum((gamma-1)**2)
            
            if verbose:
                print("    Marginal Rescaling: %s"%",".join(["%.3f"%vi for vi in gamma]))
                print("    Sum of rescaling values: %f"%diffSquareNorm)
                
            if diffSquareNorm < stoppingCriteria:
                break
                


        if verbose:
            print("  Final Rescaling: %s"%",".join(["%.3f"%vi for vi in curr_norm]))

        return self.rescalePools(curr_norm)
            
        #Step 3: Return the final normalized dataset.
        #pass


        #PRE-PROCESSING CODE

    #######################################
    # FUNCTIONS FOR FILTERING THE DATASET #
    #######################################

    #Filters the rows in the sitnal matrix, returns an new matrix
    def filterDataset(self, rowsKeep): 
        newSignal = self.signalMatrix[rowsKeep,:]
        newPreys = [ self.preyProteins[r] for r in rowsKeep ]
        
        return pooledDataset(self.mixingMatrix, newSignal, self.abProteins, newPreys, self.intercept, self.interceptNames, self.interceptOption, customIntercept=self.customIntercep)


    def filterMissingValues(self, minnonZeros, verbose=False):
        #drop the row from signalMatrix which has > maxZeros number of zero entries
        #if verbose is true, output how many rows dropped
        rowsKeep = []
        for i in range(self.signalMatrix.shape[0]):
            if(np.count_nonzero(self.signalMatrix[i]) >= minnonZeros):
                rowsKeep.append(i)
        
        if verbose:
            print("> filterMissingValues: %d entries were removed, %d remaining."%(len(self.signalMatrix) - len(rowsKeep), len(rowsKeep)))
            
        return self.filterDataset(rowsKeep)

 
    def filterByRankedValue(self, threshold, n, type="mean", verbose=False): #[“mean” or “value”])
        #calls computeRankedValue and filters out rows 
        # where the value is below a threshold, returns ad new pooledDataset. 
        preyList = self.computeRankedValue(n, type)
        rowsKeep = []
        for i in range(len(preyList)):
            if (preyList[i] > threshold):
                rowsKeep.append(i)
        
        if verbose:
            print("> filterByRankedValue: %d entries were removed, %d remaining."%(len(self.signalMatrix) - len(rowsKeep), len(rowsKeep)))
            
        return self.filterDataset(rowsKeep)

    def filterForBaits(self, verbose=False):
        """Filters the dataset to only keep rows corresponding to the baits."""
        rowsKeep = []
        for ab in self.abProteins:
            rowsKeep = rowsKeep + ab.findProteinMatches(self.preyProteins)
        
        if verbose:
            print("> filerForBaits: %d entries were removed, %d remaining. "%(len(self.signalMatrix) - len(rowsKeep)% len(rowsKeep)))

        return self.filterDataset(rowsKeep)
    

    def findMissingBaits(self):
        """Finds the bait proteins which were not detected in the experiment"""
        """i.e. finds the proteins present in the columns of the mixing matrix, but not found in the rows of the signal matrix"""

        missingBaits = []
        interceptCount = 0
        for i in range(len(self.abProteins)):
            if (self.abProteins[i].findProteinMatches(self.preyProteins)):
                continue
            elif("BG:" in self.abProteins[i].gene_symbol):
                interceptCount+=1
                continue
            else:
                missingBaits.append((i-interceptCount,self.abProteins[i]))
        return missingBaits
    
    def includeMissingBaits(self):
        """adds bait proteins missing from signal matrix but present in mixing matrix"""
        """they are added as dummy rows in the signal matrix with NaN signal value"""
        """WARNING: this only works if we have already filterByBait"""
        missingBaits = self.findMissingBaits()
        newSignal = np.copy(self.signalMatrix)
        newPreys = self.preyProteins.copy()
        for i in range(len(missingBaits)):
            newSignal = np.insert(newSignal, missingBaits[i][0], np.nan*np.ones(self.nPools), axis=0)
            newPreys.insert(missingBaits[i][0], missingBaits[i][1])
        return pooledDataset(self.mixingMatrix, newSignal, self.abProteins, newPreys, self.intercept, self.interceptNames, self.interceptOption, customIntercept=self.customIntercep)
   
        
                                  
    ######################################
    # FUNCTIONS FOR PLOTTING THE DATASET #
    ######################################


    def distributionOfRankedValue(self,n, type="mean"): #[“mean” or “value”])
        #calls computeRankedValue and plots the the distribution (probably using a density histogram 
        # on a log scale - this would be used to adjust the threshold.
        preyList =   self.computeRankedValue(n,type)
        modList = [val for val in preyList if val !=0]
        plt.hist(np.log10(modList), bins=10)
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.show()

    
    def plotMixingAndSignalForList(self, colors=None, baitsToPlot=None):
        """Plotting function for figure 2 specifically, with different colors, and list of baits"""
        """Also denotes correlation coefficient for each bait"""

      #If no baits specified, plot them all
        if baitsToPlot is None:
            baitsToPlot = self.abProteins  
      #Computes the mean
        unnormMean = np.average(self.signalMatrix, axis=0)
        bgSignal = unnormMean/np.mean(unnormMean)
        
        #Normalize mixing matrix columns
        B = np.zeros((np.shape(self.mixingMatrix)[0], np.shape(self.mixingMatrix)[1]))
        for i in range(np.shape(self.mixingMatrix)[1]):
            B[:,i] = self.mixingMatrix[:,i]/np.average(self.mixingMatrix, axis=1)

        n_cols = 3
        n_rows = int(np.ceil(1.*len(baitsToPlot)/n_cols)) 

        #extract bait entries from signal matrix (normalized)
        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True) 
        for i in range(n_cols*n_rows):
            row, col = divmod(i, n_cols)
            ax = axes[row, col]
            if i<len(baitsToPlot):
                baitMatches = baitsToPlot[i].findProteinMatches(self.preyProteins)
                if len(baitMatches)>0:                
                    vSignal = self.signalMatrix[baitMatches[0],:]
                    ax.plot(vSignal / np.average(vSignal), '-o', color=colors['pool'], label='Pool Signal of the Bait Protein')
                
                #baits to plot will always be present in mixing matrix
                abMatches = baitsToPlot[i].findProteinMatches(self.abProteins)
                ax.plot(B[:, abMatches[0]], '-o', color = colors['mixing'], label='Pool-wise Ab Mixing Profile')
                if(len(baitMatches)>0):
                    corr_coef = np.corrcoef(vSignal, B[:, abMatches[0]])[0, 1]
                    ax.text(0.05, 0.95, f'r = {corr_coef:.2f}', transform=ax.transAxes, verticalalignment='top', horizontalalignment='left', fontsize=10, bbox=dict(facecolor='white', alpha=0.0, edgecolor='none'))
                ax.set_title(baitsToPlot[i].gene_symbol, fontsize=11)
                ax.tick_params(labelsize=7)
                if row == n_rows - 1:
                    ax.set_xlabel('Pool', fontsize=8)
                if col == 0:
                    ax.set_ylabel('Normalized MS intensity', fontsize=6)

            else:
                #Removes unused axes.
                fig.delaxes(ax)

        (f_height, f_width) = (n_rows, n_cols) 
        fig.set_size_inches(f_width*6, f_height*3)
        fig.legend(['Pool Signal of the Bait Protein', 'Pool-wise Ab Mixing Profile'], loc='lower center', ncol=2)
    
        plt.show()
        
        
    def plotMixingAndSignalLines(self):
        
        #Computes the mean
        unnormMean = np.average(self.signalMatrix, axis=0)
        bgSignal = unnormMean/np.mean(unnormMean)

        #Normalize mixing matrix columns
        B = np.zeros((np.shape(self.mixingMatrix)[0], np.shape(self.mixingMatrix)[1]))
        for i in range(np.shape(self.mixingMatrix)[1]):
            B[:,i] = self.mixingMatrix[:,i]/np.average(self.mixingMatrix, axis=1)

        n_cols = 5
        n_rows = int(np.ceil(1.*len(self.abProteins)/n_cols))

        #extract bait entries from signal matrix (normalized)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 12), sharex=True, sharey=True)

        for i in range(n_cols*n_rows):
            row, col = divmod(i, n_cols)
            ax = axes[row, col]
            if i<len(self.abProteins):
                baitMatches = self.abProteins[i].findProteinMatches(self.preyProteins)
                if len(baitMatches)>0:                
                    vSignal = self.signalMatrix[baitMatches[0],:]
                    ax.plot(vSignal / np.average(vSignal), '-o', label='Pool Signal')
                
                ax.plot(B[:, i], '-o', label='Mixing matrix')
                ax.plot(bgSignal, '--', label='Background')
                ax.set_title(self.abProteins[i].gene_symbol, fontsize=11)
                ax.tick_params(labelsize=7)
                #ax.grid(True)
                if row == n_rows - 1:
                    ax.set_xlabel('Index', fontsize=8)
                if col == 0:
                    ax.set_ylabel('Value', fontsize=8)

            else:
                print("Not plotting %s"%self.abProteins[i].gene_symbol)
                #Removes unused axes.
                fig.delaxes(ax)

        fig.tight_layout()
        fig.legend(['Pool signal', 'Mixing matrix', 'Background'], loc='upper center', ncol=2)
    
        plt.show()