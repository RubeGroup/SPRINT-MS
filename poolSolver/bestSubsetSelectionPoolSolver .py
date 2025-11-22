from .poolSolver import proteinSolver, poolSolver
from .subsetSelectionProteinSolvers import subsetSelectionProteinSolvers
from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import nnls
from itertools import combinations





class NNLSBestSubsetSelectionProteinSolver(subsetSelectionProteinSolvers):

    projection_matrices  = None
    orthogonal_matrices  = None
    predictor_list       = None
    baitExclusionIndices = None
    matrixList           = None
    maxSize              = None

    def __init__(self, dm, stopping_criteria, model, maxSize=3):
        super().__init__(dm, stopping_criteria, model)      #attribute of parent protein solver class
        self.maxSize = maxSize                              
        self.matrixList = self.computeMatrices()       #precompute all projection, orthogonal
                                                                #matrices, up till maxSize number of predictors
        
        
        # For each number of variables p, and for each the possible set of variables,
        # pre-compute (X^TX)^-1X^T, and X(X^TX)^-1X^T.
        # (X^TX)^-1X^T:   [#added predictors][model inddex][iAB][iPool]
        # X(X^TX)^-1X^T:  [#added predictors][model inddex][iPool][iPool]

    #Lists all combinations
    def listAllCombinations(n):
        out = [[]]
        for i in range(n):
            out = out + [ elem + [i] for elem in out ]
        return out

    #Lists all combinations containing k elements
    def listCombinations(n,k):
        if k>0:
            out = [ [i] for i in range(n)]
            while len(out[0]) < k:
                out = [comb + [iNew]  for comb in out for iNew in range(comb[-1]+1,n)]
        else:
            out = []
        return out
    
    def computeMatrices(self):
        
        X = self.design_matrix
        matrixList = {}
        #Lists the indices of the background and bait columns
        bgIndices   = tuple(range(self.nIntercepts))
        baitIndices = tuple(range(self.nIntercepts, X.shape[1]))
        self.projection_matrices = {}
        self.orthogonal_matrices = {}
        self.predictor_list      = {}
        self.predictor_list      = {}
        self.baitExclusionIndices= {}
        
        bgCombinations = NNLSBestSubsetSelectionProteinSolver.listAllCombinations(self.nIntercepts)
        
        #Loops over the number of predictors
        for k in range(0, self.maxSize+1):
            self.projection_matrices[k] = []
            self.orthogonal_matrices[k] = []
            self.predictor_list[k]      = []
            self.baitExclusionIndices[k]= [ [] for i in range(X.shape[1])]
                            
                #Lits all combinations of k baits
            baitCombinations = [ [baitIndices[i]for i in comb] for comb in NNLSBestSubsetSelectionProteinSolver.listCombinations(len(baitIndices),k) ]
                
            #Lists all combinations
            if len(baitCombinations)>0:
                if len(bgCombinations)>0:
                    allCombs = [bgComb+baitComb for bgComb in bgCombinations for baitComb in baitCombinations ]
                else:
                    allCombs = baitCombinations
            else:
                if len(bgCombinations)>0:
                    allCombs = bgCombinations
                else:
                    allCombs = []
                
            #Loops over all cobminations of background and baits
            for i, comb in enumerate(allCombs):
                #Skips the empty combination 
                if len(comb)==0:
                    continue
                #Computes the orthogonal and projection matrices
                curr_X = X[:,comb]
                XtX = np.dot(curr_X.T, curr_X)
                try:
                    XtX_inv = np.linalg.inv(XtX)
                except np.linalg.LinAlgError:
                    # In case of singular matrix, skip this combination
                    raise("Singular matrix")

                projection_matrix = np.dot(XtX_inv, curr_X.T)
                orthogonal_matrix = np.identity(X.shape[0]) - np.dot(curr_X, projection_matrix)
                        
                #Saves the data to the lists
                self.predictor_list[k].append(comb)
                self.projection_matrices[k].append(projection_matrix)
                self.orthogonal_matrices[k].append(orthogonal_matrix)
                    
                #Identifies the predictors that are not present in this combination
                for l in  range(X.shape[1]):
                    if l not in comb:
                        self.baitExclusionIndices[k][l].append(i)
                            
            #Converts the lists to sets
            for l in  range(X.shape[1]):
                self.baitExclusionIndices[k][l] = set(self.baitExclusionIndices[k][l])

    def solveProtein(self, yMeasured, fullOutput=False):
        
        X = self.design_matrix
        b = yMeasured
        tss = np.sum((b - np.mean(b))**2)
        
        #List containing the best model for each size (number of baits)
        modelList = []
        
        #List containing the smallest RSS when each predictor is held out.
       
        bestExclusionRSS = [ tss for i in range(X.shape[1])]
        
        #List containing the re-fitted models when each predictor is held out.
        bestExclusionBetas = []
        
        #Loops over the number of predictors
        for k in range(0, self.maxSize+1):
            #Extracts the relevant lists of combinations and matrices
            pms   = self.projection_matrices[k]
            oms   = self.orthogonal_matrices[k]
            

            #Comptes the RSS for all combinations with non-negative betas.
            betaList = [np.dot(pm, yMeasured) for pm in pms ]
            iFiltered = [  i for  i,beta in enumerate(betaList) if (beta>=0).all() ]   
            
            #Contiues if all models are invalid (meaning the betas are negative ie not non-negative)
            if len(iFiltered)==0:
                print("Skipping %d"%k)
                continue         
            RSSList = [np.dot(yMeasured.T, np.dot(yMeasured,oms[i])) for i in iFiltered]
            
            #Identifies the model with the smallest RSS
           
            #Handle exception when predlist is empty, check if ALL betas are empty
            try:
                min_filtered_index = RSSList.index(min(RSSList))
                min_index          = iFiltered[min_filtered_index]
            except ValueError as e:
                if "min() arg is an empty sequence" in str(e):
                    raise("Empty list, all coefficients negative")
                
            best_RSS        = float(RSSList[min_filtered_index])
            best_predictors = self.predictor_list[k][min_index]
                
            #Identifies the smallest possible RSS for each excluded predictor.
            RSS_exclusion = []
            beta_exclusion = []
            #Loops over the predictors
            for l in range(X.shape[1]):
                previousBest = bestExclusionRSS[l]
                #Checks if the predictor l is in the best model
                if l in best_predictors:
                    #Extracts all RSS for the models not included predictor l
                    tempRSS = [RSSList[i] for (i, iFilt) in enumerate(iFiltered) if iFilt in self.baitExclusionIndices[k][l]]
                    if len(tempRSS) > 0:
                        currentBest = float(min(tempRSS))
                        
                        #Saves the index of the best model hold-out model
                        if fullOutput:
                            min_exclusion_index = [iFilt for (i, iFilt) in enumerate(iFiltered) if iFilt in self.baitExclusionIndices[k][l] and float(RSSList[i])== currentBest][0]
                            min_exclusion_betas = betaList[min_exclusion_index]
                            min_exclusion_cols  = self.predictor_list[k][min_exclusion_index]
                            beta_vectors=np.zeros(X.shape[1])
                            for (i,ic) in enumerate(min_exclusion_cols):
                                beta_vectors[ic] = min_exclusion_betas[i]
                            beta_exclusion.append(beta_vectors)
                        else:
                            beta_exclusion.append(None)
                    else:
                        currentBest = None
                        if fullOutput:
                            beta_exclusion.append(None)
                else:
                    currentBest = best_RSS
            
                #Updates the best RSS
                if previousBest is None:
                    newBest = currentBest
                else:
                    if currentBest is None:
                        newBest = previousBest
                    else:
                        newBest = min(currentBest,previousBest)
                #Saves the new best
                RSS_exclusion.append(newBest)
                bestExclusionRSS[l] = newBest
                

            #Saves information about the best model
            modelList.append({
                "RSS":     best_RSS,
                "indices": best_predictors, 
                "p": k,
                "betas":   betaList[min_index],
                "RSS_exclusion": RSS_exclusion,
                "beta_exclusion": beta_exclusion
                })
            
        #Creates the final coefficient vector.
        bestModel = self.selectBestModel(modelList, tss)
        coeffs = np.zeros(X.shape[1])
        np.put(coeffs, 
               bestModel["indices"], 
               bestModel["betas"])
        
        F =  (np.array(bestModel["RSS_exclusion"])-bestModel["RSS"])/(bestModel["RSS"]/(len(X)-bestModel["p"]))
        
        if fullOutput:
            return {"beta": coeffs, "F": F, "modelList": modelList, "indices": bestModel["indices"]}
        else:
            return {"beta": coeffs, "F": F}