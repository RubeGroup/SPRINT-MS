from .poolSolver import proteinSolver, poolSolver
from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import nnls
from itertools import combinations


class subsetSelectionProteinSolvers(proteinSolver):
    """Abstract Class for Solving a protein using forward or best subset selection"""

    model = None
    #data = None
    nIntercepts = None

    def __init__(self, trainingData, stopping_criteria, model):
        super().__init__(trainingData.getDesignMatrix(), trainingData.nIntercepts)      #attribute of parent protein solver class
        self.model = model
        self.stopping_criteria = stopping_criteria  #how strictly does the algorithm compare two consecutive model
        

    def R2(self,y_pred, y_obs):
  
        y_avg = np.mean(y_obs)
        tss = np.sum((y_obs - y_avg )**2)
        rss = np.sum((y_obs - y_pred)**2)   
        if tss > 1e-10:
            return (1 - rss/tss)
        else:
            return 0.0
        
    #Identifies the best model.
    def selectBestModel(self, modelList, tss):
        #By default, the baseline model is the best.
        bestModel = modelList[0]
        n = self.design_matrix.shape[0]
        
        #Loops over all non-trivial models.
        for i in range(1,len(modelList)):
        
            if self.model == 'abs': #Stop once R^2 reaches the stopping criteria
                #Terminates if R2
                if tss<1e-10 or (1-modelList[i]["RSS"]/tss) > self.stopping_criteria:
                    #Quitting without saving
                    break
                else:
                    #Saving and continuing
                    bestModel = modelList[i]

            elif self.model == 'diff': #stop once the improvement increase in R^2 
                #Trivial model returned if TSS is too small.
                if tss<1e-10:
                    break
                
                deltaR2 = (modelList[i-1]["RSS"]-modelList[i]["RSS"])/tss
                if deltaR2< (self.stopping_criteria):
                    #Quitting without saving
                    break
                else:
                    #Saving and continuing
                    bestModel = modelList[i]
            
            elif self.model == 'diff_global': #Returns the largest model for which the delta R2 exceeds stopping_criteria value
                
                #Trivial model returned if TSS is too small.
                if tss<1e-10:
                    break

                deltaR2 = (modelList[i-1]["RSS"]-modelList[i]["RSS"])/tss
                if deltaR2 > self.stopping_criteria:
                    #Saving
                    bestModel = modelList[i]

            elif self.model == 'F': #Stops once the F-value falls below a threshold.
                #Trivial model returned if TSS is too small.
                if tss<1e-10:
                    break
                
                #modify to add the number of predictors difference from modelList
                (RSS1, RSS2) = ( modelList[i-1]["RSS"], modelList[i]["RSS"] )
                (df1,  df2)  = ( modelList[i-1]["p"],   modelList[i]["p"] )
                F = ((RSS1-RSS2)/(df2-df1))/(RSS2/(n - df2))
                                
                #delRSS = modelList[i-1]["RSS"]-modelList[i]["RSS"]
                #delPreds = len(modelList[i]["indices"]) - len(modelList[i-1]["indices"])
                #F = (delRSS/delPreds)/(modelList[i]["RSS"]/(n - len(modelList[i]["indices"])))
                #(())/(modelList[i]["RSS"]/(modelList[i]["A"].shape[0] - modelList[i]["A"].shape[1]))
                if F< (self.stopping_criteria):
                    #Quitting without saving
                    break
                else:
                    #Saving and continuing
                    bestModel = modelList[i]
                    
            elif self.model == 'F_global': #Returns the largest model for which the F-value exceeds the stopping_criteria value
                #Trivial model returned if TSS is too small.
                if tss<1e-10:
                    break

                #delRSS = modelList[i-1]["RSS"]-modelList[i]["RSS"]
                #delPreds = len(modelList[i]["indices"]) - len(modelList[i-1]["indices"])
                #F = (delRSS/delPreds)/(modelList[i]["RSS"]/(n - len(modelList[i]["indices"])))    
                (RSS1, RSS2) = ( modelList[i-1]["RSS"], modelList[i]["RSS"] )
                (df1,  df2)  = ( modelList[i-1]["p"],   modelList[i]["p"] )
                F = ((RSS1-RSS2)/(df2-df1))/(RSS2/(n - df2))
                if F > self.stopping_criteria:
                    #Saving
                    bestModel = modelList[i]

            else:
                raise Exception("Invalid stopping criteria: %s"%self.model)
        return bestModel
    

        

class NNLSForwardSelectionProteinSolver(subsetSelectionProteinSolvers):

    maxSize = None

    def __init__(self, trainingData, stopping_criteria, model, maxSize=5):
        self.maxSize = maxSize
        super().__init__(trainingData, stopping_criteria, model)      #attribute of parent protein solver class


    def solveProtein(self, yMeasured):

        A = self.design_matrix #pdata.getDesignMatrix()
        curr_indices = list(range(self.nIntercepts))

        curr_A = A[:,curr_indices] #np.concatenate((A[:,:1], np.zeros((A.shape[0], A.shape[1]-1))),axis=1) 
        #this should have the same shape as A but the missing columns should be filled with zeros. 
        b = yMeasured
        x            = nnls(curr_A,b)[0]
        temp_yHat    = np.dot(curr_A,x)
        curr_rss   = np.sum((temp_yHat-b)**2)
                
        tss = np.sum((b - np.mean(b))**2)

        #Computes the initial R2.
        curr_yHat = np.dot(curr_A,nnls(curr_A,b)[0])
        #curr_R2   = self.R2(curr_yHat, b)
        
        #Saves performance of baseline model
        modelList = [{
            "RSS": curr_rss,
            "indices": curr_indices,
            #"A": curr_A,
            "p": 0
            }]
        
        #Set of all unused predictors
        unused_indices = set(range(A.shape[1])) - set(curr_indices)
        
        #Adds one predictor at a time. 
        while len(unused_indices)>0:
#        for i in range(A.shape[1]-1):
            rss_dict = {}
            
            #Tries all possible additions to the current model
            for kTest in unused_indices:

                #Fits the model and predictions, saves relevant information.
                temp_A          = np.concatenate((curr_A,A[:,kTest][:,np.newaxis]), axis=1) #np.zeros((A.shape[0], A.shape[1]-1))),axis=1) #A[:,k]
                x               = nnls(temp_A,b)[0]
                temp_yHat       = np.dot(temp_A,x)
                rss_dict[kTest] = np.sum((temp_yHat-b)**2)
                
            #Identifies the model with the smallest RSS
            new_predictor  = int(min(rss_dict, key =  rss_dict.get))
            curr_indices   = curr_indices       + [new_predictor]
            unused_indices = unused_indices - set([new_predictor])
            curr_rss  = rss_dict[new_predictor]
            curr_A    = np.concatenate((curr_A,A[:,new_predictor][:,np.newaxis]), axis=1)
            #Saves information about the best model
            modelList.append({
                "RSS": curr_rss, 
                "indices": curr_indices,
                #"A": curr_A,
                "p": (len(curr_indices) - self.nIntercepts)
                })

            #Stops if the set of predictors has reached the maximum size.
            if self.maxSize!=None and len(curr_indices) >=  self.nIntercepts + self.maxSize:
                break
            
        bestModel = self.selectBestModel(modelList, tss)
        
        #Creates the final coefficient vector.
        final_indices = bestModel["indices"]
        coeffs = np.zeros(A.shape[1])
        np.put(coeffs, bestModel["indices"], nnls(A[:,bestModel["indices"]],b)[0])

        return {"beta": coeffs, "F": None}

 
 
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
        #bestExclusionRSS = [ None for i in range(X.shape[1])]
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
            #new_predictor = int(min(rss_dict, key =  rss_dict.get))
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
       