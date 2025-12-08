from .parallelProteinSolver import proteinSolver, poolSolver
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

             
                (RSS1, RSS2) = ( modelList[i-1]["RSS"], modelList[i]["RSS"] )
                (df1,  df2)  = ( modelList[i-1]["p"],   modelList[i]["p"] )
                F = ((RSS1-RSS2)/(df2-df1))/(RSS2/(n - df2))
                if F > self.stopping_criteria:
                    #Saving
                    bestModel = modelList[i]

            else:
                raise Exception("Invalid stopping criteria: %s"%self.model)
        return bestModel
    

        

       