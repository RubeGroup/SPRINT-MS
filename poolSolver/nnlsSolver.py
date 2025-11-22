from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import nnls
from scipy.optimize import fmin_cobyla, minimize

from poolData.protein import protein
from poolData.proteinProteinMatrix import proteinProteinMatrix
from poolSolver.poolSolver import proteinSolver



class nnlsSolver(proteinSolver):
    """Solver that computes the non-negative least squares solution to the design matrix."""

    def __init__(self, trainingData):
        super().__init__(trainingData.getDesignMatrix(), trainingData.nIntercepts)

    def solveProtein(self, yMeasured):
        #Computes the non-negative least squares solution
        beta, rnorm = nnls(self.design_matrix, yMeasured)
        
        #Counts the number of non-zero coefficients.
        nNonZero = np.count_nonzero(beta)
        
        #Computes the F statistic
        F = np.zeros(len(beta))
        (dmRows, dmCols) = self.design_matrix.shape

        #Computes the F statistic if a) the residual is non-zero and b) the number 
        # of non-zero parameters is smaller than than the number of data points.
        if rnorm > 0 and nNonZero < dmRows:
            for i in range(dmCols):
                #print(f"Solving {i}")
                maskedMatrix = np.delete(self.design_matrix, i, axis=1)
                beta_temp, rnorm_temp = nnls(maskedMatrix, yMeasured)
                bestModelRSS = rnorm**2
                exclusionRSS = rnorm_temp**2
                F[i] =  (exclusionRSS-bestModelRSS)/(bestModelRSS/(dmRows-nNonZero))
        
        return {"beta": beta, "F": F}