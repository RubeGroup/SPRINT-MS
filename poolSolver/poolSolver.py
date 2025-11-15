from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import nnls
from scipy.optimize import fmin_cobyla, minimize

from poolData.protein import protein
from poolData.proteinProteinMatrix import proteinProteinMatrix

class proteinSolver(ABC):
    """Abstract Class For Solving A Protein"""

    design_matrix = None
    nIntercepts = None
    
    def __init__(self,design_matrix, nIntercepts):
        self.design_matrix = design_matrix
        self.nIntercepts = nIntercepts

    @abstractmethod
    def solveProtein(self, arg1):
        pass

class poolSolver(ABC):
    #pooledDataset object 
    data = None

    def __init__(self, data):
        self.data = data

    @abstractmethod
    def solveExperiment(self):
        pass

class parallelPoolSolver(poolSolver):
    #ProteinSolver object
    ps = None


    def __init__(self, data, ps):
        super().__init__(data)
        self.ps=ps

    def solveExperiment(self):
        """Loops over experiments,runs ps.solveProtein, returns a output"""
        signalMatrix = self.data.signalMatrix
        totalCoefficients = self.data.mixingMatrix.shape[1]
        nProteins = signalMatrix.shape[0]
        betaMatrix = np.zeros((nProteins, self.ps.design_matrix.shape[1]))
        FMatrix = None
        for i in range(nProteins):
            sol = self.ps.solveProtein(signalMatrix[i])
            betaMatrix[i] = sol["beta"]
            if "F" in sol:
                if FMatrix is None:
                    FMatrix = np.zeros((nProteins, self.ps.design_matrix.shape[1]))
                FMatrix[i] = sol["F"]
                
        #return beta  
        
        
        return proteinProteinMatrix(betaMatrix, self.data.preyProteins,  self.data.interceptNames + self.data.abProteins, None, sigMatrix=FMatrix)
    #return both intercept and interactions separately  
    


class correlationSolver(proteinSolver):
    """Solver that computes the correlation coefficients between the design matrix and the measured signal."""
    def __init__(self, trainingData):
        super().__init__(trainingData.getDesignMatrix(), trainingData.nIntercepts)     

    def solveProtein(self, yMeasured):
        
        #Standardizes the design matrix
        stds = np.std(self.design_matrix, axis=0)
        stds[stds == 0] = 1  # Prevent division by zero for zero-variance columns
        dm_standardized = (self.design_matrix - np.mean(self.design_matrix, axis=0)) / stds

        #Standardizes the the signal
        row_stds = np.std(yMeasured)
        # Prevent division by zero for zero-variance rows
        if row_stds==0:
            row_stds=1
        y_standardized = (yMeasured - np.mean(yMeasured)) / row_stds

        #Computes the correlation coefficients
        result = y_standardized @ dm_standardized / len(yMeasured)
        return {"beta": result, "F": None}

 
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
        #print(f'rnorm = {rnorm}, nNonZero={nNonZero}, dmRows={dmRows}, dmCols={dmCols}')
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
                #print(f"F ={F[i]}, bestModelRSS={bestModelRSS}, exclusionRSS={exclusionRSS}, masked={maskedMatrix.shape}")
        
        return {"beta": beta, "F": F}
    
class leastSquaresSolver(proteinSolver):
    """Solver that computes the least squares solution to the design matrix."""

    def __init__(self, trainingData):
        super().__init__(trainingData.getDesignMatrix(), trainingData.nIntercepts)

    def solveProtein(self, yMeasured):
        #Computes the least squares solution
        res = minimize(lambda x: np.sum((self.design_matrix @ x - yMeasured) ** 2), 
                       np.zeros(self.design_matrix.shape[1]), 
                       method='L-BFGS-B')
        return {"beta": res.x, "F": None}
    
class pseudoInverseLeastSquaresSolver(proteinSolver):
    """Solver that computes the pseudo-inverse least squares solution to the design matrix."""
  
    def __init__(self, trainingData):
        super().__init__(trainingData.getDesignMatrix(), trainingData.nIntercepts)

    def solveProtein(self, yMeasured):
        #Computes the pseudo-inverse least squares solution
        pinv_design_matrix = np.linalg.pinv(self.design_matrix)
        beta = pinv_design_matrix @ yMeasured
        return {"beta": beta, "F": None}

