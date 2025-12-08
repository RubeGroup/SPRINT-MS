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