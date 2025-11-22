from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import nnls
from scipy.optimize import fmin_cobyla, minimize

from poolData.protein import protein
from poolData.proteinProteinMatrix import proteinProteinMatrix
from poolSolver.poolSolver import proteinSolver


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