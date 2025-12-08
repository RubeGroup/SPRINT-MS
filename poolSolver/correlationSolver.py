from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import nnls
from scipy.optimize import fmin_cobyla, minimize

from poolData.protein import protein
from poolData.proteinProteinMatrix import proteinProteinMatrix
from poolSolver.parallelProteinSolver import proteinSolver



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