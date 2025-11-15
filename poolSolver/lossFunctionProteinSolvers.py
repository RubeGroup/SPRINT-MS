
from abc import ABC, abstractmethod
import numpy as np
from .poolSolver import proteinSolver
from scipy.optimize import fmin_cobyla, minimize

class dataLoss(ABC):

    def __init__(self,design_matrix):
        self.design_matrix = design_matrix
    
    def yHat(self,arg1):
        return np.dot(self.design_matrix, arg1)
        

    @abstractmethod
    def loss(self):
        pass

    @abstractmethod
    def grad(self):
        pass

class lossFunctionProteinSolver(proteinSolver):
    
    """Abstract Class For Solving A Protein Using Loss Function"""
    lf = None  #a log ratio loss object for the data loss function part
    interceptSize = None
    
    def __init__(self,lf):
        super().__init__(lf.design_matrix)
        self.lf = lf        #logRatioLoss object

    @abstractmethod
    def totalLoss(self):
        pass
    
    @abstractmethod
    def totalGrad(self):
        pass

    def checkFiniteDifference(self, v, eps, *args):
        #General for any loss-function based solver
        j = np.zeros(np.size(v))
        for i in range(j.size):
            v1 = np.copy(v)
            v2 = np.copy(v)
            v1[i] = v1[i] + eps
            v2[i] = v2[i] - eps
            j[i] = (self.lf.loss(v1,*args) - self.lf.loss(v2,*args))/(2*eps)
        return j
        

    @abstractmethod
    def solveProtein(self,arg1):
        pass

#create a class for involving the intercept, having a get design matrix method
class linearLoss(dataLoss):
    #define functions in here
    def __init__(self,design_matrix, epsiln):
        super().__init__(design_matrix)
        self.epsiln = epsiln #offset to avoid divde by zero 
     
    def loss(self,beta,*args):
        return (np.linalg.norm(args[0]- self.yHat(beta),2)**2)/15

    def grad(self,beta,*args):
        yh = self.yHat(beta)
        yp = args[0]
        return -2*np.dot(np.transpose(np.log(np.divide(yp+self.epsiln,yh+self.epsiln))),(self.design_matrix/yh[:,None])) 
    
class logRatioLoss(dataLoss):
    #define functions in here
    def __init__(self,design_matrix, epsiln):
        super().__init__(design_matrix)
        self.epsiln = epsiln #offset to avoid divde by zero 
     
    def loss(self,beta,*args):
        return np.linalg.norm(np.log(args[0]+self.epsiln) - np.log(self.yHat(beta)+self.epsiln),2)**2

    def grad(self,beta,*args):
        yh = self.yHat(beta)
        yp = args[0]
        return -2*np.dot(np.transpose(np.log(np.divide(yp+self.epsiln,yh+self.epsiln))),(self.design_matrix/yh[:,None])) 
	
        #return gradient contribution of the data loss    

class L1RegularizedProteinSolver(lossFunctionProteinSolver):

    def __init__(self,lf, lambdaL1, lambdaL2, method, fitLog,constraints, interceptSize):
        super().__init__(lf)
        self.lambdaL1 = lambdaL1	#L1, L2 parameters
        self.lambdaL2 = lambdaL2    
        self.fitLog = fitLog    #boolean to fit in log space (true) or not (false)
        self.method = method    #sklearn minimization algorithm
        #self.epsiln = epsiln    #add to denominator to avoid divide by zero error
        self.constraints = constraints  #boolean to use constraints or not
        self.interceptSize = interceptSize
        #self.d = d      #design matrix


    def totalLoss(self, v, *args):
        #add regularization to data loss from lf
        if self.fitLog:
            beta = np.exp(v) #np.concatenate((v[0, np.newaxis],np.exp(v[1:])))
        else:
            beta = v
        b_intercept = beta[:self.interceptSize]
        b_ppi = beta[self.interceptSize:]
        value = self.lf.loss(beta, *args) + self.lambdaL1*np.linalg.norm(b_ppi,1) + \
            self.lambdaL2*np.linalg.norm(b_ppi,2)**2
        if np.isnan(value):
            print("NaN:")
            print(v)
        return value
    
    def totalGrad(self,v,*args):
        #add gradient contribution from regularization terms
        if self.fitLog:
            beta = np.exp(v)
            return self.lf.grad(beta,*args)*beta + \
                  np.concatenate((np.zeros(1),(self.lambdaL1*np.ones(v.size-1)*beta[1:]))) + \
                    2*np.concatenate((np.zeros(1),(self.lambdaL2*beta[1:]))) #2*self.lambdaL2*
        else:
            return self.lf.grad(v,*args) + self.lambdaL1*np.ones(v.size) + 2*self.lambdaL2*v
        return grad
    
    def solveProtein(self, yMeasured):
        #optimization step for signal from one protein
        seed = np.random.rand(self.lf.design_matrix.shape[1])
        if(self.constraints):
            const = ({'type': 'ineq', 'fun': self.lf.yHat, 'args': ()})
        else:
            const = None
            
        while (np.min(np.dot(self.lf.design_matrix,seed))<= 0):
            seed = np.random.rand(self.lf.design_matrix.shape[0]+1)
        return minimize(self.totalLoss, seed, method = self.method,args=(yMeasured), 
                        jac = self.totalGrad, constraints = const,options={'maxiter':10000}).x
    
