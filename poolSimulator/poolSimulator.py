from abc import ABC, abstractmethod
import numpy as np
from poolData.proteinProteinMatrix import proteinProteinMatrix
from poolData.protein import protein
#from scipy.optimize import nnls
#from scipy.optimize import fmin_cobyla, minimize
from sklearn.metrics import roc_curve, auc
import numpy.random as random
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

from poolData.pooledDataset import pooledDataset
import matplotlib.pyplot as plt
#from poolData.proteinProteinMatrix import proteinProteinMatrix

class poolSimulator():
    """class for simulating a pool experiment."""

    mixingMatrix = None
    baitProteins = None
    nBaits = None
    nPools = None
    groundTruthPPIs = None
    syntheticData = None

    def __init__(self,mixingMatrix, baitProteins):
        (self.nPools, self.nBaits) = mixingMatrix.shape
        self.mixingMatrix = mixingMatrix
        self.baitProteins = [ protein("Bait %d"%(i+1), ["Bait %d"%(i+1)]) for i in range(self.nBaits)]
        
    def simulateExperiment(self, noiseSD=0.1, noiseType="normal", seed=42):
        
        rng = np.random.default_rng(seed=seed)

        """Simulate an experiment by generating a noisy PPI matrix."""
        if self.groundTruthPPIs is None:
            raise ValueError("Ground truth PPI matrix not set. Call buildDiagonalPPIMatrix or another method to set it.")
        
        if noiseType == "normal" or noiseType == "normal_nonneg":
            # Simulate noise - n
            noise = rng.normal(0, noiseSD, (self.groundTruthPPIs.nRows, self.nPools))
            signalMatrix = np.matmul(self.groundTruthPPIs.matrix, self.mixingMatrix.T)+noise   
            
            if noiseType == "normal_nonneg":
                # Ensure non-negativity
                signalMatrix = np.maximum(signalMatrix, 0)  
        else:
            raise ValueError("Unsupported noise type. Use 'normal' or 'normal_nonneg'.")
        
        self.syntheticData = pooledDataset(self.mixingMatrix, signalMatrix, self.baitProteins, self.groundTruthPPIs.rowProteins, None, None, "ones", None)
        
        return
    
    ######################################################
    # Methods for building the ground truth PPI matrices #
    ######################################################

    def buildExactPPIMatrix(self,nPreys,betaValues, ppiProb=0.7):
        """builds a PPI matrix where each prey has beta values only drawn from betaValues"""
        ppiMatrix=np.zeros((nPreys, self.nBaits))
        self.groundTruthAnnotations = []
        preyProteins = [ protein("Prey %d"%(i+1), ["Prey %d"%(i+1)]) for i in range(nPreys)]
        
        for iPrey in range(nPreys):
            if np.random.rand() < ppiProb:
                iBaits = np.random.choice(self.nBaits, size=len(betaValues), replace=False)
                ppiMatrix[iPrey,iBaits] = betaValues
                self.groundTruthAnnotations.append((preyProteins[iPrey], self.baitProteins[iBait]) for iBait in iBaits)
        
        self.groundTruthPPIs = proteinProteinMatrix(ppiMatrix, preyProteins, self.baitProteins, None,  None)

    def buildDiagonalPPIMatrix(self, strength=1):
        """Builds a diagonal PPI matrix where each bait has a single prey."""
        ppiMatrix = np.identity(self.nBaits)
        preyProteins = [ protein("Bait %d"%(i+1), ["Bait %d"%(i+1)]) for i in range(self.nBaits)]
        
        self.groundTruthPPIs = proteinProteinMatrix(ppiMatrix, preyProteins, self.baitProteins, None,  None)
        self.groundTruthAnnotations = [(p,p) for p in self.baitProteins]
        
    def buildBandedPPIMatrix(self, nMod, strength=1):
        
        #Creates an empty matrix
        ppiMatrix = np.array([ [ strength if j>=i and (j-i)%nMod==0 else 0 for i in range(self.nBaits)] for j in range(self.nBaits) ])
        preyProteins = [ protein("Bait %d"%(i+1), ["Bait %d"%(i+1)]) for i in range(self.nBaits)]
        self.groundTruthPPIs = proteinProteinMatrix(ppiMatrix, preyProteins, self.baitProteins, None,  None)
        self.groundTruthAnnotations = [(p,p) for p in self.baitProteins]


    def buildDiagonalVaryingStrengthMatrix(self, strength=1):
        
        #Creates an empty matrix
        ppiMatrix = np.array([ [ strength*i/(self.nBaits-1) if j==i  else 0 for i in range(self.nBaits)] for j in range(self.nBaits) ])
        
        preyProteins = [ protein("Bait %d"%(i+1), ["Bait %d"%(i+1)]) for i in range(self.nBaits)]
        self.groundTruthPPIs = proteinProteinMatrix(ppiMatrix, preyProteins, self.baitProteins, None,  None)
        self.groundTruthAnnotations = [(p,p) for p in self.baitProteins]



    def buildRandomPPIMatrix(self, nPreys, nPPIs, strength=1):
        """Builds a random PPI matrix with a specified number of preys and PPIs."""
        ppiMatrix=np.zeros((nPreys, self.nBaits))
        self.groundTruthAnnotations = []
        preyProteins = [ protein("Prey %d"%(i+1), ["Prey %d"%(i+1)]) for i in range(nPreys)]
        
        for iPrey in range(nPreys):
            for iBait in random.choice(range(0,self.nBaits),nPPIs, replace=False):
                ppiMatrix[iPrey, iBait] = strength
                self.groundTruthAnnotations.append((preyProteins[iPrey], self.baitProteins[iBait]))
        
        self.groundTruthPPIs = proteinProteinMatrix(ppiMatrix, preyProteins, self.baitProteins, None,  None)

    def buildBernoulliPPIMatrix(self, nPreys, p, strength=1):
        """Builds a random PPI matrix with a specified number of preys and PPIs."""
        

        ppiMatrix = np.random.binomial(1, p, (nPreys, self.nBaits))*strength

        self.groundTruthAnnotations = []
        preyProteins = [ protein("Prey %d"%(i+1), ["Prey %d"%(i+1)]) for i in range(nPreys)]
        
        for iPrey in range(nPreys):
            for iBait in range(self.nBaits):
                if ppiMatrix[iPrey, iBait] > 0:
                    self.groundTruthAnnotations.append((preyProteins[iPrey], self.baitProteins[iBait]))
        
        self.groundTruthPPIs = proteinProteinMatrix(ppiMatrix, preyProteins, self.baitProteins, None,  None)


    def splitInferredPPIsByTruth(self, inferredPPIs, splitSigMatrix=False):
        """Splits the inferred PPI matrix into two lists: one for true PPIs and one for non-PPIs."""
        truth = self.groundTruthPPIs
        
        #Creates mapping from the columns in the ground-truth matrix to the inferred matrix. 
        colsTrueToInferred = dict([ (iTruth, iPred) for (iTruth,pTruth) in enumerate(truth.columnProteins) for (iPred, pPred) in enumerate(inferredPPIs.columnProteins) if pTruth.checkMatchQ(pPred) ])
        (ppiList, nonPpiList) = ([], [])
        for iCol in range(truth.nCols):
            for iRow in range(truth.nRows):
                
                if splitSigMatrix:
                    value = inferredPPIs.sigMatrix[iRow, colsTrueToInferred[iCol]]
                else:
                    value = inferredPPIs.matrix[iRow, colsTrueToInferred[iCol]]
                    
                if truth.matrix[iRow, iCol] > 0:
                    ppiList.append(value)
                else:
                    nonPpiList.append(value)
                
        return (np.array(ppiList), np.array(nonPpiList))
    
    def plotInferredPPIHistogram(self, inferredPPIs, plotType="histogram", nBins=10, plotSig=False):
        
        """Plots the histogram or eCDF of the inferred PPI values compared to the ground truth PPI values."""
        
        (ppiList, nonPpiList) = self.splitInferredPPIsByTruth(inferredPPIs, splitSigMatrix=plotSig)


        cPPI = (0.0901961, 0.172549, 0.313725)
        cNonPPI = (0.568627, 0.223529, 0.423529)
        
        if plotType == "eCDF":
            plt.step(np.sort(ppiList),  np.arange(1, len(ppiList)+1) / len(ppiList), where='post')
            plt.step(np.sort(nonPpiList),  np.arange(1, len(nonPpiList)+1) / len(nonPpiList), where='post')
            plt.ylabel('eCDF')
            plt.legend(['Truth: PPI', 'Truth: Non-PPI'])
            
            if plotSig:            
                plt.xlabel('F')
                plt.title('Distribution of inferred F-values')
            else:
                plt.title('Distribution of inferred beta-values')
                plt.xlabel('Value')

        elif  plotType == "histogram":
            plt.hist([ppiList, nonPpiList], bins=nBins, color = [cPPI, cNonPPI], density=True)  # 'bins' sets the number of bars
            plt.ylabel('Probability Density')
            plt.legend(['Truth: PPI', 'Truth: Non-PPI'])
            
            if plotSig:            
                plt.xlabel('F')
                plt.title('Distribution of inferred F-values')
            else:
                plt.title('Distribution of inferred beta-values')
                plt.xlabel('Value')
            
        elif plotType == "ROC":
            # Calculate ROC curve            
            y_true = np.concatenate((np.ones(len(ppiList)), np.zeros(len(nonPpiList))))
            y_scores = np.concatenate((ppiList, nonPpiList))
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            roc_auc = auc(fpr, tpr)

            plt.plot(fpr, tpr, color='blue', label=f'ROC curve (AUC = {roc_auc:.2f})')
            plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
            
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            
            if plotSig:
                plt.title('ROC Curve (using inferred F-values)')
            else:
                plt.title('ROC Curve (using inferred beta-values)')
            plt.legend(loc='lower right')


        elif plotType == "PRC":
            y_true = np.concatenate((np.ones(len(ppiList)), np.zeros(len(nonPpiList))))
            y_scores = np.concatenate((ppiList, nonPpiList))
            precision, recall, _ = precision_recall_curve(y_true, y_scores)
            average_precision = average_precision_score(y_true, y_scores)

            plt.plot(recall, precision, color='blue', label=f'PRC (AP = {average_precision:.2f})')
            plt.xlabel('Recall')
            plt.ylabel('Precision')

            if plotSig:
                plt.title('PRC Curve (using inferred F-values)')
            else:
                plt.title('PRC Curve (using inferred beta-values)')
            plt.legend(loc='lower right')
            
        else:
            raise ValueError("Invalid plot type. Choose from 'histogram', 'eCDF', 'ROC', or 'PRC'.")

        
        plt.show()
    