import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import matplotlib.colors as colors


class proteinProteinMatrix:
    """General class for storing a data matrix with rows/columns corresponding to proteins """
    #preyProteins = rowProteins
    #abProteins = columnProteins 

    matrix = None
    sigMatrix = None
    rowProteins = None
    columnProteins = None
    nRows = None
    nCols = None
    
    def __init__(self, matrixData, preyProteins, abProteins, intlist, sigMatrix = None):
        
        
        self.matrix = matrixData
        self.sigMatrix = sigMatrix #Additional matrix storing significance information.
        self.rowProteins = preyProteins
        self.columnProteins = abProteins
        self.intlist = intlist      #an interactorList object containing all ineractions
        (self.nRows, self.nCols) = self.matrix.shape
        
        #Defines colors
        self.cMix=(0.0901961, 0.172549, 0.313725)
        self.cPPI=(0.207843, 0.478431, 0.462745)
        self.cSig=(0.835294, 0.67451, 0.223529)
        self.cNeg=(0.568627, 0.223529, 0.423529)
        self.cmapPPI = colors.LinearSegmentedColormap.from_list("", [self.cNeg,"white",self.cPPI])
        self.cmapPPI.set_bad((0.75, 0.75, 0.75), alpha=1.0)
        self.cmapSig = colors.LinearSegmentedColormap.from_list("", ["white",self.cSig])
        self.cmapSig.set_bad((0.75, 0.75, 0.75), alpha=1.0)
        self.cmapMix = colors.LinearSegmentedColormap.from_list("", ["white",self.cMix])
        

    def removeColumns(self, indices):
        """Removes the columns specified by indices from the matrix and updates the column proteins."""
        newMatrix = np.delete(self.matrix, indices, axis=1)
        newColumnProteins = [p for i,p in enumerate(self.columnProteins) if i not in indices]
        if self.sigMatrix is not None:
            newSigMatrix = np.delete(self.sigMatrix, indices, axis=1)
        else:
            newSigMatrix = None 
        
        return proteinProteinMatrix(newMatrix, self.rowProteins, newColumnProteins, self.intlist, sigMatrix=newSigMatrix)        

    def removeRows(self, indices):
        """Removes the rows specified by indices from the matrix and updates the row proteins."""
        newMatrix = np.delete(self.matrix, indices, axis=0)
        newRowProteins = [p for i,p in enumerate(self.rowProteins) if i not in indices]
        if self.sigMatrix is not None:
            newSigMatrix = np.delete(self.sigMatrix, indices, axis=0)
        else:
            newSigMatrix = None 
        return proteinProteinMatrix(newMatrix, newRowProteins, self.columnProteins, self.intlist, sigMatrix=newSigMatrix)
        


    def findRowMatches(self, protein):
        #given protein object, search through coefficient matrix and return row with match
        return  [i for i in range(len(self.rowProteins)) if self.rowProteins[i].checkMatch(protein)]
        
    
    def findColumnMatches(self, protein):
        #given protein object, search through coefficient matrix and return column with match
        return [i for i in range(len(self.columnProteins)) if self.columnProteins[i].checkMatch(protein)]
    
    def detectInteractors(self,bait, sort=None):            #given a bait, returns the row numbers where it's known interacting proteins are present and a list of those proteins
        plotRows = []
        proteinList = []
        for i in range(len(self.columnProteins)):
            if(self.columnProteins[i].gene_symbol == bait):
                baitProtein = self.columnProteins[i]
                col_index = i+1
        for i in range(len(self.intlist.proteinPairList)):
                if(baitProtein.checkMatch(self.intlist.proteinPairList[i][0])):
                    if self.findRowMatches(self.intlist.proteinPairList[i][1]):
                        matched_rows = [t for t in self.findRowMatches(self.intlist.proteinPairList[i][1])]
                        #matched_proteins = [(baitProtein,self.rowProteins[t[1]]) for t in self.findRowMatches(self.intlist.proteinPairList[i][1])]
                        matched_proteins = self.intlist.proteinPairList[i]
                        plotRows.extend(matched_rows)
                        proteinList.append(matched_proteins)
                elif(baitProtein.checkMatch(self.intlist.proteinPairList[i][1])):
                    if self.findRowMatches(self.intlist.proteinPairList[i][0]):
                        matched_rows = [t for t in self.findRowMatches(self.intlist.proteinPairList[i][0])]
                        matched_proteins = self.intlist.proteinPairList[i]
                        #matched_proteins = [(baitProtein,self.rowProteins[t[1]]) for t in self.findRowMatches(self.intlist.proteinPairList[i][0])]
                        plotRows.extend(matched_rows)
                        proteinList.append(matched_proteins)
                else:
                    continue
        if(sort):
            plotRows = sorted(plotRows, key=lambda i: self.matrix[i, col_index], reverse=True)
        return plotRows, proteinList

    def extractInteractions(self, bait, top):
        plotRows = []
        for i in range(len(self.columnProteins)):
            if(self.columnProteins[i].gene_symbol == bait):
                baitProtein = self.columnProteins[i]
                col_index = i+1
        plotRows = [i for i in range(self.matrix.shape[0])]
        plotRows = sorted(plotRows, key= lambda i:self.matrix[i, col_index], reverse=True)
        return plotRows[:top]
        #extract top "top" interactions of given bait protein


    def plotHeatmap(self,plotRows=None, plotDiagonal=False,proteinProteinAnnotations=None):
        #if plotDiagonal is True for each row protein, check if there is a same protein in the column
        #annotations will contain the proteins to be highlighted
        diagonalPairs = []
        ncols = self.matrix.shape[1]
        baits = [protein.gene_symbol for protein in self.columnProteins]
        numberIntercepts = ncols - len(baits)
        if(plotDiagonal):
            for i in range(len(plotRows)):
                colMatches = self.findColumnMatches(self.rowProteins[plotRows[i]])
                if len(colMatches)>0:
                    bait = [t for t in colMatches]
                    diagonalPairs.extend([(i,col+numberIntercepts) for i in [i] for col in bait])


        output = self.matrix[plotRows,:] #np.stack(matrixPlot, axis=0)
        
        
        for i in range(numberIntercepts):
            baits.insert(i,'Intercept')
        preys = [self.rowProteins[i].gene_symbol for i in plotRows] 
        rowcollist = []
        for i in range(len(proteinProteinAnnotations)):
            globalRows = [r for r in self.findRowMatches(proteinProteinAnnotations[i][0])]
            localRows = [k for k, elem in enumerate(plotRows) if elem in globalRows]
            localCols = [c+numberIntercepts for c in self.findColumnMatches(proteinProteinAnnotations[i][1])]
            if(localRows and localCols):
                to_annotate = [(row,col) for row in localRows for col in localCols]
                rowcollist.extend(to_annotate)
            globalRows = [r for r in self.findRowMatches(proteinProteinAnnotations[i][1])]
            localRows = [k for k, elem in enumerate(plotRows) if elem in globalRows]
            localCols = [c+numberIntercepts for c in self.findColumnMatches(proteinProteinAnnotations[i][0])]
            if(localRows and localCols):
                to_annotate = [(row,col) for row in localRows for col in localCols]
                rowcollist.extend(to_annotate)
        #rowcollist.extend(diagonalPairs)
        rowColList = set(rowcollist)
        
      
        
        fig, ax = plt.subplots()
        beta_heatmap = ax.imshow(output, cmap='Blues') 
        plt.colorbar(beta_heatmap)
        ax.set_xticks(np.arange(output.shape[1]))
        ax.set_xticklabels(baits, rotation=45, ha='left', fontsize=8)

        ax.set_yticks(np.arange(output.shape[0]))
        ax.set_yticklabels(preys, fontsize=10)
        # Move x-axis labels to the top
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')


        for (row, col) in rowColList:
            # Find the corresponding row in the heatmap
            # Draw rectangle around the corresponding cell
            rect = Rectangle((col -0.5, row - 0.5), 1, 1, fill=False, edgecolor='red', lw=2)
            ax.add_patch(rect)
        
        for (row,col) in diagonalPairs:
            rect = Rectangle((col -0.5, row - 0.5), 1, 1, fill=False, edgecolor='yellow', lw=2)
            ax.add_patch(rect)
        return rowColList, diagonalPairs, output
    
    def writeToFile(self,outFile):
        with open(outFile, "w+") as f:
            #Writes the header
            f.write("\t\t%s\n"%("\t".join([ p.gene_symbol for p in self.columnProteins])))
            f.write("\t\t%s\n"%("\t".join([ ",".join(p.uniprot_set) for p in self.columnProteins])))
    
            #Writes the content
            for i in range(len(self.matrix)):
                rowSymbol = self.rowProteins[i].gene_symbol
                rowAnnotation = ",".join(self.rowProteins[i].uniprot_set)
                if self.sigMatrix is None: #Only write the values
                    matrixRow =  self.matrix[i]
                    rowData = "\t".join([ "%.3f"%mi for mi in matrixRow])
                else:
                    matrixRow =  self.matrix[i]
                    sigRow    =  self.sigMatrix[i]
                    rowData = "\t".join([ "%.3f,%.2f"%(matrixRow[j], sigRow[j]) for j in range(len(matrixRow))])
                f.write("%s\t%s\t%s\n"%(rowSymbol, rowAnnotation, rowData))
    
    def findMissingBaits(self):
        """Finds the bait proteins which were not detected in the experiment"""
        """i.e. finds the proteins present in the columns of the mixing matrix, but not found in the rows of the signal matrix"""

        missingBaits = []
        interceptCount = 0
        for i in range(len(self.columnProteins)):
            if (self.columnProteins[i].findProteinMatches(self.rowProteins)):
                continue
            elif("BG:" in self.columnProteins[i].gene_symbol):
                interceptCount+=1
                continue
            else:
                missingBaits.append((i-interceptCount,self.columnProteins[i]))
        return missingBaits
    
    def includeMissingBaits(self):
        """insert a dummy row of NaNs into the beta matrix, which stand in for missing bait entries"""
        """also insert missing bait proteins into self.rowProteins"""
        """this only works if we are only solving for the set of bait proteins"""
        missingBaits = self.findMissingBaits()
        newMatrix = np.copy(self.matrix)
        newSigMatrix = np.copy(self.sigMatrix)
        newRowProteins = self.rowProteins.copy()
        for i in range(len(missingBaits)):
            newMatrix = np.insert(newMatrix, missingBaits[i][0], np.nan*np.ones(len(self.columnProteins)), axis=0)
            newSigMatrix = np.insert(newSigMatrix, missingBaits[i][0], np.nan*np.ones(len(self.columnProteins)), axis=0)
            newRowProteins.insert(missingBaits[i][0], missingBaits[i][1])
        return proteinProteinMatrix(newMatrix, newRowProteins, self.columnProteins, self.intlist, sigMatrix=newSigMatrix)
        #np.nan*np.ones(len(self.columnProteins))
        

    
    
    def plotCombinedMixSignalBetaPlot(self, trainingData=None, proteinIndices=None, mixingIndices=None, imageScale=0.5, xPad=0, yPad=0, highlightBaits=True, ppiAnnotations=[], specialOrientation=False, sigThreshold=10, ax=None):
    

        
        #Gather relevant information information
        betaMatrix=self.matrix
        preyProteins=self.rowProteins
        abProteins=self.columnProteins
        
        #Plots all rows if none are specified
        if proteinIndices is None:
            proteinIndices = range(len(self.rowProteins))
        nProt = len(proteinIndices)
    
        nCols=1
        #Plots all proteins if no indices are specificed
        if proteinIndices is None:
            proteinIndices = range(len(self.matrix))

        #Determines what to plot and layout
        if trainingData is None:
            #Not showing signal matrix (or mixing matrix)
            (plotSignal, plotMixing)=(False, False)
    
            f, aBeta = plt.subplots(1, 1, gridspec_kw={'width_ratios': [1], 'height_ratios': [1]})
            (aMix, aSignal) = (None, None)
            (f_height, f_width) = (nProt+yPad, self.nCols+xPad)

        else:
            #Processes the signal matrix
            signalMatrix = trainingData.signalMatrix
            nPools = signalMatrix.shape[1]
            
            # = np.insert(signalMatrix, 1, np.nan*np.ones(nPools), axis=0)
            if mixingIndices is None:
                #Showing signal matrix, but not mixnig matrix
                (plotSignal, plotMixing) = (True, False)
                f, (aSignal, aBeta) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [nPools, self.nCols], 'height_ratios': [1]})
                aMix = None
                (f_height, f_width) = (nProt+yPad, self.nCols+nPools+xPad)
            else:
                #Showing both signal matrix and mixing matrix
                nMixPlot=len(mixingIndices)
                mixingMatrix=trainingData.mixingMatrix

                (plotSignal, plotMixing) = (True, True)
                (nCols, nRows) = (3,1)
                (widthRatios, heightRatios) = ([nPools, self.nCols], [nProt,nMixPlot])
                (f_height, f_width) = (nProt+nMixPlot+yPad, self.nCols+nPools+xPad)
                
                
                #if we want special plot with three matrices side-by-side
                if(specialOrientation):
                    f, ((aMix, aSignal, aBeta)) = plt.subplots(1,3, gridspec_kw={'width_ratios': [nPools, nPools, self.nCols], 'height_ratios': [nMixPlot]}) #, nProt]})
                else:
                    f, ((aMix, aNone),(aSignal, aBeta)) = plt.subplots(2, 2, gridspec_kw={'width_ratios': [nPools, self.nCols], 'height_ratios': [nMixPlot, nProt]})
                    f.delaxes(aNone)
                
        plt.subplots_adjust(wspace=0, hspace=0)
        if(specialOrientation):
         plt.subplots_adjust(wspace=0.01, hspace=0)   
        


        #Formats AB names row names
        colNames = []
        for p in abProteins:
            temp = p.gene_symbol.split(';')
            if len(temp)>1:
                colNames.append("%s;..."%temp[0])
            else:
                colNames.append(temp[0])

        #Formats prey names row names
        preyNames = []
        for ip in proteinIndices:
            temp = preyProteins[ip].gene_symbol.split(';')
            if len(temp)>1:
                preyNames.append("%s;..."%temp[0])
            else:
                preyNames.append(temp[0])
                
        
        #Plots the beta coefficients
        #if plotMixing:
        #    self.plotBetaHeatmapWithMissingBaits(aBeta, missingProts, proteinIndices=proteinIndices, ppiAnnotations=ppiAnnotations, preyNames=None, colNames=colNames, highlightBaits=highlightBaits, sigThreshold=sigThreshold)

        if plotSignal:
            self.plotBetaHeatmap(aBeta, proteinIndices=proteinIndices, ppiAnnotations=ppiAnnotations, preyNames=None, colNames=colNames, highlightBaits=highlightBaits, sigThreshold=sigThreshold)
        
        else:
            self.plotBetaHeatmap(aBeta, proteinIndices=proteinIndices, ppiAnnotations=ppiAnnotations, preyNames=preyNames, colNames=colNames, highlightBaits=highlightBaits, sigThreshold=sigThreshold)

        #Plots the mixing matrix
        if plotMixing:
            self.plotMixingMatrix(aMix, specialOrientation, trainingData, mixingIndices=mixingIndices, nPools=nPools, nMixPlot=nMixPlot)

        #Plots the signal matrix
        if plotSignal:
            self.plotSignalMatrix(aSignal, specialOrientation, signalMatrix, proteinIndices = proteinIndices, preyNames=preyNames, nPools=nPools, nProt= nProt)
            
    
        f.set_size_inches(f_width*imageScale, f_height*imageScale)
        #f.savefig("exp.pdf")
        plt.show()
        
       
        
        #plt.close()
        #f.tight_layout()

    def plotMixingMatrix(self, ax, specialOrientation, trainingData, mixingIndices=None, nPools=None, nMixPlot=None ):
        aMix = ax
        mixingMatrix=trainingData.mixingMatrix
        aMix.imshow(mixingMatrix.T[mixingIndices,:]/np.max(mixingMatrix), cmap=self.cmapMix, vmin=0, vmax=1.0) 
        aMix.get_xaxis().set_visible(False)
        if(specialOrientation):
            aMix.get_xaxis().set_visible(True)  
            aMix.set_xticks(np.arange(nPools))
            aMix.set_xticklabels(range(1,nPools+1), rotation=0, ha='left', fontsize=12)
        aMix.set_yticks(np.arange(nMixPlot))
        #Note abProtesin includes the intercept
        aMix.set_yticklabels([ trainingData.abProteins[ip].gene_symbol for ip in mixingIndices], rotation=0, ha='right', fontsize=12)
    
    def plotSignalMatrix(self, ax, specialOrientation, signalMatrix, proteinIndices=None, preyNames = None, nPools = None, nProt = None):
        """Plots the signal matrix for the protein indices"""
        aSignal = ax
        aSignal.imshow(signalMatrix[proteinIndices,:], cmap=self.cmapSig, vmin=0, vmax=1.25) 
        aSignal.set_yticks(np.arange(nProt))
        aSignal.set_yticklabels(preyNames, rotation=0, ha='right', fontsize=12)
        if(specialOrientation):
            aSignal.set_yticks([])
            aSignal.set_yticklabels([])
        aSignal.set_xticks(np.arange(nPools))
        aSignal.set_xticklabels(range(1,nPools+1), rotation=0, ha='left', fontsize=12)



    def plotBetaHeatmap(self, ax, proteinIndices=None, ppiAnnotations=[], preyNames=None, colNames=None, highlightBaits=True, sigThreshold=None):
        """Plots a heatmap of the beta coefficients for the specified proteins."""
        
        betaMatrix=self.matrix
        #Plots all rows if none are specified
        if proteinIndices is None:
            proteinIndices = range(len(self.rowProteins))

        nProt = len(proteinIndices)
        
        #Identifies 'diagonal' bait proteins.
        diagonalPairs = []
        baits = [protein.gene_symbol for protein in self.columnProteins]
        if(highlightBaits):
            for (rowInPlot, iProt) in enumerate(proteinIndices):
                baitIndices = self.findColumnMatches(self.rowProteins[iProt])
                if len(baitIndices)>0:
                    diagonalPairs.extend([(rowInPlot,col) for i in [rowInPlot] for col in baitIndices])

        
        #Searches for all matches with the with set of annotation proteins 
        rowcollist = []
        plotProteins = [ self.rowProteins[i] for i in proteinIndices]
        for currentAnnotation in ppiAnnotations:
            
            #Checks if p1=row and p2=col
            localRows = [i for i in range(len(plotProteins)) if plotProteins[i].checkMatch(currentAnnotation[0])]
            localCols = [c for c in self.findColumnMatches(currentAnnotation[1])]
            if(localRows and localCols):
                rowcollist.extend([(row,col) for row in localRows for col in localCols])

            #Checks if p2=row and p1=col
            localRows = [i for i in range(len(plotProteins)) if plotProteins[i].checkMatch(currentAnnotation[1])]
            localCols = [c for c in self.findColumnMatches(currentAnnotation[0])]
            if(localRows and localCols):
                rowcollist.extend([(row,col) for row in localRows for col in localCols])
        rowcolset = set(rowcollist)

        #Creates heatmap for betas
        beta_heatmap = ax.imshow(betaMatrix[proteinIndices,:], cmap=self.cmapPPI, vmin=-1, vmax=1.0) 
        ax.set_xticks(np.arange(self.nCols))
        if colNames is None:
            ax.get_xaxis().set_visible(False)
        else:
            ax.set_xticklabels(colNames, rotation=90, ha='left', fontsize=12)
            
        if preyNames is None:
            ax.get_yaxis().set_visible(False)
        else:
        #if plotSignal is False:
            ax.set_yticks(np.arange(nProt))
            ax.set_yticklabels(preyNames, rotation=0, ha='right', fontsize=12)
            
        
        #Adds the significance markers
        if self.sigMatrix is not None and sigThreshold is not None:
            for (iPlotRow, iSigRow) in enumerate(proteinIndices):
                for iPlotCol in range(self.nCols):
                    if self.sigMatrix[iSigRow,iPlotCol] > sigThreshold:
                        ax.text(x=iPlotCol, y=iPlotRow, s="\u2020", ha="center", va="center")
            
        #Add the annotations.
        for (row, col) in rowcolset:
            # Find the corresponding row in the heatmap
            # Draw rectangle around the corresponding cell
            rect = Rectangle((col -0.5, row - 0.5), 1, 1, fill=False, edgecolor='black', lw=2, linestyle = 'dashed')
            ax.add_patch(rect)
        
        #Adds the 'diagonal' annotations    
        for (row,col) in diagonalPairs:
            rect = Rectangle((col -0.5, row - 0.5), 1, 1, fill=False, edgecolor='black', lw=2)
            ax.add_patch(rect)
        
        
    def plotBetaVersusF(self, ppiAnnotations=[], yAxis="linear"):
        n_rows, n_cols=5,5

        colProts = self.columnProteins
        n_cols = 6
        n_rows = int(np.ceil(1.*len(colProts)/n_cols))
        scaleFactor=3

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*scaleFactor, n_rows*scaleFactor), sharex=True, sharey=True)

        for i in range(n_cols*n_rows):
            row, col = divmod(i, n_cols)
            ax = axes[row, col]
            if i<len(self.columnProteins):
        
                #Gets all protein corresponding tho all PPIs with the column
                matchingPPIs = [ an[1] for an in ppiAnnotations if colProts[i].checkMatch(an[0])] + [ an[0] for an in ppiAnnotations if colProts[i].checkMatch(an[1])]
                #Checks if each row has an annotation
                cols = []
                for p in self.rowProteins:
                    if p.checkMatch(colProts[i]):
                        cols.append("red")
                    elif len(p.findProteinMatches(matchingPPIs))>0:
                        cols.append("black")
                    else:
                        cols.append("grey")
                        
                        
                
                #Builds scater plots
                if yAxis=="linear" or yAxis=="log":
                    ax.scatter(self.matrix[:,i], self.sigMatrix[:,i], c=cols)
                elif yAxis=="logFplus1":
                    ax.scatter(self.matrix[:,i], self.sigMatrix[:,i]+1, c=cols)
                else:
                    raise ValueError("Argument yAxis must take value 'linear', 'log', or 'logFplus1'")
                
                ax.set_title(colProts[i].gene_symbol, fontsize=11)
                if yAxis=="log" or yAxis=="logFplus1":
                    ax.set_yscale('log')
                if row==n_rows-1 and col==0:
                    ax.set_xlabel('Coefficient', fontsize=8)
                    if yAxis=="logFplus1":
                        ax.set_ylabel('F+1', fontsize=8)
                    else:
                        ax.set_ylabel('F', fontsize=8)
            else:
                fig.delaxes(ax)
                    
        plt.show()

      



        