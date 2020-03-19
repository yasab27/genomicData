#!/usr/bin/env python
# coding: utf-8

# The purpose of this notebook is simply to create pathway and module correlation matrices for all all four organsisms studied. We will be accomplishing this by use of the `genomicAbundanceSample` class.

# In[2]:


# Begin by importing all necessary dependencies for this class. We will require numerical computation ability as well 
# as the ability to parse HTML documents.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import scipy.spatial as spat

class genomicAbundanceSample():
    """
        This class models microbiome genomic abundance for a particular subpopulation for various KEGG ortholog named
        gene across multiple populations. For example, a sample of all humans, all Americans, all macacs sampled
        maybe modeled. 
        
        Params:
            -dataString: The string which refers to the .xlxs document within the directory storing the
            abundance data within the text.
            -htmlParserString: The string which refers to the .html document with the directory that provides
            the list for categorizing the individual genes into larger functional modules or pathways. 
    """
    
    def __init__(self, name, dataString, htmlParserString, bounds=None):
        """
            Initialize a new genomicAbudance sample with the above parameters. Additionally, generate automatically
            dataframe representations of both the pure KO genetic abundance data, the module/pathway abundance data,
            and a similarity matrix (using correlation distance). The bounds measure gives the columns of the
            excel to consider when constructing the matrix. 
            
            Params:
                -bounds: a two element array which is [lowerBound, upperBound]
        """
        self.name = name
        self.dataString = dataString
        self.htmlParserString = htmlParserString
        
        self.KOAbundances = self.generateKODataFrame(dataString)
        self.modules = self.generateModuleDictionary(htmlParserString)
        self.moduleAbundances = self.generateModuleDataFrame(bounds)
        self.correlationVector = self.generateModuleCorrelationVector()
        self.moduleSimilarity = self.generateModuleSimilarityMatrix()
        
    def generateKODataFrame(self,dataString):
        """
            Given a certain string reference to a .xlxs file in the directory, generate a pandas dataframe
            showing the abundance data for all the KO genes within the database. 
        """
        # Read in the file from the database, using the first row as the header. 
        df = pd.read_excel(dataString, header = 0)
        # Rename the first column to the KEGG column and then set this as the index. 
        df.rename(columns={"Unnamed: 0":"KEGG"}, inplace=True)
        df = df.set_index("KEGG")
        
        return df
    
    def generateModuleDictionary(self, htmlParserString):
        """
            Given a reference to an html document from the database, generate a dictionary showing how genes break
            down in terms of modules. 
        """
        
        # Parse the HTML document containing our downloaded data relating to the mass genome search. 
        with open(htmlParserString, "r") as file:
            # Read in the file
            contents = file.read()
            # Construct a BS object to allow for easy parsing of the document. 
            soup = BeautifulSoup(contents, 'html')
        
        # Now we will grab the section of the HTML dedicated purely to the module classification. We first query for 
        # the div with class box2 which will contain all the information about this unique module return,
        # and then we subselect the sole list element within it which contains information about the genomic list. 
        moduleSection =  soup.find("div", {"class": "box2"}).find_all("li")

        # Store all module data
        modules = []

        # Now parse each list element and create a hashmap where keys correspond to the module name and the elements  
        # contain the actual gene id
        for module in moduleSection:

            # Create an object containing all elements for the dataset. 
            moduleObject = {}

            # First grab the ID of the module:
            moduleID = list(module.children)[0].getText()

            # Next grab the actual semantic name of the module
            moduleName = list(module.children)[1]

            # Create a new module object containing these elements 
            moduleObject["id"] = moduleID
            moduleObject["name"] = moduleName.replace("\n","").replace("\xa0","")[:-1] 
            # stray HTML tag and second indexing eliminates the last character (.

            # Now we add all of the genes in these modules to its objects
            geneList = module.find("div",{"class":"object"}).find_all("a")
            moduleObject["genes"] = []

            for gene in geneList:

                geneObject = {}
                # Acquire the gene ID and gene semantic identifier
                geneObject["geneID"] = gene.string[3:]
                geneObject["geneName"] = gene.nextSibling.replace("\n","").replace("\xa0","")
                moduleObject["genes"].append(geneObject)

            # Append the completed module object to the list of all modules
            modules.append(moduleObject)
        
        return modules
    
    def generateModuleDataFrame(self, bounds):
        """
            Generate the dataframe of module abundances by summing over the abundances for individual genes 
            within those modules as given by the module dictionary and KO gene abundances. 
        """
        # Start by creating a new pandas dataframe
        dfModules = pd.DataFrame()

        # Traverse through every found module and sum up all of the values of abundance for each of their genes
        # across every gene. This gives a rough measure of module abundance. 
        for module in self.modules:

            # Get a list of all the IDs of the genes within this module by traversing the list of gene object
            # within the module and then appending them to a larger list. 
            geneIDs = []
            for gene in module["genes"]:
                geneIDs.append(gene["geneID"])

            # Sum across all genes within module to generate one vector containing all the information about the
            # module abundance. Each entry in the vector corresponds to an individual being measured (i.e. abundance
            # of a specific module within one person's microbiome).

            # We accomplish this by first locating all genes of interset within the dataframe, summing them, combining
            # the summed vector into a pd series and then appending it to the larger dataframe. 
            summedAbundance = pd.Series(np.sum(self.KOAbundances.loc[geneIDs,:] ,axis = 0), name = module["id"])
            dfModules = dfModules.append(summedAbundance)

        # Drop all rows with 0 abundance
        dfModules = dfModules[dfModules.values.sum(axis=1) != 0]
        
        # Check to see if there are bounds
        if bounds is None:
            return dfModules
        else:
            return dfModules.iloc[:,bounds[0]:bounds[1]]
    
    def generateModuleSimilarityMatrix(self):
        """
            Given the module abundance matrix, construct a similarity matrix between points using 1-r 
            as the similarity matrix. 
        """
        similarity = spat.distance.squareform(spat.distance.pdist(self.moduleAbundances, "correlation"))

        simDF = pd.DataFrame(similarity, columns= self.moduleAbundances.index, index=self.moduleAbundances.index)
        return simDF
    
    def generateModuleCorrelationVector(self):
        """
            Given the module abundance matrix, construct a vector of the correlation corresponding to the upper
            triangle of the correlation matrix. The value at position (i,j) is the correlation between module i
            and module j.
        """
        correlationVector = 1-spat.distance.pdist(self.moduleAbundances, "correlation")
        return correlationVector


# In[7]:


# Start by createing all module maps
humanAbundances = genomicAbundanceSample("Human", "humans/raw/humanKOAbundances.xlsx","humans/raw/humanModules.html")


# In[11]:


correlationHuman = 1- humanAbundances.moduleSimilarity
correlationHuman.to_csv("humans/humanModuleCorrelation.csv")


# In[12]:


macacAbundances = genomicAbundanceSample("Macac", "macac/raw/macacKOAbundances.xlsx","macac/raw/macacModules.html")
correlationMacac = 1- macacAbundances.moduleSimilarity
correlationMacac.to_csv("macac/macacModuleCorrelation.csv")


# In[13]:


miceAbundances = genomicAbundanceSample("Mice", "mice/raw/mouseKOAbundances.xlsx","mice/raw/mouseModules.html")
correlationMouse = 1- miceAbundances.moduleSimilarity
correlationMouse.to_csv("mice/mouseModuleCorrelation.csv")


# In[14]:


pigAbundances = genomicAbundanceSample("Pig", "pig/raw/pigKOAbundances.xlsx","pig/raw/pigModules.html")
correlationPig = 1- pigAbundances.moduleSimilarity
correlationPig.to_csv("pig/pigModuleCorrelation.csv")


# In[16]:


# Now for pathways
humanAbundances = genomicAbundanceSample("Human", "humans/raw/humanKOAbundances.xlsx","humans/raw/humanPathways.html")
correlationHuman = 1- humanAbundances.moduleSimilarity
correlationHuman.to_csv("humans/humanPathwayCorrelation.csv")


# In[18]:


macacAbundances = genomicAbundanceSample("Macac", "macac/raw/macacKOAbundances.xlsx","macac/raw/macacPathways.html")
correlationMacac = 1- macacAbundances.moduleSimilarity
correlationMacac.to_csv("macac/macacPathwayCorrelation.csv")


# In[19]:


miceAbundances = genomicAbundanceSample("Mice", "mice/raw/mouseKOAbundances.xlsx","mice/raw/mousePathways.html")
correlationMouse = 1- miceAbundances.moduleSimilarity
correlationMouse.to_csv("mice/mousePathwayCorrelation.csv")


# In[21]:


pigAbundances = genomicAbundanceSample("Pig", "pig/raw/pigKOAbundances.xlsx","pig/raw/pigPathways.html")
correlationPig = 1- pigAbundances.moduleSimilarity
correlationPig.to_csv("pig/pigPathwayCorrelation.csv")


# Now we'll also construct correlation matrices for abundances of bacterial genus.

# In[ ]:




