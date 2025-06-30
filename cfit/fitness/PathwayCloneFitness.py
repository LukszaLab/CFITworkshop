
import requests 
from cfit.fitness.CloneFitness import CloneFitness
#import random


class PathwayCloneFitness(CloneFitness):
    def __init__(self, DB_ID="KEGG:04151"): # Where DB_ID is a Keg / WP / GO name, such as Keg:04151 for Pi3k-AKT Signaling Pathway 
    
        self.DB_ID = DB_ID 
    
        
        r = requests.post(
                            url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
                            json={
                            'organism':'hsapiens',
                            'target':'UCSC',
                            'query':[DB_ID],
                            'sources' :["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP"], #only look into database terms 
                                    'user_threshold':1e-8, #reduce the significance threshold,
                                    'significance_threshold_method':'bonferroni', #use bonferroni correction instrad of the default 'g_SCS'.
                            'no_evidences':True, #skip lookup for evidence codes. Speeds up queries, if there is no interest in evidence codes.
                            'no_iea':True, #Ignore electonically annotated GO annotations
        
                              'domain_scope':'custom',#use the genes in the probe as the statistical background.
                       'background':'AFFY_HG_U133A'
                                },
                        headers={
                       'User-Agent':'FullPythonRequest'
                       }
        )
        
        
        
            
        self.genes = []
        
        for d in r.json()["result"]:
            # Loop through every individual dictionary in results list
            if d["name"] not in self.genes:
                self.genes.append(d["name"])
        #
        
                
    def compute_node_fitness(self, node, sampleTree, sample=None):

        # Return the number of items common between Pathway genes and (affected genes + node mutations)
        affected_node_genes = set() # affected genes set + node mutations list (?)
        
        for item in sample.affectedGenes:
        	affected_node_genes.add(item)
 
        for item in node.mutations:
        	affected_node_genes.add(item.gene)
        
        fitness = 0
        for gene in affected_node_genes:
        	if gene in self.genes:
        		fitness += 1
        
        #fitness = len([mut.gene for mut in node.mutations if mut.gene in self.genes]) 
        node.fitness_components[self.name] = fitness
            
        
        

