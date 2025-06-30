from cfit.fitness.CloneFitness import CloneFitness


class DGCloneFitness(CloneFitness):
    '''
    Positive fitness component of driver genes

    '''

    def __init__(self, genes=["KRAS"]):
        '''
        Constructor

        :param gene: str
            name of the gene

        '''
        super(DGCloneFitness, self).__init__("DG")
        self.genes = genes

    def compute_node_fitness(self, node, sampleTree, sample):
        '''

        Computes the values of the fitness component for the node and sets fitness_components[component_name] attribute
        of the node.

        :param node: cfit.tree.node.Node
            the clone for which the fitness component is evaluated.

        :param sampleTree: cfit.tree.SampleTree
            the clonal structure to which the clone belongs in the sample

        :param sample: cfit.patient.Sample
            the sample to which the clone belongs
        '''
        fitness = len([mut.gene for mut in node.mutations if mut.gene in self.genes])
        node.fitness_components[self.name] = fitness
