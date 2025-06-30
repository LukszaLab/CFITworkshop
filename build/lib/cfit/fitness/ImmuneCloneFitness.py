from cfit.fitness.CloneFitness import CloneFitness


class ImmuneCloneFitness(CloneFitness):
    '''
    Immune selection component, accounts for the negative selection on neoantigens in the clone.

    Attributes:
        aggrfun: function
            how the neaontigen qualities are aggregated over multiple neoantigens in the clone.

    '''

    def __init__(self, aggrfun=max, **kwargs):
        '''
        Constructor

        :param aggrfun: function
            defaults to max, a function to aggregate neoantigen qualities

        :param kwargs: dict
            named parameters to the aggregating function

        '''
        super(ImmuneCloneFitness, self).__init__("immune")
        self.aggrfun = aggrfun
        self.aggrparams = kwargs

    def compute_node_fitness(self, node, sampleTree, sample):
        '''

        Computes the values of the fitness component for the node and sets fitness_components[component_name] attribute
        of the node. The simple model, with no neaontigen competition

        :param node: cfit.tree.node.Node
            the clone for which the fitness component is evaluated.

        :param sampleTree: cfit.tree.SampleTree
            the clonal structure to which the clone belongs in the sample

        :param sample: cfit.patient.Sample
            the sample to which the clone belongs

        '''
        #        qualities = [neo.quality*sample.mutations[neo.mid].expression for neo in node.neoantigens]
        qualities = [sampleTree.neoantigenQualities[neo.id] for neo in node.neoantigens]
        if len(qualities) == 0:
            fitness = 0.0
        else:
            fitness = -self.aggrfun(qualities, **self.aggrparams)
        node.fitness_components[self.name] = fitness

    def compute_node_fitness_dynamic(self, node, sampleTree, sample, **kwargs):
        '''

        Computes the values of the fitness component for the node and sets fitness_components[component_name] attribute
        of the node.

        :param node: cfit.tree.node.Node
            the clone for which the fitness component is evaluated.

        :param sampleTree: cfit.tree.SampleTree
            the clonal structure to which the clone belongs in the sample

        :param sample: cfit.patient.Sample
            the sample to which the clone belongs

        :param kwargs: dict
            parameters to the aggreating function

        '''
        #        qualities = [neo.quality*sample.mutations[neo.mid].expression for neo in node.neoantigens]
        qualities = [sampleTree.neoantigenQualities[neo.id] for neo in node.neoantigens]
        if len(qualities) == 0:
            fitness = 0.0
        else:
            fitness = -self.aggrfun(qualities, **kwargs)

        node.fitness_components[self.name] = fitness
