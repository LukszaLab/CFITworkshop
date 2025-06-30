from cfit.CoreObject import CoreObject


class CloneFitness(CoreObject):
    '''
    Abstract class for a fitness model component.

    Attributes:
        name: str
            name of the component
    '''

    def __init__(self, name):
        self.name = name

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
        pass
