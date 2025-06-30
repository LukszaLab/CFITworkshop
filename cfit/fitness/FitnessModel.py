import numpy as np

from cfit.CoreObject import CoreObject


class FitnessModel(CoreObject):
    '''
    Class implementing fitness model with independent, additive components

    Attributes:
        components: dict: str -> cfit.fitness.CloneFitness
            maps component name to CloneFitness objects. eg. "immune" -> ImmuneCloneFitness
        weights: dict: str -> float

    '''

    def __init__(self):
        self.components = {}
        self.weights = {}

    def reset_fitness_model_components(self):
        self.components = {}
        self.weights = {}

    def add_component(self, fitnessComponent, name, weight):
        '''

        Adds a fitness model component

        :param fitnessComponent: cfit.fitness.CloneFitness
            the component object

        :param name: str
            name of the component

        :param weight: float
            weight of the component in the fitness function

        '''
        self.components[name] = fitnessComponent
        fitnessComponent.name = name

        self.weights[name] = weight

    def compute_node_fitness(self, node, sampleTree, sample, components=None, recompute_components=True):
        '''
        Evaluates fitness for clone 'node', sets node.fitness attribute

        :param node: cfit.tree.node.Node
            the clone for which the fitness is evaluated

        :param sampleTree: cfit.tree.SampleTree
            the clonal structure of the tree in the sample

        :param sample: cfit.patient.Sample
            the sample to which the clone belongs

        :param components: list
            list of str, component names to be included

        :param recompute_components: bool
            whether to recompute individual components
        '''

        if components is None:
            components = list(self.components.keys())

        if recompute_components:  # otherwise only the weights have changed
            for cname in components:
                self.components[cname].compute_node_fitness(node, sampleTree, sample)
        fitness = sum([node.fitness_components[name] * self.weights[name] for name in components])
        node.fitness = fitness


    def normalize_fitness(self, nodes, components=None):
        '''

        :param nodes: list/set
            collection of tuples of objects (cfit.tree.node.Node, cfit.tree.SampleTree, cfit.patient.Sample)

        :param components: list
            list of str, component names to be included
        '''

        self.logger("normalizing  fitness")
        if components is None:
            components = list(self.components.keys())
        for cname in components:
            for (node, sampleTree, sample) in nodes:
                self.components[cname].compute_node_fitness(node, sampleTree, sample)
            vals = [node.fitness_components[cname] for (node, sampleTree, sample) in nodes]
            std = np.std(vals)
            self.logger(cname + " std = " + str(std))
            self.logger("number of nodes: " + str(len(nodes)))
            for (node, sample) in nodes:
                node.fitness_components[cname] /= std
        fitness = sum([node.fitness_components[name] * self.weights[name] for name in components])
        node.fitness = fitness

    def get_parameters(self):
        p = {}
        for comp in self.components:
            p[comp + "_component_weight"] = self.weights[comp]
            # pc = self.components[comp].get_parameters()
            # for pname in pc:
            #    p[comp+ "_" + pname] = pc[pname]
        return p
