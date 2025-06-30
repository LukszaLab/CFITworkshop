#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 08:06:46 2021

@author: zacharysethna
"""

from cfit.fitness.CloneFitness import CloneFitness


class VariableImmuneCloneFitness(CloneFitness):
    '''
    Immune selection component, accounts for the negative selection on neoantigens in the clone.

    Attributes:
        quality_function: function
            defines neoantigen quality function, arguments are Neoantigen and Sample objects.

        aggrfun: function
            how the neaontigen qualities are aggregated over multiple neoantigens in the clone.

    '''

    def __init__(self, quality_function=lambda neo, sampleTree, sample: neo.quality,
                 aggrfun=max, **kwargs):
        '''
        Constructor

        :param quality_function:
            function called on Neoantigen, SampleTree and Sample objects

        :param aggrfun: function
            defaults to max, a function to aggregate neoantigen qualities

        :param kwargs: dict
            named parameters to the aggregating function
        '''

        super(VariableImmuneCloneFitness, self).__init__("immune")
        self.quality_function = quality_function
        self.aggrfun = aggrfun
        self.aggrparams = kwargs

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

        qualities = [self.quality_function(neo, sampleTree, sample) for neo in node.neoantigens]
        # qualities = [neo.quality for neo in node.neoantigens]
        if len(qualities) == 0:
            fitness = 0.0
        else:
            fitness = -self.aggrfun(qualities, **self.aggrparams)

        node.fitness_components[self.name] = fitness
