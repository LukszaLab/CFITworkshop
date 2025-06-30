from cfit.CoreObject import CoreObject


class NeoantigenQuality(CoreObject):
    '''
    General class for neantigen quality computation


    Required methods:

    quality(neoantigen, args, **kwargs): float

    initialize_neoantigens(Analysis)
    '''

    def compute_quality(self, neo, args, **kwargs):
        '''

        Computes quality of the neoantigen, modifies neo.quality attribute

        :param neo: cfit.tree.mutation.Neoantigen

        :param args: list

        :param kwargs: dict

        :return: float
        '''

        pass

    def compute_neoantigen_sample_quality(self, neo, sample, args, **kwargs):
        '''

        Computes quality of the neoantigen, modifies sample.neoantigenQualities attribute

        :param neo: cfit.tree.mutation.Neoantigen

        :param sample: cfit.patient.Sample

        :param args: list

        :param kwargs: dict

        :return: float
        '''

        pass

    def initialize_neoantigens(self, anl, args, **kwargs):
        '''
        Initializes neoantigen quality related data

        :param anl: cfit.util.Analysis
        '''
        neos = anl.get_neoantigens()
        for neo in neos:
            self.set_quality_data(neo)

    def set_quality_data(self, neo):
        '''
        Set neoantigen specific attributes, needed for the quality evaluation

        :param neo: cfit.tree.mutation.Neoantigen

        '''
        pass

    def get_parameters(self):
        '''
        Dictionary of parameters specific to the model
        :return: dict
            str -> value
            name of the parameter -> value
        '''
        return {}

