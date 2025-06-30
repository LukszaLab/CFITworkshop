from cfit.tree.mutation.Neoantigen import Neoantigen


class FrameShiftNeoantigen(Neoantigen):
    '''
    Class implements non-SNP neoantigens
    '''

    def __init__(self, params):
        '''
        Constructor
        '''
        Neoantigen.__init__(self, params)
