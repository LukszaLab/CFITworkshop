'''
Created on Mar 30, 2015

@author: mluksza
'''

from collections import defaultdict

from cfit.CoreObject import CoreObject


class Epitope(CoreObject):
    '''
    IEDB epitope class
    '''

    def __init__(self, row, eid):
        '''
        Constructor

        :param row: str
            str(seq.seq) + "\t" + seq.description

        :param eid: int
            epitope identifier, extracted from the fasta file header,
            >eid|***
        '''

        self.__id = eid
        self.shift = 0
        [seq, desc] = row.split("\t")
        seq = seq.split("+SCM(")[0].replace(".", "").replace('"', "")
        self.seq = seq

        self.startpos = 0
        self.endpos = 0
        self.antigen = 'Unknown'
        self.antigenId = 'Unknown'
        self.organism = 'Unknown'
        self.organismId = -1
        self.species = 'Unknown'

        dtab = desc.split("|")
        [iedbid, antigen, antigenId, organism, organismId] = [int(dtab[0]), "Unknown", -1, "Unknown", -1]
        shift = 0
        org2id = defaultdict(lambda: 0)
        try:
            ddtab = desc.split("|")
            if len(ddtab) == 6:
                [iedbid, shift, antigen, antigenId, organism, organismId] = ddtab
            elif len(ddtab) == 5:
                [iedbid, antigen, antigenId, organism, organismId] = ddtab
            elif len(ddtab) == 4:
                [iedbid, antigen, organism, organismId] = ddtab
            if len(organismId) == 0:
                organism = "Unknown"
                organismId = -1
            org2id[organism] = organismId
        except:
            self.logger("exception", 0)
            pass
        self.iedbid = iedbid
        self.shift = shift
        self.antigen = antigen
        self.antigenId = antigenId
        self.organism = organism
        try:
            self.organismId = int(organismId)
        except ValueError:
            self.organismId = org2id[organism]

        self.name = self.get_eheader()

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, id):
        self.__id = id

    def get_eheader(self):
        # eheader=self.antigen+"|"+self.antigenId+"|"+self.organism+"|"+self.organismId
        eheader = str(self.organismId)
        eheader = eheader.replace(" ", "_")
        return eheader
