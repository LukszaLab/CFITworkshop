'''
Created on Nov 28, 2018

@author: mluksza
'''
import json
import time
import sys

class CoreObject(object):
    '''
    classdocs
    '''

    VERBOSE = 3
    TSTART = time.time()

    def toJSON(self):
        '''

        :return: dict
        '''
        js = {}
        return js

    def writeJSON(self, of):
        '''

        :param of: file handle

        :return:
        '''
        js = self.toJSON()
        json.dump(js, of, indent=True)

    def __logger(self, msg, level=0, warn=False, prefix=""):
        """
        Print logger message msg to stdout.

        Parameters
        -----------

         :param msg : str
            String to print on the screen

         :param level : int
            Log-level. Only the messages with a level higher than the
            current verbose level will be shown.

         :param warn : bool
            Warning flag. If True, the message will be displayed
            regardless of its logger-level.

        """
        if level < self.VERBOSE or (warn and level <= self.VERBOSE):
            msg = str(msg)
            dt = time.time() - self.TSTART
            outstr = '\n' if level < 2 else ''
            outstr += format(dt, '4.2f') + '\t'
            outstr += level * '-'
            outstr += prefix + ": " + msg
            print(outstr)  # , file=sys.stdout)

    def logger(self, msg, level=0, warn=False):
        '''
        Prints the msg, preceeded by time and the name of the class

        :param msg: str
            message

        :param level: int
            verbosity level

        :param warn: bool

        '''

        sys.stdout = sys.__stdout__
        prefix = str(type(self)).strip()
        prefix = prefix.replace("<class '", "")
        prefix = prefix.replace("'>", "")
        prefix = prefix.split(".")[-1]
        self.__logger(msg, level=level, warn=warn, prefix=prefix)
