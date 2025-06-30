'''
Created on Nov 22, 2018

@author: mluksza
'''

import time
import sys

class Log(object):
    '''
    classdocs
    '''

    VERBOSE = 3
    TSTART = time.time()

    @staticmethod
    def logger(msg, level=0, warn=False):
        """
        Print logger message msg to stdout.

        Parameters
        -----------

         msg : str
            String to print on the screen

         level : int
            Log-level. Only the messages with a level higher than the
            current verbose level will be shown.

         warn : bool
            Warning flag. If True, the message will be displayed
            regardless of its logger-level.

        """


        if level < Log.VERBOSE or (warn and level <= Log.VERBOSE):
            msg = str(msg)
            dt = time.time() - Log.TSTART
            outstr = str('\n' if level < 2 else '')
            outstr += str(format(dt, '4.2f')) + '\t'
            outstr += level * '-'
            outstr = str(outstr)
            outstr += msg
            sys.stdout = sys.__stdout__
            print(outstr)
