'''
Created on Feb 22, 2016

@author: mluksza
'''
from cfit.CoreObject import CoreObject


class GeneralNode(CoreObject):
    '''
    Class implementing the basic node in a tree functionalities.

    Attributes:
        __id: int
            node identifier

        name: str:
            name of the node "internal"

        parent: cfit.tree.node.GeneralNode
            the link to the parent on the tree

        children: list
            the list of children nodes

    '''

    def __init__(self, nid):
        '''
        Constructor method

        :param nid: int
            node identifier
        '''
        self.__id = nid
        self.name = "internal"
        self.parent = self
        self.children = []

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, nid):
        self.__id = nid

    def copyAttributes(self, cnode):
        '''
        Initializes the node based on another node.

        :param cnode: cfit.tree.node.GeneralNode
        :return:
        '''
        self.__id = cnode.__id
        self.name = cnode.name
        self.parent = cnode.parent
        self.children = cnode.children

    def path_to_root(self):
        '''
        Return path of nodes to the root.
        :return: list
            list of GeneralNode objects
        '''
        if self.parent.id == self.id:
            return [self]
        rpath = self.parent.path_to_root()
        rpath.append(self)
        return rpath

    def get_subtree_nodes(self):
        '''
        Returns all the nodes in the subtree.
        :return: set
            set of nodes
        '''

        nodes = set()
        for cnode in self.children:
            lnodes = cnode.get_subtree_nodes()
            for lnode in lnodes:
                nodes.add(lnode)
        nodes.add(self)
        return nodes

    def add_child(self, c):
        '''
        Adds a child node.
        :param c: cfit.tree.node.GeneralNode
        :return:
        '''
        children = set(self.children)
        if self.__id != c.__id:
            children.add(c)
        children = list(children)
        children.sort(key=lambda node: node.__id)
        self.children = children
        c.set_parent(self)

    def set_parent(self, p):
        '''
        Sets the parent node

        :param p: GeneralNode
        :return:
        '''
        self.parent = p

    def tree_depth(self):
        if self.parent.id == self.id:  # root node
            depth = 0
        else:
            depth = self.parent.tree_depth() + 1
        return depth
