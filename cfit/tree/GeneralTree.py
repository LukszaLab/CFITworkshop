'''
Created on Feb 22, 2016

@author: mluksza
'''
from cfit.CoreObject import CoreObject


class GeneralTree(CoreObject):
    '''
    Class for representing the basic tree functionality.

    Attributes:
        root: cfit.tree.Node
            The root node of the tree

        nodes: dict
            dictionary mapping node identifiers to Node objects, storing all nodes on the tree.

    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.root = -1
        self.nodes = {}

    def connect_nodes(self, nid, pid):
        '''
        Connect node nid to parent node pid

        :param nid: int
            Identifier of the node

        :param pid: int
            Identifier of the parent node
        '''
        self.nodes[pid].add_child(self.nodes[nid])

    def get_nodes(self):
        '''
        Returns all nodes in the tree.

        :return: list
            list of nodes on the tree
        '''

        return self.root.get_subtree_nodes()

    def set_root(self, rnode):
        '''
        Sets the given node to be the root of the tree

        :param rnode: cfit.tree.node.Node

        '''
        self.nodes[rnode.id] = rnode
        self.root = rnode
        self.root.parent = self.root

    def remove_node(self, nid):
        '''
        Removes a node of a given identifier from the tree. Reconnects children to grandparents.

        :param nid: int
        '''

        if nid != 0 and nid in self.nodes:
            rnode = self.nodes[nid]
            for cnode in rnode.children:
                rnode.parent.add_child(cnode)
            del self.nodes[nid]

    def get_pre_order(self):
        '''

        :return: list
        '''
        order = []
        queue = [self.root]
        while len(queue) > 0:
            top = queue[0]
            children = top.children
            queue = queue[1:]
            order.append(top)
            for child in children:
                queue.append(child)
        return order

    def get_post_order(self):
        '''

        :return: list
            list of GeneralNode objects
        '''

        order = []
        stack = [self.root]

        while len(stack) > 0:
            top = stack[-1]
            stack = stack[:-1]
            order.append(top)
            for cnode in top.children:
                stack.append(cnode)
        order.reverse()
        return order
