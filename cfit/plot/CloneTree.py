
from cfit.plot.CloneNode import CloneNode

class CloneTree():

    """
        CloneTree object stores tree data structure with nodes CloneNodes.
    
    Attributes:

        root: class CloneNode
            clone 0 with no mutations 

        parents: list
            parent vector representation of tree 

        size: int
            number of clones in tree

        leaves: list
            list of CloneNode objects which are leaves

        height: int
            height of tree
        
        leaf2ord: dict (leaf CloneNode.id --> int)
           dictionary used to look up ordering of leaves within laderrized CloneTree
        
        nid2node: dict (CloneNode.id --> CloneNode object)
            used to look up CloneNode by its id

        nid2depth: dict (CloneNode.id --> CloneNode.depth)
            initialized from parent vector prior to setting depth attribute to individual CloneNodes

    """


    def __init__(self, parent_vector):

        """
        Params
        ----------
        parent_vector : list
            parent vector representation of tree topology
            entry at index i is parent id of clone i+1

        """


        self.root = CloneNode(id=0, depth=0)
        self.parents = parent_vector
        self.size = len(parent_vector) + 1
        
        self.nid2node = {i+1: CloneNode(id=i+1) for i in range(self.size - 1)}
        self.nid2node[0] = self.root

        self.nid2depth = dict()
        self.nid2depth[0] = 0

        self.leaves = list()
        self.height = 0
        self.leaf2ord = dict()

        self.__init_clone_depths()

        self.__init_clone_infos()

    
    def __init_clone_depths(self):
        '''
        Constructor
            initializes depths of all clones within tree from the parent vector
            used to later to initialize CloneNodes of tree
        '''

        h = 0
        for i in range(1, self.size):
            d = 1
            pid = self.parents[i-1]
            while pid != 0:
                d += 1
                pid = self.parents[pid-1]
            self.nid2depth[i] = d
            if d >= h:
                h = d
        self.height = h


    def __init_clone_infos(self):
        '''
        Constructor
            fully initializes CloneNodes for all nodes in tree
        '''

        for nid in self.nid2node.keys():
            if nid == 0:
                self.nid2node[nid].set_pid(0)
                self.nid2node[nid].set_parent(self.root)
                self.nid2node[nid].set_y(0)
            else:
                self.nid2node[nid].set_pid(self.parents[nid-1])
                self.nid2node[nid].set_depth(float(self.nid2depth[nid]/2))
        for node in self.nid2node.values():
            node.set_parent(self.nid2node[node.pid])
            parent = self.nid2node[node.pid]
            if parent != node:
                parent.add_child(node)
        for node in self.nid2node.values():
            node.set_subtree_height()
        for node in self.nid2node.values():
            node.children_ordered = sorted(node.children, key=lambda x: x.subtree_h, reverse=True)
            if len(node.children) == 0:
                self.leaves.append(node)
        self.order_leaves()
        self.set_ys_of_clones()



    def get_path_between(self, node1, node2):
        '''
        Method to determine shortest path between two nodes of tree by using least common ancestor approach
        
        :param node1: CloneNode
            starting node
        :param node2: CloneNode
            ending node
        
        :return: list
            list of CloneNode.ids with path from node1 to node2
            '''
        
        if node1 == node2:
            return []
        parent1, parent2 = node1, node2
        seen = {parent1.id, parent2.id}
        LCA =  -1
        while LCA == -1: # determine least common ancestor
                parent1 = parent1.parent
                parent2 = parent2.parent
                if parent1.id == parent2.id:
                    LCA = parent1.id
                if parent1.id in seen:
                    LCA = parent1.id
                if parent2.id in seen:
                    LCA = parent2.id
                seen.add(parent1.id)
                seen.add(parent2.id)
        #print(LCA)
        path_ford = [node1]
        path_back = [node2]
        while not (path_ford[-1].id == LCA and path_back[-1].id == LCA):
            if not path_ford[-1].id == LCA:
                path_ford.append(node1.parent)
                node1 = node1.parent
            if not path_back[-1].id == LCA:
                path_back.append(node2.parent)
                node2 = node2.parent
        #print(path_ford, path_back)
        for i in range(len(path_back)-2, -1, -1):
            path_ford.append(path_back[i])
        return path_ford
    
        

    def get_nodes_at_depth(self, d):
        '''
        Method to to retrive all nodes at specified depth within tree
        
        :param d: int
            depth of node(s)
            
        :return: list
            list of CloneNodes
        '''

        result = []
        for node in self.nid2node.values():
            if node.depth == d:
                result.append(node)
        return result
    

    def get_siblings_of(self, node):
        '''
        Method to retrive all siblisngs of specified node
        
        :param node: CloneNode
        
        :return: set
            set of CloneNodes
        '''

        parent = node.parent
        siblings = set()
        for node in self.nodes:
            if node.is_root():
                continue
            if node.parent == parent:
                siblings.add(node)
        return siblings
    

    def postorder_traversal(self):
        '''
        Return post-order traversal (ensure all children processed before parent, similar to bfs) of tree
        
        :return: list
            list of CloneNodes
        '''

        result = []
        to_visit = [self.root]
        while len(to_visit) > 0:
            curr_node = to_visit.pop()
            result.append(curr_node)
            for child in curr_node.children:
                to_visit.append(child)
        return result[::-1]
    

    def levelorder_traversal(self):
        '''
        Return level-order (process nodes level by level, from root to leaves) traversal of tree
        
        :return: list
            list of CloneNodes
        '''
       
        path = self.root
        for d in range(self.height):
            path.extend(self.get_nodes_by_depth(d))
        return path
    

    def bfs_traversal(self):
        '''
        Return bfs traversal of tree (more efficient than level_order traversal)
        
        :return: list
            list of CloneNodes
        '''

        to_visit = [self.root]
        path = []
        while len(to_visit) > 0:
            curr_node = to_visit.pop()
            path.append(curr_node)
            for child in curr_node.children:
                    to_visit.append(child)
        return path
    
  

    def order_leaves(self): 
        '''
        Method to order all leaves in ladderized tree, used to evenly space out 
        leaves in the vertical direction when plotting 

        '''

        done = False
        d = -1
        curr_node = self.root
        while not done:
            d += 1
            while not curr_node in self.leaves and curr_node.children_ordered:
                curr_node = curr_node.children_ordered.pop(0)
            if curr_node in self.leaves:
                self.leaf2ord[curr_node.id] = d
            else:
                d -= 1
            if curr_node.is_root():
                done = True
            else:
                curr_node = curr_node.parent


    def set_ys_of_clones(self):
        '''
        Method to determine y-coords of all nodes within tree, working from leaves to root
        '''

        for leaf in self.leaves:
            leaf.set_y(self.leaf2ord[leaf.id] * 30) # 30 is circle size of ndoe in tree plot
        for node in self.postorder_traversal():
            parent = node.parent
            y = 0
            for child in parent.children:
                y += child.y
            y /= len(parent.children)
            parent.set_y(y)





    

