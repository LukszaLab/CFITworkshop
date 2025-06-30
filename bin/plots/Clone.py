class Clone():
    """
        Node object stores item and references to its children. parent may have any number of children, but each child has only one parent."
    """

    def __init__(self, id, depth=None):
        
        """
        Params
        ----------
        id : int
                clone id to be displayed on node
        depth : int
                distance from root, used as x-coordinate when plotting node
                set to 0 when initializing root node
        """

        self.id = id
        self.depth= depth
        self.parent = None
        self.pid = None # to be set for all nodes --> int id of parent node (nid2node dictionary part of Tree object to look up actual Node)
        self.y = 0 # to be computed later
        self.children = list()
        self.children_ordered = list()
        self.subtree_h = 0

    def __repr__(self):
        return "<Node id={0}>".format(self.id)
    
    def __eq__(self, other):
        return (self.id == other.id)
    
    def __lt__(self, other):
        return (self.depth < other.depth)
    
    def __hash__(self):
        return hash(str(self.id))
    
    def set_depth(self, d):
        self.depth = d
    
    def set_depth1(self):
        if self.is_root():
            return 
        else:
            d = 1
            curr_node = self.parent
            while not curr_node.is_root():
                print(curr_node.id, d)
                d += 1
                curr_node = curr_node.parent
            self.depth = d
    

    def set_subtree_height(self):
        if self.is_leaf():
            self.subtree_h = 0
        else:
            max_depth = -1
            to_visit = [self]
            while to_visit:
                curr_node = to_visit.pop()
                if curr_node.depth > max_depth:
                    max_depth = curr_node.depth
                for child in curr_node.children:
                    to_visit.append(child)
            self.subtree_h = (max_depth - self.depth)


    def get_num_leaves_subtree(self):
        to_visit = [self]
        leaves = 0
        while to_visit:
            curr_node = to_visit.pop()
            if curr_node.is_leaf():
                leaves += 1
            for child in curr_node.children:
                to_visit.append(child)
        return leaves 


    def add_child(self, node):
        self.children.append(node)
    
    def add_children(self, children=[]):
        for child in children:
            self.add_child(child)
    
    def set_parent(self, parent):
        self.parent = parent

    def set_pid(self, pid):
        self.pid = pid
          
    def set_y(self, y):
        self.y = y


    def is_leaf(self):
        return len(self.children) == 0
    
    def is_root(self):
        return self.id == self.parent.id
    
    
    def to_list(self):
        self.children = list(self.children)


    def get_leaves(self):
        if not self.children:
            return [self]
        else:
            leaves = []
            for child in self.children:
                leaves += child.get_leaves()
            return

