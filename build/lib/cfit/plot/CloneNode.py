class CloneNode():
    """
        CloneNode object stores item and references to its children. 
        Parent clone may have any number of children, but each child has only one parent.
    
    Attributes:

        id: int
            clone number

        depth: int
            distance from clone 0 (root), sed as x-coord to plot CloneNode in scatter plot for horizontal tree

        parent: class CloneNode
            pointer to parent CloneNode object

        pid: int
            id of parent CloneNode
        
        y: int
            y-coord used to plot CloneNode in scatter plot for horizontal tree 
        
        chidlren: list
            list of children CloneNode objects 

        chidlren_orderd: list
            sorted list of CloneNodes, from node with largest subtree_height to 
            node tith smallest subtree_height

        subtree_height: int
            height of subtree with CloneNode as root, used to sort children for children_ordered attirbute
    
    """


    def __init__(self, id, depth=None):
        
        """
        Constructor

        :param id: int
                clone id to be displayed on node

        :param depth: int
                distance from root, used as x-coordinate when plotting node
                set to 0 when initializing root node
        """

        self.id = id
        self.depth= depth
        self.parent = None
        self.pid = None 
        self.y = 0 # to be computed later
        self.children = list()
        self.children_ordered = list()
        self.subtree_h = 0


    def __repr__(self):
        return "<Clone id={0}>".format(self.id)
    
    def __eq__(self, other):
        return (self.id == other.id)
    
    def __lt__(self, other):
        return (self.depth < other.depth)
    
    def __hash__(self):
        return hash(str(self.id))
    
    def set_depth(self, d):
        '''
        :param d: int
        '''
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
    
    def set_parent(self, parent):
        '''
        :param parent: CloneNode
        '''

        self.parent = parent


    def set_pid(self, pid):
        '''
        :param pid: int
        '''

        self.pid = pid
          

    def set_y(self, y):
        '''
        :param y: int
        '''
        self.y = y
    

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


    def add_child(self, node):
        '''
        :param node: CloneNode
        '''
        self.children.append(node)
    

    def add_children(self, children=[]):
        '''
        :param chidlren: list
            list of CloneNode objects
            '''
        for child in children:
            self.add_child(child)


    def is_leaf(self):
        return len(self.children) == 0
    
    def is_root(self):
        return self.id == self.parent.id
  

