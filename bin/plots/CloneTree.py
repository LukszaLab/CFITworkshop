from Clone import Clone

class CloneTree():

    def __init__(self, root, nodes, leaves, h):
        """
        Params
        ----------
        root : class Node
                Node object with None parent that is root of tree
        
        nodes: set 
                set of Node objects with necessary child/parent connections for tree

        h: int
                height of tree (i.e. depth of deepest node in nodes)
        """


        self.root = root 
        self.nodes = nodes 
        self.leaves = leaves
        self.height = h 
        self.leaf2ord = dict()
        self.nid2node = dict()
        

    def get_path_between(self, node1, node2):
        if node1 == node2:
            return []
        parent1, parent2 = node1, node2
        seen = {parent1.id, parent2.id}
        LCA =  -1
        while LCA == -1:
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
        

    def get_node_by_id(self, id):
        for node in self.nodes:
            if node.id == id:
                return node
        else:
            return None
        

    def get_nodes_at_depth(self, d):
        result = []
        for node in self.nodes:
            if node.depth == d:
                result.append(node)
        return result
    

    def get_siblings_of(self, node):
        parent = node.parent
        siblings = set()
        for node in self.nodes:
            if node.is_root():
                continue
            if node.parent == parent:
                siblings.add(node)
        return siblings
    

    def postorder_traversal(self):
        result = []
        to_visit = [self.root]
        while len(to_visit) > 0:
            curr_node = to_visit.pop()
            result.append(curr_node)
            for child in curr_node.children:
                to_visit.append(child)
        return result[::-1]
    

    def levelorder_traversal(self):
        path = self.root
        for d in range(self.height):
            path.extend(self.get_nodes_by_depth(d))
        return path
    

    def bfs_traversal(self):
        to_visit = [self.root]
        path = []
        while len(to_visit) > 0:
            curr_node = to_visit.pop()
            path.append(curr_node)
            for child in curr_node.children:
                    to_visit.append(child)
        return path
    

    def to_list(self):
        for node in self.nodes:
            node.to_list()
  

    def order_leaves(tree): # ladderize clonetree (leaves ordered from deepest to most shallow)
        leafid2ord = dict()
        done = False
        d = -1
        curr_node = tree.root
        while not done:
            d += 1
            while not curr_node in tree.leaves and curr_node.children_ordered:
                curr_node = curr_node.children_ordered.pop(0)
            if curr_node in tree.leaves:
                leafid2ord[curr_node.id] = d
            else:
                d -= 1
            if curr_node.is_root():
                done = True
            else:
                curr_node = curr_node.parent
        tree.leaf2ord = leafid2ord


    def set_ys_of_clones(self):
        for leaf in self.leaves:
            leaf.set_y(self.leaf2ord[leaf.id] * 30)
        for node in self.postorder_traversal():
            parent = node.parent
            y = 0
            for child in parent.children:
                y += child.y
            y /= len(parent.children)
            parent.set_y(y)





    

