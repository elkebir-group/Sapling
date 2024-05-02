class Tree:
    def __init__(self, edges):
        self.n = len(edges)
        edges = list(edges)
        edges.sort()
        self.edges = tuple(edges)
        self.vertices = list(set([i for e in self.edges for i in e]))
        self.vertices.sort()
        label = {v:i for i,v in enumerate(self.vertices[1:])}
        label[-1] = -1
        self.edges_relabel_from_zero = tuple((label[e[0]],label[e[1]]) for e in self.edges)
        self.children = [[] for _ in range(self.n+1)]
        
        for e in self.edges_relabel_from_zero:
            self.children[e[0]].append(e[1])
        
        self.subtree_size = [None for _ in range(self.n)]
        stack = [-1]
        r_order = []
        while len(stack)>0:
            top = stack.pop(-1)
            r_order.append(top)
            for i in self.children[top]:
                stack.append(i)
        for i in r_order[:0:-1]:
            self.subtree_size[i]=1+sum([self.subtree_size[j] for j in self.children[i]])

    def expand(self, i, allow_GL_multi,filters=None):
        """
        add a node i to tree
        filters: optional arg, a set of allowed edges
        allow_GL_multi: allow multiple children from GL, default:False
        """
        expand_trees = []
        for p in self.vertices:
            if (filters is not None) and ((p,i) not in filters):
                continue
            if p==-1 and not allow_GL_multi:
                continue
            __tmp = list(self.edges)
            __tmp.append((p,i))
            expand_trees.append(Tree(__tmp))
        for e in self.edges:
            if (filters is not None) and ((e[0],i) not in filters or (i,e[1]) not in filters):
                continue
            __tmp = list(self.edges)
            __tmp.remove(e)
            __tmp.append((e[0],i))
            __tmp.append((i,e[1]))
            expand_trees.append(Tree(__tmp))
        return expand_trees
        