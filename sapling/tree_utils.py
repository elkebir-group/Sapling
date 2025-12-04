from itertools import chain, combinations
from collections import defaultdict

def powerset(iterable):
    """
    Generates the power set of an iterable, which is the set of all possible subsets.

    The power set includes all combinations of the input iterable's elements, ranging from
    the empty set to the full set. Each subset is returned as a tuple.

    Args:
        iterable (iterable): An iterable (e.g., list, tuple, or set) for which to generate
                             the power set.

    Returns:
        itertools.chain: An iterable of tuples, where each tuple represents a subset of the
                        input iterable. The subsets are ordered by increasing size.

    Example:
        >>> list(powerset([1, 2, 3]))
        [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

class Tree:
    """
    A class representing a tree structure, used for modeling hierarchical relationships
    such as phylogenetic trees or other tree-based data structures.

    The tree is defined by its edges, and it supports two modes of expansion:
    - `big_expand`: Generates all possible expansions by adding a new node and reconfiguring
      the tree structure.
    - `small_expand`: Generates a limited set of expansions by adding a new node in specific
      positions.

    Attributes:
        edges (tuple): A tuple of edges representing the tree structure.
        vertices (list): A sorted list of unique vertices in the tree.
        edges_relabel_from_zero (tuple): A relabeled version of edges where vertices are
                                         reindexed starting from 0.
        children (list): A list of lists representing the children of each vertex.
        subtree_size (list): A list storing the size of the subtree rooted at each vertex.
        n (int): The number of edges in the tree.
    """
    def __init__(self, edges, use_big_expand):
        """
        Initializes a Tree object.

        Args:
            edges (list): A list of tuples representing the edges of the tree.
            use_big_expand (bool): If True, the tree will use `big_expand` for expansion;
                                   otherwise, it will use `small_expand`.
        """
        self.n = len(edges)
        if use_big_expand:
            self.expand = self.big_expand
        else:
            self.expand = self.small_expand
            
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

    def big_expand(self, i, poly_clonal_root):
        """
        Expands the tree by adding a new node `i` and generating all possible configurations
        by reconfiguring the children of existing nodes.

        Args:
            i (int): The new node to be added to the tree.
            poly_clonal_root (bool): If True, allows multiple children from the root node (-1).

        Returns:
            list: A list of Tree objects representing all possible expansions.
        """
        children_dict = defaultdict(list)
        for (u,v) in self.edges:
            children_dict[u].append(v)

        expand_trees = []
        for p in self.vertices:
            # enumerate all subsets of p's children
            if p == -1 and not poly_clonal_root: continue
            for subset in powerset(children_dict[p]):
                __tmp = list(self.edges)
                __tmp.append((p,i))
                for j in subset:
                    __tmp.remove((p,j))
                for j in subset:
                    __tmp.append((i,j))
                expand_trees.append(Tree(__tmp, True))
        return expand_trees
    
    def small_expand(self, i, poly_clonal_root,filters=None):
        """
        Expands the tree by adding a new node `i` in specific positions, optionally filtered
        by a set of allowed edges.

        Args:
            i (int): The new node to be added to the tree.
            poly_clonal_root (bool): If True, allows multiple children from the root node (-1).
            filters (set, optional): A set of allowed edges. If provided, only edges in this
                                    set will be considered for expansion.

        Returns:
            list: A list of Tree objects representing the valid expansions.
        """
        expand_trees = []
        for p in self.vertices:
            if (filters is not None) and ((p,i) not in filters):
                continue
            if p==-1 and not poly_clonal_root:
                continue
            __tmp = list(self.edges)
            __tmp.append((p,i))
            expand_trees.append(Tree(__tmp, False))
        for e in self.edges:
            if (filters is not None) and ((e[0],i) not in filters or (i,e[1]) not in filters):
                continue
            __tmp = list(self.edges)
            __tmp.remove(e)
            __tmp.append((e[0],i))
            __tmp.append((i,e[1]))
            expand_trees.append(Tree(__tmp, False))
        return expand_trees
        