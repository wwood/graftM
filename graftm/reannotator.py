from skbio.tree._tree import TreeNode

class TreeParaphyleticException(Exception): pass

class Reannotator:
    def _reroot_tree_by_old_root(self, old_tree, new_tree):
        '''reroot the new tree so that it matches the old tree's root, if 
        possible.'''
        
        # make a list of the left and right leaf names that are in the new tree
        if len(old_tree.children) != 2: raise Exception("Expected a binary tree, root node has number of children != 2")
        new_tip_names = set([tip.name for tip in new_tree.tips(include_self=False)])
        old_left_tip_names = [tip.name for tip in old_tree.children[0].tips(include_self=False) if tip.name in new_tip_names]
        old_right_tip_names = [tip.name for tip in old_tree.children[1].tips(include_self=False) if tip.name in new_tip_names]
        
        # find the LCA of the lefts and the rights, giving up at the root.
        left_lca = new_tree.lowest_common_ancestor(old_left_tip_names)
        right_lca = new_tree.lowest_common_ancestor(old_right_tip_names)
        
        # if both didn't LCA before hitting the root, tree paraphyletic
        # take the first one where the LCA was hit, reroot here.
        # find the LCA of the other in the other half of the tree.
        # if failed before getting to the root, then tree paraphyletic
        # reroot on one side of the internal branch between the two LCAs
        if left_lca == new_tree:
            if right_lca == new_tree:
                raise TreeParaphyleticException("Tree paraphyletic case #1")
            else:
                new_tree = self._reroot_workaround(new_tree, right_lca)
                new_lca = new_tree.lowest_common_ancestor(old_left_tip_names)
        else:
            new_tree = self._reroot_workaround(new_tree, left_lca)
            new_lca = new_tree.lowest_common_ancestor(old_right_tip_names)
            
        print new_tree.ascii_art()
        if new_lca.parent is None:
            raise TreeParaphyleticException("Tree paraphyletic case #2")
        print str(new_tree)
        far_node = self._find_longest_internal_branch_node(new_tree, new_lca)
        new_tree = self._reroot_workaround(new_tree, far_node)
        
        return new_tree
        
    def _find_longest_internal_branch_node(self, root_node, node):
        '''return the node that has the longest branch length between the given
        node and the root'''
        max_length = -1
        max_node = None
        while node.parent is not None:
            print "testing %s" % str(node)
            if node.length > max_length:
                max_node = node
                max_length = node.length
            node = node.parent
        return max_node
    
    def _reroot_workaround(self, tree, node):
        '''reroot the tree by adding a new node and removing the old root
        node, preserving distances'''
        if node == tree:
            raise Exception("Unexpected rooting: cannot reroot by current root")
        elif node in tree.children:
            # trivial case, just move the root across the branch
            sisters = [n for n in tree.children if n != node]
            if len(sisters) != 1: raise Exception("Unexpectedly found num sisters != 1")
            sister = sisters[0]
            new_root = TreeNode(children=[node, sister], parent=None)
            sister.length = sister.length + node.length
            sister.parent = new_root
            node.length = 0.0
            node.parent = new_root
            return new_root
        else:
            # regular case, new root is down the tree somewhere
            print tree.ascii_art()
            an_old_child = tree.children[0]
            OLD_CHILD_NAME = 'ochild'
            an_old_child.name = OLD_CHILD_NAME
            print "Rooting at %s" % str(node)
            
            # create new root and reroot there
            new_root = TreeNode(length=0.0, children=[node], parent=node.parent)
            new_root.parent.children = [n for n in new_root.parent.children if n != node]+[new_root]
            new_root = tree.root_at(new_root)
            
            # remove leftover dummy root
            old_child_node = new_root.find(OLD_CHILD_NAME)
            nodes_to_remove = [n for n in [old_child_node,old_child_node.parent] if len(n.children) == 1]
            if len(nodes_to_remove) != 1: raise Exception("Unexpected number of 1 child nodes")
            n = nodes_to_remove[0]
            child = n.children[0]
            grandchildren = child.children
            if len(grandchildren) != 2:
                raise Exception("Unexpectedly found num grandchildren of dummy not != 2")
            for g in grandchildren:
                g.parent = n
            n.length = n.length + child.length
            n.children = grandchildren
            n.name = None
            print new_root.ascii_art()
    
            return new_root
        
