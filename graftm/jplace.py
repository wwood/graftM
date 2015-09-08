from skbio.tree import TreeNode
from StringIO import StringIO
import re

class NewickFormatError(Exception): pass

class JPlaceTree:
    def __init__(self, newickish_string):
        self.tree, self._edge_indices = self._jplace_tree_to_tree_node(StringIO(newickish_string))
        
    def child_node(self, edge_index):
        '''Return the node that is on the leaf side of the given edge index'''
        return self.edge_indices[edge_index]
    
    def _jplace_tree_to_tree_node(self, fh, convert_underscores=True):
        '''Most of this code copied from skibio.io.format._newick_to_tree_node
        except that the extra edge identifiers e.g. {0} are parsed.
        
        Returns
        -------
        list of two elements:
        1. a TreeNode parsed as if the tree was regular newick format
        2. list of pointers to TreeNodes. The position in the list is the 
        edge number, and the TreeNode pointed to is the TreeNode on the child
        side of the edge.
        '''
        
        tree_stack = []
        current_depth = 0
        last_token = ''
        next_is_distance = False
        root = TreeNode()
        tree_stack.append((root, current_depth))
        edge_indices = []
        edge_regex = re.compile('^([\d\.]+)\{(\d+)\}$')
        for token in self._tokenize_newick(fh, convert_underscores=convert_underscores):
            # Check for a label
            if last_token not in '(,):':
                if not next_is_distance:
                    tree_stack[-1][0].name = last_token if last_token else None
                else:
                    next_is_distance = False
            # Check for a distance
            if token == ':':
                next_is_distance = True
            elif last_token == ':':
                child_node = tree_stack[-1][0]
                er = edge_regex.match(token)
                if er is None:
                    raise NewickFormatError("Edge label appears to be missing from token %s" % token)
                try:
                    child_node.length = float(er.group(1))
                except ValueError:
                    raise NewickFormatError("Could not read length as numeric type"
                                            ": %s." % token)
                edge_index = int(er.group(2))
                if edge_indices[edge_index]:
                    raise NewickFormatError("It appears the edge index %i from token %s was found multiple times" % (edge_index, token))
                edge_indices[edge_index] = child_node
            elif token == '(':
                current_depth += 1
                tree_stack.append((TreeNode(), current_depth))
            elif token == ',':
                tree_stack.append((TreeNode(), current_depth))
            elif token == ')':
                if len(tree_stack) < 2:
                    raise NewickFormatError("Could not parse file as newick."
                                            " Parenthesis are unbalanced.")
                children = []
                # Pop all nodes at this depth as they belong to the remaining
                # node on the top of the stack as children.
                while current_depth == tree_stack[-1][1]:
                    node, _ = tree_stack.pop()
                    children.insert(0, node)
                parent = tree_stack[-1][0]
                if parent.children:
                    raise NewickFormatError("Could not parse file as newick."
                                            " Contains unnested children.")
                # This is much faster than TreeNode.extend
                for child in children:
                    child.parent = parent
                parent.children = children
                current_depth -= 1
            elif token == ';':
                if len(tree_stack) == 1:
                    return root, edge_indices
                break
    
            last_token = token
    
        raise NewickFormatError("Could not parse file as newick."
                                " `(Parenthesis)`, `'single-quotes'`,"
                                " `[comments]` may be unbalanced, or tree may be"
                                " missing its root.")
        
    def _tokenize_newick(self, fh, convert_underscores=True):
        '''directly copied from skbio. Only copied because it is an underscore
        method and therefore potentially susceptible to change without notice
        there'''
        structure_tokens = set('(),;:')
        not_escaped = True
        label_start = False
        last_non_ws_char = ''
        last_char = ''
        comment_depth = 0
        metadata_buffer = []
        # Strategy:
        # We will iterate by character.
        # Comments in newick are defined as:
        # [This is a comment]
        # Nested comments are allowed.
        #
        # The following characters indicate structure:
        #      ( ) , ; :
        #
        # Whitespace is never allowed in a newick label, so an exception will be
        # thrown.
        #
        # We use ' to indicate a literal string. It has the highest precedence of
        # any operator.
        for line in fh:
            for character in line:
                # We will start by handling the comment case.
                # This code branch will probably never execute in practice.
                # Using a comment_depth we can handle nested comments.
                # Additionally if we are inside an escaped literal string, then
                # we don't want to consider it a comment.
                if character == "[" and not_escaped:
                    # Sometimes we might not want to nest a comment, so we will use
                    # our escape character. This is not explicitly mentioned in
                    # any format specification, but seems like what a reasonable
                    # person might do.
                    if last_non_ws_char != "'" or comment_depth == 0:
                        # Once again, only advance our depth if [ has not been
                        # escaped inside our comment.
                        comment_depth += 1
                if comment_depth > 0:
                    # Same as above, but in reverse
                    if character == "]" and last_non_ws_char != "'":
                        comment_depth -= 1
                    last_non_ws_char = character
                    continue
                # We are not in a comment block if we are below here.
    
                # If we are inside of an escaped string literal, then ( ) , ; are
                # meaningless to the structure.
                # Otherwise, we are ready to submit our metadata token.
                if not_escaped and character in structure_tokens:
                    label_start = False
                    metadata = ''.join(metadata_buffer)
                    # If the following condition is True, then we must have just
                    # closed a literal. We know this because last_non_ws_char is
                    # either None or the last non-whitespace character.
                    # last_non_ws_char is None when we have just escaped an escape
                    # and at the first iteration.
                    if last_non_ws_char == "'" or not convert_underscores:
                        # Make no modifications.
                        yield metadata
                    elif metadata:
                        # Underscores are considered to be spaces when not in an
                        # escaped literal string.
                        yield metadata.replace('_', ' ')
                    # Clear our buffer for the next metadata token and yield our
                    # current structure token.
                    metadata_buffer = []
                    yield character
                # We will now handle escaped string literals.
                # They are inconvenient because any character inside of them is
                # valid, especially whitespace.
                # We also need to allow ' to be escaped by '. e.g. '' -> '
                elif character == "'":
                    not_escaped = not not_escaped
                    label_start = True
                    if last_non_ws_char == "'":
                        # We are escaping our escape, so it should be added to our
                        # metadata_buffer which will represent some future token.
                        metadata_buffer.append(character)
                        # We do not want a running chain of overcounts, so we need
                        # to clear the last character and continue iteration from
                        # the top. Without this, the following would happen:
                        # ''' ' -> '' <open literal>
                        # What we want is:
                        # ''' ' -> '<open literal> <close literal>
                        last_non_ws_char = ''
                        last_char = ''
                        continue
    
                elif not character.isspace() or not not_escaped:
                    if label_start and last_char.isspace() and not_escaped:
                        raise NewickFormatError("Newick files cannot have"
                                                " unescaped whitespace in their"
                                                " labels.")
                    metadata_buffer.append(character)
                    label_start = True
    
                # This is equivalent to an `else` however it prevents coverage from
                # mis-identifying the `continue` as uncalled because cpython will
                # optimize it to a jump that is slightly different from the normal
                # jump it would have done anyways.
                elif True:
                    # Skip the last statement
                    last_char = character
                    continue
    
                last_char = character
                # This line is skipped in the following cases:
                #    * comment_depth > 0, i.e. we are in a comment.
                #    * We have just processed the sequence '' and we don't want
                #      the sequence ''' to result in ''.
                #    * We have encountered whitespace that is not properly escaped.
                last_non_ws_char = character