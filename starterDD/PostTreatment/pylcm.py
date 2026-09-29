# -*- coding: utf-8 -*-
# __author__ = B. Daré
# __date__ = 2026-04-17

import numpy as np          

class Pylcm:
    """
    
    This class allows the user to represent LCM data structures, used in the 
    codes DRAGON and DONJON. A LCM data structure can be interpreted as an 
    associative table, containing heterogeneous elements by associating a 
    12-character key and a content, or as a list, containing heterogeneous 
    contents identified by their index.

    A LCM data structure can contain other LCM data structures (associative 
    tables or lists), or elementary objects:
    
    - integer
    - single or double precision real number
    - character string
    - logical number (0 or 1)
    - complex number
        
    Note that an elementary content in a LCM data structure can be an 
    elementary object alone, or a set of several elementary objects of the 
    same type.

    Those properties allow consequently tree structures.

    """

    def __init__(self, name='Pylcm_object', level=0, type='table', size=np.inf, parent=None):
        if len(name) > 12:
            raise ValueError("Name must be a 12-character string")
        self.name = name # 12-character name
        self.content = {}
        self.level = level # level in the global tree
        self.type = type # 'table' or 'list'
        self.size = size # np.inf if table, len(list) if list
        self.parent = parent # Pylcm parent, None if root


    def add(self, key: str, object):
        if self.type == 'list':
            self.content[int(key)] = object
        else:
            self.content[key] = object


    def get(self, key: str):
        if self.type == 'list':
            return self.content[int(key)]
        else:
            return self.content[key]


    def print_tree(self, path: str, name: str=None):
        """
        
        This method creates a txt file at the specified *path* representing the
        tree structure of a ``Pylcm`` object (and that of any other child 
        ``Pylcm`` object contained).

        Args:
            *path* (str): The path to the folder where the output file will 
                          be created.
            *name* (str): Optional, name of the output file. If None, the 
                          function will give the name of the ``Pylcm`` object.

        Returns:
            ``None``
        
        """

        # Define the name of the file
        if name is None:
            name = self.name

        # Create the output file
        with open(f"{path}/{name}.txt", 'w') as f:
            f.write(self.name.ljust(14,'_'))
            self._print_level(f, self.level, [])
                        

    def _print_level(self, f, root_level, state_list):
        """
        
        Internal method used by the method ``print_tree`` to recursively 
        generate the tree structure in the output file by browsing the ``Pylcm`` 
        object.

        Args:
            *f*: The file object where to write.
            *root_level*: The level of the root, in other words the level of the
                          first ``Pylcm`` object on which the function 
                          ``print_tree`` is called.
            *state_list*: List of boolean values used to determine if at a 
                          specified level (linked to the index in the list) all 
                          the branches of the current node have been explored 
                          (True) or not (False). It is used for display purpuse 
                          only.
        
        Returns:
            ``None``
        
        """

        if (self.level - root_level) < len(state_list):
            state_list[self.level - root_level] = False
        else:
            state_list.append(False)
        children = list(self.content.keys())
        for i in range(len(children)):
            if i == 0:
                ch = '_'
            else:
                ch = '|'
                f.write(' '*14)
                for j in range(self.level - root_level):
                    if not state_list[j]:
                        f.write('|'.ljust(17,' '))
                    else:
                        f.write(' '*17)
            if i == len(children)-1:
                state_list[self.level - root_level] = True
            if isinstance(self.get(children[i]), Pylcm):
                f.write(ch + f"__{children[i]}".ljust(16,'_'))
                self.get(children[i])._print_level(f, root_level, state_list)
            else:
                f.write(ch + f"__{children[i]}\n")
                # f.write(' '*14)
                # for j in range(self.level - root_level +1):
                #     if not state_list[j]:
                #         f.write('|'.ljust(17,' '))
                #     else:
                #         f.write(' '*17)
                # f.write('\n')


    @staticmethod
    def _ascii_reader(path: str) -> list:
        """
        
        This method reads an ASCII file at the specified *path* containing the 
        information of a LCM data structure. It formats the data so as to be 
        used by the method ``generate_from_ascii``.

        Args:
            *path* (str): The path to the ASCII file where the output file will 
                          be created.

        Returns:
            ``list``: A list of strings containing the node information 
                      "*tuple* <- *content*" from an ASCII file. The list is 
                      reversed to work as a stack.
        
        """    
        with open(path, 'r') as f:
            file = f.read()
            nodes_list = file[2:].split('->') # Separate the nodes information with '->'. A node is in the form "<tuple> <- <content>".
            nodes_list.reverse() # reverse the list to work as a stack
        return nodes_list


    @classmethod
    def generate_from_ascii(cls, path: str,**kwargs):
        """

        This method generates a ``Pylcm`` object (associative table or list) and 
        its content from an ASCII file located at the given *path*. In the case 
        of a tree structure, for example a multicompo directory from a DRAGON 
        run, the function will automatically and recursively generate the 
        associative tables and lists (as ``Pylcm`` objects) contained in the 
        tree.
        
        Args:
            *path* (str): The path to the ASCII file to read.
            ***kwargs* (dict): Dictionary of optional arguments used for the 
                               internal recursive calls. Could be the optional
                               instance attributes.

        Returns:
            ``Pylcm``: The generated ``Pylcm`` object.

        Details on the structure of the ASCII file:
        -------------------------------------------
        The ASCII file contains the elements of a LCM data structure.
            
        Each element or 'node' (an associative table, a list, or any elementary 
        content) is represented in the ASCII file by a pair tuple/content 
        implemented:
            
        "-> *level*  *name_type*  *content_type*  *content_size* <-
         *content*"
            
        Where:
            - *level* (int): the level of the node in the global tree 
                             structure, 0 is the root, 1 the first children 
                             and so on.
            - *name_type* (int): if I understand correctly, 12 means a 
                                 12-character key and 0 means no key.
            - *content_type* (int): 
                                * 0 --> associative table, 
                                * 1 --> integer, 
                                * 2 --> single precision real number, 
                                * 3 --> character*4 data, 
                                * 4 --> double precision real number, 
                                * 5 --> logical number (0 or 1), 
                                * 6 --> complex number, 
                                * 10 --> list, 
                                * 99 --> undefined.
            - *content_size* (int): -1 for a associative table, N for a list 
                                    of size N or for N elementary contents of 
                                    the same type, 0 if the content doesn't 
                                    exist.

        To "close" an associative table (or to come back to the parent 
        structure), the ASCII file contains an empty node with a tuple like 
        this:
            
        "-> -*current level* 0 0 0 <-" 
        
        """
        if len(kwargs.get('name', 'Pylcm_object')) > 12:
            raise ValueError("Name must be a 12-character string")
        pylcm_object = cls(name=kwargs.get('name', 'Pylcm_object'), level=kwargs.get('level', 0), type=kwargs.get('type', 'table'), size=kwargs.get('size', np.inf), parent=kwargs.get('parent', None))
        nodes_list = kwargs.get('nodes_list', Pylcm._ascii_reader(path)) # For the recursive construction, check if the nodes_list has been already extracted.
        # The ASCII file extraction and formatting is performed with the function ascii_reader.
        # The nodes_list is a list of strings containing the node information "<tuple> <- <content>" from an ASCII file.
        # Pay attention to the fact that the list needs to be reversed to work as a stack.

        # Formatting the node information
        node = nodes_list.pop()
        info, content = node.replace('\n', '').split('<-')
        level, name_type, content_type, content_size = map(int, info.split())
        iteration = 1

        while level >= 0 and iteration <= pylcm_object.size:
            # Content type cases
            # If necessary, the cases logical number (=5) and complex number (=6) can be implemented as well
            if content_type == 0: # if associative table
                pylcm_object.add(content.strip(), cls.generate_from_ascii(path, nodes_list=nodes_list, level=level, parent=pylcm_object))
                if name_type == 12:
                    pylcm_object.get(content.strip()).name = content.strip()
            elif content_type == 10: # if list
                pylcm_object.add(content.strip(), cls.generate_from_ascii(path, nodes_list=nodes_list, level=level, type='list', size=content_size, parent=pylcm_object))
                if name_type == 12:
                    pylcm_object.get(content.strip()).name = content.strip()
            elif content_type == 1: # if integer
                if content_size == 1:
                    key, object = content.split()
                    pylcm_object.add(key, int(object))
                else:
                    elements = content.split()
                    key, object = elements[0], np.array(list(map(int, elements[1:])))
                    pylcm_object.add(key, object)
            elif content_type == 2 or content_type == 4: # if float (single precision = 2) or double (double precision = 4)
                if content_size == 1:
                    key, object = content.split()
                    pylcm_object.add(key, float(object))
                else:
                    elements = content.split()
                    key, object = elements[0], np.array(list(map(float, elements[1:])))
                    pylcm_object.add(key, object)
            elif content_type == 3: # if string
                key, object = content.split('         4'*content_size)
                pylcm_object.add(key.strip(), object.strip())
            elif content_type == 99: # if undefined
                pylcm_object.add(content.strip(), None)

            iteration += 1
            if iteration <= pylcm_object.size:
                node = nodes_list.pop()
                info, content = node.replace('\n', '').split('<-')
                level, name_type, content_type, content_size = map(int, info.split())
            
        
        return pylcm_object
    

    def __repr__(self):
        if self.type == 'list':
            return f"Pylcm(name={self.name!r}, type={self.type!r}, size={self.size!r}, parent={self.parent!r})"
        else:
            return f"Pylcm(name={self.name!r}, type={self.type!r}, parent={self.parent!r})"