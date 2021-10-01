#q1b.py - attempt 2
"""
UPGMA - Unweighted pair group method using arithmetic averages
Takes as input the distance matrix of species as a numpy array
Returns tree either as a dendrogram or in Newick format

Resources reffered to and used :

1. https://stackoverflow.com/questions/58532758/read-csv-file-into-list-and-convert-the-strings-into-integers-python 
2. http://etetoolkit.org/treeview/ 
3. https://www.geeksforgeeks.org/writing-data-from-a-python-list-to-csv-row-wise/ 
4. https://stackoverflow.com/questions/56025046/2d-list-to-csv-by-column 
5. https://en.wikipedia.org/wiki/UPGMA 
6. https://en.wikipedia.org/wiki/Phylogenetic_tree 
7. https://en.wikipedia.org/wiki/Child_node 
8. https://en.wikipedia.org/wiki/Tree_structure 
9. https://medium.com/@sharma.ravit/upgma-method-designing-a-phylogenetic-tree-9a708de18419#:~:text=%20UPGMA%20Method%3A%20Designing%20a%20Phylogenetic%20Tree%20,1%E2%80%932%20until%20the%20tree%20is%20complete%20More%20 
10. https://bit.ly/39T1tmb 
11. https://telliott99.blogspot.com/2010/11/upgma-in-python.html 
12. https://telliott99.blogspot.com/2010/11/upgma-in-python-2.html 
13. https://telliott99.blogspot.com/2010/11/upgma-in-python-3-sarich-data.html 
14. https://telliott99.blogspot.com/2010/03/clustering-with-upgma.html 
15. https://telliott99.blogspot.com/2010/04/visualizing-upgma-clustering.html 
16. https://www.geeksforgeeks.org/draw-a-tree-using-arcade-library-in-python/ 
17. https://www.geeksforgeeks.org/binarytree-module-in-python/ 
18. https://github.com/itsjinendrajain/Coding-Ninjas-Problem-Solving-Using-Python/tree/main/8.Searching%20%26%20Sorting 
19. https://www.geeksforgeeks.org/tree-sort/ 
20. http://www.bx.psu.edu/~dcking/man/newicktree.html#:~:text=The%20Newick%20Standard%20for%20representing%20trees%20in%20computer-readable,represented%20by%20the%20following%20sequence%20of%20printable%20characters%3A 
21. https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python 
22. https://pure.mpg.de/rest/items/item_3258810_1/component/file_3273695/content 
23. http://etetoolkit.org/docs/2.3/tutorial/tutorial_trees.html 

"""

#%% import modules
import numpy as np
import itertools
import csv
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
import pandas as pd

def merge(list1, list2):
      
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list

# Python program to get transpose
# elements of two dimension list
def transpose(l1, l2):
 
    # iterate over list l1 to the length of an item
    for i in range(len(l1[0])):
        # print(i)
        row =[]
        for item in l1:
            # appending to new list with values and index positions
            # i contains index position and item contains values
            row.append(item[i])
        l2.append(row)
    return l2

# class for nodes
class Node:
    """
    Data structure to store node of a UPGMA tree
    """

    def __init__(self, left=None, right=None, up_height=0.0, down_height=0.0):
        """
        Creating a node.
        For a single taxon, set taxon name as self.left, leave right as none.
        For an operational taxonomic unit(OTU) set left and right to child nodes.

        Parameters
        ----------
        left : default = none, taxon label
        right : default = none, taxon label
        up_height : float, default = 0.0, dist to parent node, if any
        down_height : float, default = 0.0, dist to child node, if any
        """
        self.left = left
        self.right = right
        self.uh = up_height
        self.dh = down_height

    def leaves(self) -> list:
        """
        Method to find the taxa under any given node, effectively equivalent to
        finding leaves of a binary tree. Only lists original taxa and not OTUs.

        Returns a list of node names, not nodes themselves.
        """
        if self == None:
            return []
        if self.right == None:
            return [self.left]
        leaves = self.left.leaves() + self.right.leaves()
        return leaves

    def __len__(self) -> int:
        """
        Method to define len() of a node.

        Returns the number of original taxa under any given node.
        """
        return sum(1 for taxa in self.leaves())

    def __repr__(self) -> str:
        """
        Method to give readable print output
        """
        return "-".join(self.leaves())


# class for UPGMA
class UPGMA:
    def __init__(self, dist_matrix: np.ndarray, taxa: list):
        """
        Initialize an UPGMA class.
        Takes a nxn distance matrix as input. A list of n taxon id is required
        in the same order as the distance matrix row/column

        Parameters
        ----------
        dist_matrix : numpy array, distance matrix of species
        taxa : list of int or str to identify taxa
        """
        self.distances = dist_matrix
        self.taxa = taxa
        self.build_tree(self.distances, self.taxa)

    def build_tree(self, dist_matrix: np.ndarray, taxa: list) -> Node:
        """
        Method to construct a tree from a given distance matrix and taxa list.

        Parameters
        ----------
        dist_matrix : np.ndarray of pairwise distances
        taxa : list of taxa id. Elements of lists have to be unique

        Returns the root node for constructed tree.
        """
        # creating node for each taxa
        nodes = list(Node(taxon) for taxon in taxa)

        # dictionary row/column id -> node
        rc_to_node = dict([i, j] for i, j in enumerate(nodes))
        
        # dictionary taxa -> row/column id
        taxa_to_rc = dict([i, j] for j, i in enumerate(taxa))
        
        # make copy of dist matrix to work on
        work_matrix = dist_matrix
        # set all diagonal elements to infinity for ease of finding least distance
        work_matrix = np.array(work_matrix, dtype=float)
        np.fill_diagonal(work_matrix, np.inf)

        # loop
        while len(nodes) > 1:
            # finding (row, col) of least dist
            least_id = np.unravel_index(work_matrix.argmin(), work_matrix.shape, "C")
            least_dist = work_matrix[least_id[0], least_id[1]]
            # nodes corresponding to (row, col)
            node1, node2 = rc_to_node[least_id[0]], rc_to_node[least_id[1]]
            
            # add OTU with node1 and node2 as children. set heights of nodes
            new_node = Node(node2, node1)
            nodes.append(new_node)
            node1.uh = least_dist / 2 - node1.dh
            node2.uh = least_dist / 2 - node2.dh
            new_node.dh = least_dist / 2
            nodes.remove(node1)
            nodes.remove(node2)
           
            # create new working distance matrix
            work_matrix = self.update_distance(dist_matrix, nodes, taxa_to_rc)
            
            # update row/col id -> node dictionary
            rc_to_node = dict([i, j] for i, j in enumerate(nodes))
        # set tree to root
        self.tree = nodes[0]

    def update_distance(
        self, dist_matrix: np.ndarray, nodes: list, taxa_to_rc: dict
    ) -> np.ndarray:
        """
        Method to make a new distance matrix with newer node list.

        Parameters
        ----------
        dist_matrix : np.ndarray of pairwise distances for all taxa
        nodes : list of updated nodes
        taxa_to_rc : dict for taxa -> row/col id

        Returns np.ndarray of pairwise distances for updated nodes
        """
        # dictionary for node -> row/col id
        node_to_rc = dict([i, j] for j, i in enumerate(nodes))
        
        rc = len(nodes)
        new_dist_matrix = np.zeros((rc, rc), dtype=float)
        for node1 in nodes:
            row = node_to_rc[node1]
            for node2 in nodes:
                node_pairs = list(itertools.product(node1.leaves(), node2.leaves()))
                col = node_to_rc[node2]
                new_dist_matrix[row, col] = sum(
                    dist_matrix[taxa_to_rc[i], taxa_to_rc[j]] for i, j in node_pairs
                ) / len(node_pairs)
        np.fill_diagonal(new_dist_matrix, np.inf)
        return new_dist_matrix


def tree_to_newick(t) -> str:
    """
    Function to convert tree to Newick, slightly modified form of the tree.py version.
    Takes the root node of an UPGMA tree as input
    """
    if t.right == None:
        return t.left + ":" + str(t.uh)
    else:
        return (
            "("
            + ",".join([tree_to_newick(x) for x in [t.left, t.right]])
            + "):"
            + str(t.uh)
        )


#main

def main ():
#to help resolve indent errors 
    file_in ='./Ndistance.txt'
    file_in_2='./Ntaxa.txt'
    file_out="./Nnewick.txt"
    if os.access(file_in, os.R_OK):
        print ("File ", file_in ," is accessible to read")
    else:
        print ("File ", file_in ," is not accessible to read, exiting.")
        exit (1)
    reader = csv.reader(open(file_in,"r+"), delimiter=',')
    distance = list(reader)
    print(distance)
    outer_out_list = []
    for inner_list in distance:
        innet_out_list = []
        for string in inner_list:
            innet_out_list.append(float (string))
        outer_out_list.append(innet_out_list)
    print(outer_out_list)
    
    if os.access(file_in_2, os.R_OK):
        print ("File ", file_in_2 ," is accessible to read\n")
    else:
        print ("File ", file_in_2 ," is not accessible to read, exiting.")
        exit (1)    
    reader_2 = csv.reader(open(file_in_2,"r+"), delimiter=',')
    taxa = list(reader_2)
    taxa = taxa[0]
    print( taxa )
    
    distances = np.array( outer_out_list )
    x = UPGMA(distances, taxa).tree
    final_output_newick=(tree_to_newick(x))
    print(final_output_newick)
    with open(file_out,"w") as f_out:
      f_out.write(final_output_newick)
    
if __name__ == "__main__":
    main()