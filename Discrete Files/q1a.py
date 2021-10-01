
#q1a.py

'''
 After reading given input of Nucleotide.txt, I see that each new entry begins with a '>',
 followed by a line of Information on the gene, followed by some N lines of ATGC

 Identified as the FASTA format.

 I shall try to use biopython to give us the (n^2 - n)/2 values for distance matrix.

 Resources referred to 
 
 1. https://python.omics.wiki/biopython/examples/read-fasta
 2. https://www.geeksforgeeks.org/biopython-pairwise-alignment/
 3. https://stackoverflow.com/questions/20580657/how-to-read-a-fasta-file-in-python 
 4. https://vlab.amrita.edu/index.php?sub=3&brch=274&sim=1438&cnt=1 
 5. https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation#ClustalOmegaHelpandDocumentation-sequence 
 6. https://www.tutorialspoint.com/biopython/biopython_sequence_alignments.htm 
 7. https://www.geeksforgeeks.org/biopython-sequence-alignment/#:~:text=Biopython%20%E2%80%93%20Sequence%20Alignment%20Last%20Updated%20%3A%2011,identify%20the%20region%20of%20similarity%20among%20them.%20 
 8. https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#:~:text=fasta%20files%20may%20be%20submitted%20with%20corresponding%20qual%20files%2C%20too.%20these%20are%20recognized%20in%20the%20sra%20data%20processing%20pipeline%20as%20equivalent%20to%20fastq%20and%20should%20be%20specified%20as%20fastq%20when%20submitting%20the%20data%20files. 
 9. https://www.csestack.org/python-check-if-file-directory-exists/
 10. https://flexiple.com/python-create-file/ 
 11. https://www.askpython.com/python/built-in-methods/python-print-to-file#:~:text=Python%20%E2%80%93%20Print%20to%20File%201%20Method%201%3A,4%20Method%204%3A%20Use%20the%20logging%20module.%20 
 12. https://stackoverflow.com/questions/9426045/difference-between-exit0-and-exit1-in-python#:~:text=exit%20%281%29%20means%20there%20was%20some%20issue%20%2F,A%20zero%20error%20code%20means%20a%20successful%20exit. 


'''

import csv
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2

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

def main ():
#to help resolve indent errors 
    file_in ='./Nucleotide.txt'
    file_out='./Ndistance.txt'

    if os.access(file_in, os.R_OK):
        print ("File ", file_in ," is accessible to read")
    else:
        print ("File ", file_in ," is not accessible to read, exiting.")
        exit (1)
    
    if os.access(file_out, os.W_OK):
        print ("File ", file_out," is accessible to write \n")
    else:
        print ("File ", file_out," is not accessible to write, exiting.")
        exit (1)

    IdList = list ()
    NucleotideSeqList= list ()

    with open(file_out, 'w') as f_out:
    #to help resolve indent errors   
        print("Beginning to Read ",file_in,"\n")
        for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        #to help resolve indent errors 

            # remove .id from .description record (remove all before first space)
            seq_record.description=' '.join(seq_record.description.split()[1:])
            
            # do something (print or edit seq_record)       
            print('SequenceID = '  + seq_record.id)
            print('Description = ' + seq_record.description)
            print('Sequence = ' + seq_record.seq + '\n' )

            # Save the Sequence into NucleotideSeqList and Id into IdList
            NucleotideSeqList.append( seq_record.seq )
            IdList.append( seq_record.id )

        no_of_nucleotides = len(NucleotideSeqList )
        print(no_of_nucleotides,'Nucleotides scanned into our List\n')
        print("Our Lists are : ",merge ( IdList, NucleotideSeqList ))
        print(NucleotideSeqList)
        print(IdList,"\n")

        #Initialize distance matrix
        dist_Matrix=list()

        #Now to find distances, pairwise 
        for i in range(0,no_of_nucleotides,1):
            column= list ()
            for i_1 in range(0,i+1,1):
                column.append(0)
            for j in range (i+1, no_of_nucleotides, 1):
                #compare with j'th nucleotide, and find best pairing.
                alignments = pairwise2.align.globalxx(NucleotideSeqList[i], NucleotideSeqList[j])
                length = ( max(len(NucleotideSeqList[i]), len(NucleotideSeqList[j])) // 1)
                score= 0
                for alignment in alignments: 
                    # print( (alignment[2]//1) , (alignment[4]//1) , ( (alignment[2]//1)/ (alignment[4]//1) ) ) 
                    if ( (alignment[2]//1) >  score ):
                        score = (alignment[2]//1)
                norm_score = score/length
                print(score,length,norm_score)
                #add to 'column' IMPORTANT LINE - CHANGE THE "score" to "norm_score" to get normalized ans.
                column.append(norm_score)
            dist_Matrix.append(column)
        out_Matrix=list ()
        out_Matrix = transpose (dist_Matrix, out_Matrix)
        print(out_Matrix)
        writer = csv.writer(f_out)
        # f_out.write(str (no_of_nucleotides)+",\n")
        writer.writerows(out_Matrix)
    file_out_2="./Ntaxa.txt"
    with open(file_out_2, 'w') as f_out_2:
      writer = csv.writer(f_out_2)
      writer.writerow (IdList)

if __name__ == "__main__":
    main()

