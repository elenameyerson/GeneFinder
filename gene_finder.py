# -*- coding: utf-8 -*-
"""
My first software design mini-project.

@author: Elena Meyerson

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

from load import load_seq
dna = load_seq("./data/X73525.fa")




def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        I added the third test to make sure that the error message worked correctly.

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('D')
    'Error, Invalid Input'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    else:
        return 'Error, Invalid Input'



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        I added the third test to make sure that the function had the correct error message when an invalid string of DNA was entered.

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CCGCGTTCAD")
    'Error, Invalid Input'
    """
    comp = ''
    for letter in dna:
        if get_complement(letter) == 'Error, Invalid Input':
            return('Error, Invalid Input')
        else:
            comp = get_complement(letter) + comp
    return comp


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        I added the 3rd test to make sure that the function works if there is no stop codon.
        I added the 4th test to make sure that the function works if the stop codon makes up the last 3 letters.

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAAGAT")
    'ATGAGA'
    >>> rest_of_ORF("ATGTGGGAAAGAGAA")
    'ATGTGGGAAAGAGAA'
    >>> rest_of_ORF("ATGTGGGAAAGAGAATAG")
    'ATGTGGGAAAGAGAA'
    """
    num_letters = len(dna)
    x = 0
    protein = ''

    while x+2 < num_letters:
        if str(dna[x]+dna[x+1]+dna[x+2]) == 'TAG':
            return protein
            break
        elif str(dna[x]+dna[x+1]+dna[x+2]) == 'TAA':
            return protein
            break
        elif str(dna[x]+dna[x+1]+dna[x+2]) == 'TGA':
            return protein
            break
        else:
            protein = protein + dna[x]+dna[x+1]+dna[x+2]
            x = x+3
    return dna


def find_all_ORFs_oneframe(dna,start):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. thfrom load import load_seq
ey start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I added the second test to make sure that the code outputs nothing if there are no start codons in the string.
        I added the third test to make sure that the code doesnt ouput nested ORFs as well.

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC",0)
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe('ATAGCTAGCTAGCTGAACTTAG',0)
    []
    >>> find_all_ORFs_oneframe("ATGATGGAATGTAGATAGATGTGCCC",0)
    ['ATGATGGAATGTAGA', 'ATGTGCCC']
    """

    num_letters = len(dna)
    x = start
    ORF_list = []

    while x+2 < num_letters:
        if str(dna[x]+dna[x+1]+dna[x+2]) == 'ATG':
            non_nested_orf = rest_of_ORF(dna[x:])
            next_start = len(non_nested_orf)
            ORF_list = ORF_list + [non_nested_orf]
            x = x + next_start
        else:
            x = x+3
    return ORF_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    I added the second test to make sure that everything works if there are no ORFs.

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATCAGAATTAG")
    []
    """

    full_ORF_list = []

    for i in range(0,3):
        full_ORF_list = full_ORF_list + find_all_ORFs_oneframe(dna, i)
    return full_ORF_list




def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    all_ORFs = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return all_ORFs



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    list_of_ORFs = find_all_ORFs_both_strands(dna)
    current_longest = []
    for sequence in list_of_ORFs:
        if len(sequence)>=len(current_longest):
            current_longest = sequence
        else:
            current_longest = current_longest
    return current_longest



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_length = 0
    for i in range(num_trials+1):
        new_dna = shuffle_string(dna)
        #print(new_dna)
        longest_string = longest_ORF(new_dna)
        #print(longest_string)
        if len(longest_string) >= longest_length:
            longest_length = len(longest_string)
        #print(longest_length)
    return longest_length


#longest_ORF_noncoding('ATGACTTACACCCCACCTCTAAAATAGCCATAACTAGAATGACAGACTAGCATGGACTACTAG',500)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    x=0
    protein = ''
    while x +2 < len(dna):
        codon = dna[x]+dna[x+1]+dna[x+2]
        amino_acid = aa_table[codon]
        protein = protein + amino_acid
        x = x+3
    return protein


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    #print(threshold)
    long_enough = []
    protein_list = []

    total_orfs = find_all_ORFs_both_strands(dna)
    for orf in total_orfs:
        #print(orf)
        if len(orf)>(int(threshold)):
            #print(len(orf))
            long_enough = long_enough + [orf]
            #print(long_enough)
        #else:
            #print('Not long enough')
        #print('This list is long enough:', long_enough, 'End of List!')
    for orfs in long_enough:
        #print(orfs)
        protein_list = protein_list + [coding_strand_to_AA(orfs)]

    print(protein_list)
    return protein_list


gene_finder(dna)

if __name__ == "__main__":
    import doctest
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
    #doctest.testmod()
