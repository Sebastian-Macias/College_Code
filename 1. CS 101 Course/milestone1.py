#!/usr/bin/env python
# coding: utf-8

# In[ ]:


genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
# Many codes rely on the genetic_code variable
# So its much easier to define it here once

def locate_substring(dna_snippet, dna):
    """ The idea with this is to find the locations of occurences of a
        snippet within a given dna string. Not count the occurences. Just
        find where they are hiding at.
    """
    index = 0
    occurence_list = []
    while dna.find(dna_snippet, index) != -1:
        occurence_list.append(dna.find(dna_snippet, index))
        index = 1 + dna.find(dna_snippet, index)
    return occurence_list

def hamming_dist(dna1, dna2):
    """ The goal with this code is to figure out the hamming distance between
        any two given input strings.
    """
    i=0
    dist = 0
    while i != len(dna1): #This function proceeds until all values are checked
        if dna2.find(dna1[i:i+1:],i,i+1) < 0:
            # This if statement will take a given value in the first string,
            # and search for it in the same spot in the second string. If
            # it is found, then it is added to the count. 
            dist = dist + 1
        i = i + 1 # i represents the place in line.
    return dist

def mendels_law(Hom, Het, Rec):

    total = Hom + Het + Rec;
    probability = (4 *( Hom * (Hom - 1) + 2 * Hom * Het + 2 * 
                       Hom * Rec + Het * Rec) + 3 * Het*(Het - 1)) / (
                4 * total * (total - 1));
    return probability

def fibonacci_rabbits(n, k):
    f1, f2 = 1, 1
    for i in range(n - 1):
        f2, f1 = f1, f1 + (f2 * k)
    return f2

def dna_count(dna):
    # take str representing DNA sequence
    #return integers representing number of times each nucleotide occurs
    dna = dna.upper() #make uppercase
    count_A = dna.count('A') #count As
    count_C = dna.count('C') #count Cs
    count_G = dna.count('G') #count Gs
    count_T = dna.count('T') #count Ts
    return {'A':count_A, 'C': count_C, 'G': count_G, 'T':count_T} #return dict

def dna2rna(dna):
    #replace all occurrences of "T" with "U"
    #dna transcribed to rna
    rna = ''
    for symbol in dna:
        if symbol == 'T':
            rna = rna + 'U' #replace all Ts with Us
        elif symbol == 'A':
            rna = rna + 'A' #keep As the same
        elif symbol == 'C':
            rna = rna + 'C' #keep Cs the same
        else:
            rna = rna + 'G' #keep Gs the same
    return rna

def reverse_complement(dna):
    #reverse entire dna string and then take complement of each nucleotide
    str = dna
    for symbol in str: #take the complement of each nucleotide-- A<->T, C<->G
        if symbol == 'A':
            str = str + 'T'
        elif symbol == 'C':
            str = str + 'G'
        elif symbol == 'G':
            str = str + 'C'
        else:
            str = str + 'A'
    stringlen = len(dna) -1 #take whole string
    reverse = str[:stringlen:-1] #make string be reversed
    return reverse

def count_dom_phenotype(genotypes): #define func
    genotypes = str(genotypes).replace(",", "") #cut commas
    genotypes = genotypes.replace(" ", "") #cut spaces
    ucounter = 0 #initialize uppercount #initialize
    lcounter = 0 #initialize lowercount
    for i in genotypes: #for letter in str
        ucounter = genotypes.count("1") 
        lcounter = genotypes.count("0")
        length = (ucounter+lcounter)
    return ucounter + lcounter/length

def source_rna(protein):
    """The logic of the program
    protein is 'MA'
    M  = 1, A = 4 == 4*1*3 = 12
    """
    val_list = list(genetic_code.values()) 
    #want to create a list with the genetic_code 
    #and takes the values post :
    seq = 1 #want to initialize at 1
    for letter in protein: #look for letter in protein
        count = val_list.count(letter) #count letter in protein
        seq = seq*count 
        #sequence*count multiplies number of 
        #occurances of letter with eachother
    return seq*3 #multiplies by protein length

def rna2codon(triplet):
    output = ""
    #establishing triplet length and range
    for i in range(0,len(triplet)-3,3):
        #outputing codon value for any given rna triplet
        output += genetic_code[triplet[i:i+3]]
    return output

def GC_content(dna_list):
    for dna in dna_list:
        b = 0
        l = []
        #finding the index with the highest GC concentration
        if (dna.count('C')+dna.count('G'))/len(dna) > b:
            #documenting index with highest rate
            a = dna_list.index(dna)
            #documenting GC content of string with maximum
            b = (dna.count('C')+dna.count('G'))/len(dna)
            #representing index and rate of GC content as percent
            l.extend([a,b*100])
    #returning max index and value in tuple form (a,b)
    return tuple(l)

def rna2codons(triplets): 
    # takes triplet form of RNA 
    # and runs it through rna2codon to find the proteins
    # Not a required function, but necessary for other functions
    proteins = ''
    for i in range( 0,int( len(triplets) / 3 ) ):
        proteins = proteins + rna2codon(triplets[ 3*i:3*i+3 ])
    return proteins

def splice_rna(dna,intron_list):
    #rna = dna2rna(dna)  # Convert DNA to RNA
    intronstr0 = ''.join(intron_list[0]) #take the introns in list and separate (want to make into str)
    intronstr1 = ''.join(intron_list[1])
    for i in dna:
        if intronstr0 in dna: #check if intron string in DNA
            dna = dna.replace(intronstr0, "") #replace 
        if intronstr1 in dna: #check if intron 2 string in DNA
            dna = dna.replace(intronstr1, "") #replace
    rna = dna2rna(dna) # Convert DNA to RNA
    return rna2codons(rna).replace("*","")

