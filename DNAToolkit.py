#DNA Toolkit file
"""

    This file contains the functions for the DNA?RNA Toolkit.
    Copyriright (C) 2021 manikineko.nl Under the MIT Licsensing.
"""
import  collections
import json

Nucleotides = ['A', 'C', 'G', 'T']
class DNA:
    """
    DNA Class
    """
    def __init__(self, dna_seq):
        """
        Initializes the DNA class.
        """
        self.seq = validateSeq(dna_seq)
        if self.seq == False:
            raise ValueError('Invalid sequence')
        self.freq = countNucFreq(self.seq)
        self.rna = transcribe(self.seq)
        self.revcomp = reverseComplement(self.seq)
        self.orf = findORF(self.seq)
        self.orf_all = findORF_all(self.seq)
        self.orf_all_sorted = findORF_all_sorted(self.seq)
        self.orf_all_sorted_longest = findORF_all_sorted_longest(self.seq)
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)
def validateSeq(dna_seq):
    """
    Validates the sequence to ensure that it contains only valid nucleotides.
    """
    tmpseq = dna_seq.upper();
    for nucleotide in tmpseq:
        if nucleotide not in Nucleotides:
            return False
    return tmpseq
def countNucFreq(dna_seq):
    """
    Counts the frequency of each nucleotide in the sequence.
    """
    tmpseq = dna_seq.upper();
    freq = {'A':0, 'C':0, 'G':0, 'T':0}
    for nucleotide in tmpseq:
        freq[nucleotide] += 1
    return freq
def transcribe(dna_seq):
    """
    Transcribes the DNA sequence into RNA.
    """
    tmpseq = dna_seq.upper();
    rna_seq = ''
    for nucleotide in tmpseq:
        if nucleotide == 'T':
            rna_seq += 'U'
        else:
            rna_seq += nucleotide
    return rna_seq
def reverseComplement(dna_seq):
    """
    Returns the reverse complement of the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    revcomp = ''
    for nucleotide in tmpseq:
        if nucleotide == 'A':
            revcomp += 'T'
        elif nucleotide == 'T':
            revcomp += 'A'
        elif nucleotide == 'C':
            revcomp += 'G'
        elif nucleotide == 'G':
            revcomp += 'C'
    return revcomp[::-1]
def findORF(dna_seq):
    """
    Finds the longest open reading frame in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf = ''
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            break
    for i in range(0, len(orf), 3):
        codon = orf[i:i+3]
        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
            return orf[:i]
    return orf
def findORF_all(dna_seq):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return orf_all
def findORF_all_sorted(dna_seq):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)
def findORF_all_sorted_longest(dna_seq):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[0]
def findORF_all_sorted_longest_n(dna_seq, n):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[:n]
def findORF_all_sorted_longest_n_r(dna_seq, n):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[-n:]
def findORF_all_sorted_longest_n_r_r(dna_seq, n):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[-n:][::-1]
def findORF_all_sorted_longest_n_r_r_r(dna_seq, n):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[-n:][::-1][::-1]
def findORF_all_sorted_longest_n_r_r_r_r(dna_seq, n):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[-n:][::-1][::-1][::-1]
def findORF_all_sorted_longest_n_r_r_r_r_r(dna_seq, n):
    """
    Finds all the longest open reading frames in the DNA sequence.
    """
    tmpseq = dna_seq.upper();
    orf_all = []
    for i in range(0, len(tmpseq), 3):
        codon = tmpseq[i:i+3]
        if codon == 'ATG':
            orf = tmpseq[i:]
            for j in range(0, len(orf), 3):
                codon = orf[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orf_all.append(orf[:j])
                    break
    return sorted(orf_all, key=len, reverse=True)[-n:][::-1][::-1][::-1][::-1]

