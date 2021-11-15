# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""
    pass


class ValueTooSmallError(Error):
    """Raised when RNA Parser encountered a RNA Sequence Error"""
    pass


class DNASyntaxError(Error):
    """Raised when the DNA Parser has encoutered a DNA Syntax Error"""
    pass


#DNA Sequencer:
"""
    This is a DNA sequencer. 
    It converts a CVS file to python objects.
    It can also convert a python object to a csv file.
    Copyriright (C) 2021 manikineko.nl Under the MIT Licsensing.
"""
from DNAToolkit import *
from cvstoolkit import *
import random


rawDNAStr = "";


def readDNACVS(filepath):
    """
    Reads a csv file and converts it to a python object.
    """
    global rawDNAStr;
    for csvchunk in ReadCVS(filepath):
      for DoubleCharSeq in csvchunk:
        if validateSeq(DoubleCharSeq):
          rawDNAStr += DoubleCharSeq + " ";
    DNAStr = validateSeq(rawDNAStr.replace(" ",""));
    if(DNAStr == False):
        raise DNASyntaxError("Invalid DNA sequence");
    else:
        return DNAStr
def DNASeq2DNAOBJ(DNASeq):
    """
    Converts a DNA sequence to a DNA object.
    """
    DNAObj = DNA(DNASeq);
    DNAObj.seq = DNASeq;
    return DNAObj;
    
