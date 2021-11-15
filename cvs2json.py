from os import read
from DNASeq import *
from cvstoolkit import *
def main():
    """
    Main function.
    """
    filepath = input("Enter the the DNA CVS File Path: ");
    # Read the csv file and convert it to a python object.
    DNAStr = readDNACVS(filepath);
    # Convert the String to a DNA object.
    DNAOBJ =  DNASeq2DNAOBJ(DNAStr);
    #Convert DNAOBJ to a json file.
    file = open(filepath + ".dna.json", "w");
    file.write(DNAOBJ.toJSON());

main()