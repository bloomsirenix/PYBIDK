import csv
def ReadCVS(filename):
    """
    Reads a CSV file and returns a list of lists.
    """
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        return list(reader)

def WriteCVS(filename, data):
    """
    Writes a list of lists to a CSV file.
    """
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(data)
