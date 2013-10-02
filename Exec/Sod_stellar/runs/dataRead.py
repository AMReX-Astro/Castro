import string
import numpy

def getFileParams(filename):

    # get the number of lines and fields in the data file                                                             
    # skip comment and blank lines, and also, require that                                                            
    # the time (first field) is monotonically increasing                                                              

    mf = open(filename, "r")

    numLines = 0
    numFields = -1

    oldTime = -1.0

    for line in mf:
        if (not (line.startswith("#") or line.startswith(" #") or line.lstrip() == "")):

            fields = string.split(line)
            time = float(fields[0])

            if (numFields == -1):
                numFields = len(fields)
                oldTime = time
                numLines += 1

            elif (time > oldTime):
                oldTime = time
                numLines +=1

    mf.close()

    return numLines, numFields


def getData(filename):

    # get the data from the data file                                                                                 
    # skip comment and blank lines, and also, require that                                                            
    # the time (first field) is monotonically increasing                                                              

    numLines, numFields = getFileParams(filename)

    data = numpy.zeros( (numLines, numFields), numpy.float64)

    mf = open(filename, "r")

    indexLine = 0

    oldTime = -1.0

    for line in mf:
        if (not (line.startswith("#") or line.startswith(" #") or line.lstrip() == "")):

            fields = string.split(line)
            time = float(fields[0])

            if (oldTime == -1.0):
                data[indexLine,:] = fields
                oldTime = time
                indexLine += 1

            elif (time > oldTime):
                data[indexLine,:] = fields
                oldTime = time
                indexLine +=1

    mf.close()

    return data
