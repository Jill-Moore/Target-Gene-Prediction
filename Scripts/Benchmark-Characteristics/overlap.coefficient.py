import sys
from collections import Counter

def Overlap_Coefficient(d1, d2):
    intersection=len(list(set(d1).intersection(d2)))
    minimum=min(len(d1),len(d2))
    if minimum == 0:
        return 0
    else:
        return intersection/float(minimum)
    
def Create_Data_Array(data):
    dataArray=[]
    data=open(data)
    for line in data:
        line=line.rstrip().split("\t")
        if line[2] == "1 ":
            dataArray.append(line[0]+"-"+line[1])
    data.close()
    return dataArray

d1Array=Create_Data_Array(sys.argv[1])
d2Array=Create_Data_Array(sys.argv[2])

print Overlap_Coefficient(d1Array, d2Array)
