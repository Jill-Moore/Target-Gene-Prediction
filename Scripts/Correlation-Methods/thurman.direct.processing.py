import sys, scipy
from scipy import stats


def Create_Symbol_Dict(symbols):
    symbolDict={}
    for line in symbols:
        line=line.rstrip().split("\t")
        symbolDict[line[3]]=line[6]
    return symbolDict

def Create_Stat_Dict(stats):
    statDict={}
    for line in stats:
        line=line.rstrip().split("\t")
        if line[0] not in statDict:
            statDict[line[0]]={line[1]:float(line[2])}
        elif line[1] not in statDict[line[0]]:
            statDict[line[0]][line[1]]=float(line[2])
        elif float(line[2]) > statDict[line[0]][line[1]]:
            statDict[line[0]][line[1]]=float(line[2]) 
    return statDict

symbols=open(sys.argv[1])
symbolDict=Create_Symbol_Dict(symbols)
symbols.close()

stats=open(sys.argv[2])
statDict = Create_Stat_Dict(stats)
stats.close()

pairs=open(sys.argv[3])

for line in pairs:
    line=line.rstrip().split("\t")
    gene=symbolDict[line[1]]
    els=line[0]
    if els in statDict and gene in statDict[els]:
        cor=statDict[els][gene]
    else:
        cor=-1
    print line[2], "\t", cor, "\t", els, "\t", line[1]
