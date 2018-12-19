import sys, scipy
from scipy import stats

def Calculate_Correlation(array1, array2):
    stat=stats.spearmanr(array1, array2)[0]
    return stat

def Create_Gene_Dict(genes):
    geneDict={}
    genes.next()
    for line in genes:
        line=line.rstrip().split("\t")
        geneDict[line[3]]=[float(i) for i in line[5:]]
    return geneDict

def Create_ELS_Dict(els):
    elsDictA={}
    elsDictB={}
    for line in els:
        line=line.rstrip().split("\t")
        if line[3] not in elsDictA:
            elsDictA[line[3]]=[float(line[-1]),float(line[28])]
            elsDictB[line[3]]=[float(i) for i in line[7:-1]]
        elif elsDictA[line[3]][0] < float(line[-1]):
            elsDictA[line[3]]=[float(line[-1]),float(line[28])]
            elsDictB[line[3]]=[float(i) for i in line[7:-1]]
        elif elsDictA[line[3]][1] < float(line[28]):
            elsDictA[line[3]]=[float(line[-1]),float(line[28])]
            elsDictB[line[3]]=[float(i) for i in line[7:-1]]
    return elsDictB, elsDictA

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
        statDict[line[0]]=[float(line[1]),float(line[2])]
    return statDict

genes=open(sys.argv[1])
geneDict=Create_Gene_Dict(genes)
genes.close()

symbols=open(sys.argv[2])
symbolDict=Create_Symbol_Dict(symbols)
symbols.close()

els=open(sys.argv[3])
elsDict, test =Create_ELS_Dict(els)
els.close()

stats=open(sys.argv[4])
statArray = Process_Data(stats)
stats.close()

pairs=open(sys.argv[5])

for line in pairs:
    line=line.rstrip().split("\t")
    els=line[0]
    gene=symbolDict[line[1]]
    if els in elsDict and gene in geneDict:
        cor=Calculate_Correlation(elsDict[els],geneDict[gene])
        if math.isnan(cor):
            corArray.append(0)
        else:
            corArray.append(cor)
        if statsArray[gene][1] != 0:
            Z=(cor-statsArray[gene][0])/statsArray[gene][1]
        else:
            Z=0
        p=stats.norm.sf(abs(Z))*2
        print line[2], "\t", cor, "\t", p, "\t", Z, "\t", els, "\t", line[1]
 
    else:
        print line[2], "\t", 0, "\t", 0, "\t", 0, "\t", els, "\t", line[1]
