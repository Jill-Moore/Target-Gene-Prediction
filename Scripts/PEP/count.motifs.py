import sys, numpy

def Create_Region_Dict(regions, length):
    regionDict={}
    for line in regions:
        line=line.rstrip().split("\t")
        regionDict[line[3]]=[numpy.zeros(length),int(line[2])-int(line[1])+1]
    return regionDict

def Create_Motif_Array(motifs):
    motifArray=[]
    for line in motifs:
        line=line.rstrip()
        motifArray.append(line)
    return motifArray

def Process_FIMO_Output(fimo, motifArray, regionDict):
    fimo.next()
    for line in fimo:
        line=line.rstrip().split("\t")
        if float(line[7]) < 1E-4:
            regionDict[line[2]][0][motifArray.index(line[0])] += 1
    return regionDict
    

    
motifs=open(sys.argv[1])
motifArray=Create_Motif_Array(motifs)
motifs.close()

regions=open(sys.argv[2])
regionDict=Create_Region_Dict(regions, len(motifArray))
regions.close()

fimo=open(sys.argv[3])
regionDict=Process_FIMO_Output(fimo, motifArray, regionDict)
fimo.close()

print "cREs"+"\t"+"\t".join(motifArray)

for region in regionDict:
    print region+"\t"+ "\t".join([str(i) for i in \
    numpy.ndarray.tolist(regionDict[region][0]/regionDict[region][1])])

