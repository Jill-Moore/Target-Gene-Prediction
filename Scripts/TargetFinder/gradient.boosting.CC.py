from operator import add
import random, subprocess,numpy, sys, sklearn, scipy, itertools, shutil, math, joblib
from sklearn import ensemble
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier

        
def Run_GBM(valLabels, valFeatures, outputPrefix1, outputPrefix2, repetitions, \
    i, header, outputDir1, outputDir2, version):
    
    valOutput=outputPrefix2+"-validation-"+str(i)+"."+version+"."+outputPrefix1+".txt"
    output=open(valOutput, "w")
    o=[0]*len(valLabels)
    for x in range(0,repetitions):
        gbm=joblib.load(outputDir1+'/'+outputPrefix1+'.4C.model.'+str(x)+'.'+str(i)+'.pkl')
        M=gbm.predict_proba(valFeatures)
        k=0
        for entry in M:
            o[k] += entry[1]
            k += 1
    i=0
    for entry in o:
        print >> output, valLabels[i], "\t", entry/float(repetitions)
        i += 1
    output.close()
    shutil.move(valOutput, outputDir2+"/"+valOutput)


def Run_Model(trainFeat,trainLab, outputPrefix1, outputPrefix2, cvList, \
    cvGroups, header, outputDir1, outputDir2, version):
    i=1
    for group in cvList:
        print "Running cross-validation for "+ group
        testMatrix=[]
        testY=[]
        j=0
        for entry in cvGroups:
            if entry == group:
                testMatrix.append(trainFeat[j])
                testY.append(trainLab[j])
            j+=1
        Run_GBM(testY, testMatrix, outputPrefix1, outputPrefix2, 1, group, \
            header, outputDir1, outputDir2, version)
        i+=1    
    return

def Create_Distance_Dict(distances):
    distanceDict={}
    for line in distances:
        line=line.rstrip().split("\t")
        distanceDict[line[0].rstrip()+"-"+line[1].rstrip()]=[int(line[2])]
    return distanceDict

def Create_Gene_Dict(tss):
    tssDict={}
    geneDict={}
    for line in tss:
        line=line.rstrip().split("\t")
        tssDict[line[3]]=line[4]
        if line[4] not in geneDict:
            geneDict[line[4]]=[line[3]]
        else:
            geneDict[line[4]].append(line[3])
    return tssDict, geneDict
    
def Process_ELS_Gene_Pairs(pairs):
    pairArray=[]
    cvList=[]
    for line in pairs:
        line=line.rstrip().split("\t")
        pairArray.append([line[0],line[1],int(line[2]),line[3]])
        if line[3] not in cvList:
            cvList.append(line[3])
    return pairArray, cvList

def Process_Peak_Matrix(matrix, mode, header):
    elementDict={}
    h=matrix.next().rstrip().split("\t")[1:]
    for entry in h:
        header.append(entry+"-"+mode)
    for line in matrix:
        line=line.rstrip().split("\t")
        elementDict[line[0]]=[float(i) for i in line[1:]]
    return elementDict, header

def Create_Feature_Array(data, enhancerSignals, tssSignals, optionalSignals, geneDict):
    labels=[]
    cvGroups=[]
    features=[]
    for pair in data:
        tssFeatures=[]
        for tss in geneDict[pair[1]]:
            if len(tssFeatures) > 0:
                tssFeatures=map(add, tssFeatures, tssSignals[tss])
            else:
                tssFeatures=tssSignals[tss]
        tssFeatures=[x / float(len(geneDict[pair[1]])) for x in tssFeatures]
        optionalFeatures=[]
        for i in optionalSignals:
            optionalFeatures+=i[pair[0]+"-"+pair[1]]
        features.append(enhancerSignals[pair[0]]+tssFeatures+ \
            optionalFeatures)
        labels.append(pair[2])
        cvGroups.append(pair[3])
    return features, labels, cvGroups

header=[]
trainingPairs=open(sys.argv[1])
trainingArray, cvList=Process_ELS_Gene_Pairs(trainingPairs)
trainingPairs.close()

enhancerMatrix=open(sys.argv[2])
enhancerSignals, header=Process_Peak_Matrix(enhancerMatrix, "enhancer",header)
enhancerMatrix.close()

tssMatrix=open(sys.argv[3])
tssSignals, header=Process_Peak_Matrix(tssMatrix,"tss",header)
tssMatrix.close()

optionalSignals=[]
feature1=sys.argv[4]
if "Window-Feature-Matrix" in feature1:
    windowMatrix=open(sys.argv[4])
    windowSignals, header=Process_Peak_Matrix(windowMatrix,"window",header)
    windowMatrix.close()
    optionalSignals.append(windowSignals)
elif "Distance" in feature1:
    distances=open(sys.argv[4])
    distanceDict=Create_Distance_Dict(distances)
    header.append("distance")
    distances.close()
    optionalSignals.append(distanceDict)
else:
    i=4

feature2=sys.argv[5]
if "Distance" in feature2:
    distances=open(sys.argv[5])
    distanceDict=Create_Distance_Dict(distances)
    header.append("distance")
    distances.close()
    optionalSignals.append(distanceDict)
    i=6
else:
    i=5
print i
tss=open(sys.argv[i])
tssDict, geneDict=Create_Gene_Dict(tss)
tss.close()

outputPrefix1=sys.argv[i+1]
outputPrefix2=sys.argv[i+2]

outputDir1=sys.argv[i+3]
outputDir2=sys.argv[i+4]

version=sys.argv[i+5]

print "processing feature matrices..."

trainFeat, trainLab, cvGroups, = Create_Feature_Array(trainingArray, \
    enhancerSignals, tssSignals, optionalSignals, geneDict)

print "running models..."

print outputPrefix1
print outputPrefix2
print outputDir1
print outputDir2


Run_Model(trainFeat, trainLab, outputPrefix1, outputPrefix2, cvList, cvGroups, \
    header, outputDir1, outputDir2, version)
