from operator import add
import random, subprocess,numpy, sys, sklearn, scipy, itertools, itertools, math
from sklearn import ensemble
from sklearn.ensemble import GradientBoostingClassifier


def Run_GBM(allLabels, allFeatures, valLabels, valFeatures, output,\
                      repetitions):
    output=open(output, "w")
    oob=[]
    acc=[]
    o=[0]*len(valLabels)
    FI=[0]*len(valFeatures[0])
    for x in range(0,repetitions):
	gbm = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
        gbm.fit(allFeatures,allLabels)
	predictions=gbm.predict(valFeatures)
        i=0
        correct=0
        incorrect=0
        for x in predictions:
            if x == valLabels[i]:
                correct += 1
            else:
                incorrect += 1
            i += 1
        M=gbm.predict_proba(valFeatures)
        acc.append(correct/float(correct+incorrect))
        k=0
        for entry in M:
            o[k] += entry[1]
            k += 1
        FI=[x + y for x, y in zip(FI, gbm.feature_importances_)]
    for element in FI:
        print element/float(repetitions), "\t",
    print "\n"
    i=0
    for entry in o:
        print >> output, valLabels[i], "\t", entry/float(repetitions)
        i += 1
    output.close()

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
    for line in pairs:
        line=line.rstrip().split("\t")
        pairArray.append([line[0],line[1],int(line[2])])
    return pairArray

def Process_Peak_Matrix(matrix):
    elementDict={}
    matrix.next()
    for line in matrix:
        line=line.rstrip().split("\t")
        elementDict[line[0]]=[float(i) for i in line[1:]]
    return elementDict

def Create_Feature_Array(data, enhancerSignals, tssSignals, geneDict):
    labels=[]
    features=[]
    for pair in data:
        tssFeatures=[]
        for tss in geneDict[pair[1]]:
            if len(tssFeatures) > 0:
                tssFeatures=map(add, tssFeatures, tssSignals[tss])
            else:
                tssFeatures=tssSignals[tss]
        tssFeatures=[x / float(len(geneDict[pair[1]])) for x in tssFeatures]
        features.append(enhancerSignals[pair[0]]+tssFeatures)
        labels.append(pair[2])
    return features, labels

trainingPairs=open(sys.argv[1])
trainingArray=Process_ELS_Gene_Pairs(trainingPairs)
trainingPairs.close()

validationPairs=open(sys.argv[2])
validationArray=Process_ELS_Gene_Pairs(validationPairs)
validationPairs.close()

enhancerMatrix=open(sys.argv[3])
enhancerSignals=Process_Peak_Matrix(enhancerMatrix)
enhancerMatrix.close()

tssMatrix=open(sys.argv[4])
tssSignals=Process_Peak_Matrix(tssMatrix)
tssMatrix.close()

tss=open(sys.argv[5])
tssDict, geneDict=Create_Gene_Dict(tss)
tss.close()


trainFeat, trainLab=Create_Feature_Array(trainingArray, enhancerSignals, \
                                         tssSignals, geneDict)
valFeat, valLab=Create_Feature_Array(validationArray, enhancerSignals, \
                                         tssSignals, geneDict)

output=sys.argv[6]+"-Output.txt"
Run_GBM(trainLab, trainFeat, valLab, valFeat, output, 25)

