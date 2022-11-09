from dataclasses import field
from scipy import stats
import pandas as pd
import numpy as np
from EmpiricalBrownsMethod import *
import statsmodels
from joblib import Parallel, delayed
import datetime

def BrayCurtis(u,v):
    dif = 0
    sum = 0
    for i in range(0,len(u)):
        dif += abs(u[i]-v[i])
        sum += abs(u[i]+v[i])
    if sum == 0:
        return 0
    return dif/sum

def RenormalizationReBoot(u,v):
    a = list()
    b = list()
    for i,j in zip(u,v):
        sum = i+j
        if sum ==0:
            a.append(0)
            b.append(0)
        else:
            a.append(i/sum)
            b.append(j/sum)
    return a,b

def ReBoot(df):
    for col in df.columns:
        sum = 0
        for c in df[col]:
            sum += c
        for i in range(len(df[col])):
            df[col][i] = df[col][i]/sum
    return


def PermuteTaxon(t):
    #p = t.copy()
    return np.random.permutation(t)

def BootstrapTaxon(t):
    s = list()
    for i in range(len(t)):
        s.append( t[ np.random.randint(0,len(t)) ] )
    return s


def CalculateMetrics(df, reBoot = True, thresholdSp = 0.3, thresholdBC = 0.3, thresholdPval = .05):
    rows = df.to_dict('records')
    numElts = len(rows)
    spearmanCorrelation = list()
    BrayCurtisDissimilarity = list()
    network = list()
    segment = numElts
   # segment = 10000

    for i in range(segment):
#        print( f"progress:{str(i/numElts)}%" )
#        indexI = df.index[[i]][0]
        rawTaxonI = list(rows[i].values())
        for j in range(i+1,segment):
#            indexJ = df.index[[j]][0]
            rawTaxonJ = list(rows[j].values())
            taxonI = rawTaxonI
            taxonJ = rawTaxonJ
            if reBoot:
                taxonI, taxonJ = RenormalizationReBoot(rawTaxonI, rawTaxonJ)
            spearman = stats.spearmanr(taxonI,taxonJ).correlation
            brayCurtis = BrayCurtis(taxonI,taxonJ)
#            if abs(spearman) > thresholdSp and brayCurtis < thresholdBC:
            if abs(spearman) != 1 or (min(max(taxonI),max(taxonJ)) > .1):
                network.append( [i,j,spearman,brayCurtis])
    return network

def CalculateMetricsParallel(df, reboot = True, thresholdSp = 0.7, thresholdBC = 0.3, thresholdPval = .05):
    rows = df.to_dict('records')
    numElts = len(rows)
    network = list()
    segment = numElts
#    segment = 100
    print(f"total: {numElts}")
    start = datetime.datetime.now()
    for i in range(segment):
        if i%100 == 0:
            print( f"progress:{str(i)}" )
            finish = datetime.datetime.now()
            print ( f"metrics: \t {finish-start}")
        rawTaxonI = list(rows[i].values())
        copiesTaxonI = list()
        for k in range(i+1,segment):
            copiesTaxonI.append(rawTaxonI.copy())
        metrics = Parallel(n_jobs=-1)(delayed( CalculateMetricsTaxons) (copiesTaxonI[j-i-1], list(rows[j].values())) for j in range(i+1,segment) )
        for j,k in zip(range(i+1,segment),range(len(metrics))):
            spearman, brayCurtis, minMaxCheck,pval = metrics[k]
            #if abs(spearman) > thresholdSp and brayCurtis < thresholdBC:
            if ( (minMaxCheck) or abs(spearman)!= 1 ) and brayCurtis < thresholdBC and abs(spearman) > thresholdSp and pval < thresholdPval:
                network.append( [i,j,spearman,brayCurtis])
    return network


def CalculateCorrelationTaxon(df, taxon, outFileName):
    out = open(outFileName,'w')
    out.write('taxon,correlation,pvalue\n')
# parse to dictionary for faster reading
    rows = df.to_dict('records')
    i = list(df.index).index(taxon)
    numElts = len(rows)
    rawTaxonI = list(rows[i].values())
# generate copies for parallel processing
    copiesTaxonI = list()
    for k in range(numElts):
        copiesTaxonI.append(rawTaxonI.copy())
    metrics = Parallel(n_jobs=-1)(delayed( CalculateMetricsTaxons) (copiesTaxonI[j-i-1], list(rows[j].values())) for j in range(numElts) )
    for j,k in zip(range(numElts),range(len(metrics))):
        spearman, brayCurtis, minMaxCheck,pval = metrics[k]
        out.write(f"{df.index[j]},{spearman},{pval}\n")
    

def CalculateMetricsTaxons(rawTaxonI, rawTaxonJ, reBoot = False):
    if reBoot:
        taxonI, taxonJ = RenormalizationReBoot(rawTaxonI, rawTaxonJ)
    else:
        taxonI = rawTaxonI
        taxonJ = rawTaxonJ    
    regression = stats.spearmanr(taxonI,taxonJ)
    spearman = regression.correlation
    pval = regression.pvalue
    brayCurtis = BrayCurtis(taxonI,taxonJ)
    return spearman, brayCurtis, (min(max(taxonI),max(taxonJ)) > .1),pval

def PermutationTest(df, network, bootstrap =False, numPermutations = 10000, reBoot = True):
# accumulate all pvalues for later adjustment for multiple comparison
    pvalsSpearman = list()
    pvalsBrayCurtis = list()
    shallowSample = min(10,numPermutations)

    PRUEBA = False
    if PRUEBA:
        if bootstrap:
            out = open("muestrasBootstrap.csv","w")
        else:
            out = open("muestrasPermutacion.csv","w")
# obtain all taxons of interest
    taxonIds = set()
    print(f"total number of links: {len(network)}")
    for link in network:
        i = link[0]
        j = link[1]
        taxonIds.add(i)
        taxonIds.add(j)
# generate all permutation simulations
    permutations = dict()
    for i in taxonIds:
        permutations[i] = list()
        for s in range(numPermutations):
            if bootstrap:
                permutations[i].append( BootstrapTaxon( list(df.iloc[i]) ) )
            else:
                permutations[i].append( PermuteTaxon( list(df.iloc[i]) ) )
    count = 0
    start = datetime.datetime.now()
    # link format i, j, spearman, braycurtis
    for link in network:
        count+=1
        if count%1000==0:
            print(f"current link: \t{count}")
            finish = datetime.datetime.now()
            print ( f"timelapse: \t {finish-start}")
#        i,j,sp,bc = link
        i=link[0]
        j=link[1]
        sp=link[2]
        bc=link[3]

        if PRUEBA and (i != 24 and i!= 25 and i != 63):
            continue        
    # calculate metrics for all the links
        simmulationMetrics = list()
#        for s in range(numPermutations):
#            spearman, brayCurtis = CalculateMetricsTaxons( permutations[i][s], permutations[j][s]  )
#            simmulationMetrics.append([spearman, brayCurtis ])
    # test variance of spearman on shallow sample
        shallowSimmulation = Parallel(n_jobs=-1)(delayed( CalculateMetricsTaxons) (permutations[i][s], permutations[j][s]) for s in range(shallowSample) )
    # if variance is not zero continue simmulations
        if CheckShallowSimmulation(shallowSimmulation):
            # continue simmulations using shallow sample
            for simm in shallowSimmulation:
                simmulationMetrics.append(simm)
            simmulationMetrics = Parallel(n_jobs=-1)(delayed( CalculateMetricsTaxons) (permutations[i][s], permutations[j][s]) for s in range(shallowSample, numPermutations) )
        # calculate pvalues
            simmulationsSpearman = list()
            simmulationsBrayCurtis = list()
            for simm in simmulationMetrics:
                sp = simm[0]
                bc = simm[1]
                simmulationsSpearman.append(sp)
                simmulationsBrayCurtis.append(bc)
            varSpearman = np.var(simmulationsSpearman)
            meanSpearman = np.mean(simmulationsSpearman)
        # calculate pvalue depending on the side of the distribution
            if varSpearman == 0:
                pvalueSpearman = 1
            else:
                if sp > meanSpearman:
                    pvalueSpearman =  stats.norm( meanSpearman , varSpearman).cdf(sp)
                else:
                    pvalueSpearman = 1- stats.norm( meanSpearman , varSpearman).cdf(sp)
            varBrayCurtis = np.var(simmulationsBrayCurtis)
            meanBrayCurtis = np.mean(simmulationsBrayCurtis)
            if varBrayCurtis == 0:
                pvalueBrayCurtis = 1
            else:
                if bc > meanBrayCurtis:
                    pvalueBrayCurtis = stats.norm( meanBrayCurtis , varBrayCurtis).cdf(bc)
                else:
                    pvalueBrayCurtis = 1- stats.norm( meanBrayCurtis , varBrayCurtis).cdf(bc)

    # if variance is zero, there is no point to continue the simmulations and pvalue is 1
        else:
            pvalueSpearman = 1
            pvalueBrayCurtis = 1
        if PRUEBA:
            if count > 200000:
                out.close()
                return
            out.write(f"var spearman: {varSpearman} \t meanSpearman: {meanSpearman}\n")
            out.write(f"varBrayCurtis: {varBrayCurtis}\t meanBrayCurtis: {meanBrayCurtis}\n")
            out.write(f"mean taxon {i}: {np.mean(permutations[i][0])}\t var taxon {i}: {np.var(permutations[i][0])}\n")
            out.write(f"mean taxon {j}: {np.mean(permutations[j][0])}\t var taxon {j}: {np.var(permutations[j][0])}\n")
            for simm in range(len(simmulationMetrics)):
                out.write(f"taxon {i}")
                for v in permutations[i][simm]:
                    out.write(f",{v}")
                out.write(f",taxon {j}")
                for v in permutations[j][simm]:
                    out.write(f",{v}")
                normI, normJ = RenormalizationReBoot(permutations[i][simm], permutations[j][simm])
                out.write(f",norm taxon{i}")
                for v in normI:
                    out.write(f",{v}")
                out.write(f",norm taxon {j}")
                for v in normJ:
                    out.write(f",{v}")
                out.write(f",{simmulationMetrics[simm][0]},{simmulationMetrics[simm][1]}\n")
            out.write("\n")

    # link updated format i, j, spearman, braycurtis, pvalueSpearman, pvalueBrayCurtis
        link.append( pvalueSpearman)
        link.append( pvalueBrayCurtis)
        pvalsSpearman.append(pvalueSpearman)
        pvalsBrayCurtis.append(pvalueBrayCurtis)
    if PRUEBA:
        return

# adjust for multiple comparison
#    for link in network:
#        sp = link[2]
#        bc = link[3]
#        pvalsSpearman.append(sp)
#        pvalsBrayCurtis.append(bc)
    t, adjustedSpearman = statsmodels.stats.multitest.fdrcorrection(pvalsSpearman, is_sorted=False)
    t, adjustedBrayCurtis = statsmodels.stats.multitest.fdrcorrection(pvalsBrayCurtis, is_sorted=False)

# link updated format i, j, spearman, braycurtis, (pvalueSpearman, pvalueBrayCurtis, adjustedPvalSpearman, adjustedPvalBrayCurtis ), (pvalueSpearman, pvalueBrayCurtis, adjustedPvalSpearman, adjustedPvalBrayCurtis )
    for i in range(len(network)):
        network[i].append(adjustedSpearman[i])
        network[i].append(adjustedBrayCurtis[i])
# brown metric
    #   PENDING

    return  

# checks if all the simmulations have the same value
def CheckShallowSimmulation(simm):
    simmulationsSpearman = list()
    for simm in simm:
        simmulationsSpearman.append(simm[0])
    return np.var(simmulationsSpearman) != 0


def printNetwork(network, filename):
    out = open(filename,'w')
    out.write("taxon1,taxon2,Spearman Correlation,Bray Curtis Dissimilarity,pvalue Spearman,pvalue Bray Curtis,adjusted pvalue Spearman,adjusted pvalue Bray Curtis\n")
    for link in network:
        out.write(str(link[0]))
        for i in range(1,len(link)):
            out.write(","+str(link[i]))
        out.write("\n")
        
def loadNetwork(filename):
    network = list()
    file = open(filename,'r')
    file.readline()
    for line in file:
        fields = line.replace('\n','').split(',')
        link = list()
        link.append(int(fields[0]))
        link.append(int(fields[1]))

        for i in range(2,len(fields)):
            link.append(float(fields[i]))
        network.append(link)
    return network

def printNetworkGephi(network,taxonNames,filename):
    nodes = open(filename+"_nodes.csv",'w')
    edges = open(filename+"_edges.csv",'w')
    
    nodes.write("Id,Label\n")
    for i in range(len(taxonNames)):
        nodes.write(f"{i},{taxonNames[i]}\n")
    
    edges.write("Source,Target,Correlation,Dissimilarity\n")
    for link in network:
        edges.write(f"{link[0]},{link[1]},{link[2]},{link[3]}\n")

    nodes.close()
    edges.close()
    return