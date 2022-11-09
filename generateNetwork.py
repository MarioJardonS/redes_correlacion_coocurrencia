import math
from biom import load_table
from scipy.spatial import distance
from scipy import stats
import pandas as pd
import numpy as np
import datetime
from joblib import Parallel, delayed
import statisticsFunctions
import os
import sys

if not os.path.exists('networks'):
    os.mkdir('networks')

np.random.seed(1)

start = datetime.datetime.now()
numPermutations = 100
numBootstraps   = 100
loadNetworkFlag = False

#level = 1

#for level in range(2,32):

#table = load_table('D:/data/genomics/unam/agromicrobioma/pacific_ocean.biom')
#table = load_table('/mnt/d/data/genomics/unam/agromicrobioma/pacific_ocean.biom')
table = load_table(f"maiz.txt_9.biom")

#table = load_table('~/solena/biom/solena_maiz.biom')
outName = f"maiz_species"

# https://biom-format.org/documentation/table_objects.html
numTaxons = int(table.shape[0])
numSamples = int(table.shape[1])

rawData = table.to_dataframe().sparse.to_dense()
statisticsFunctions.ReBoot(rawData)
finish = datetime.datetime.now()
print ( f"loading: \t {finish-start}")

network =list()

if loadNetworkFlag:
    network = statisticsFunctions.loadNetwork('maiz_network_Pval.csv')
else:
        #network = efficientStats.CalculateMetrics(rawData)
    network = statisticsFunctions.CalculateMetricsParallel(rawData)
    statisticsFunctions.printNetwork(network,f"networks/{outName}_raw_network.csv")

    statisticsFunctions.printNetworkGephi(network,list(rawData.index),f"networks/{outName}_network")

sys.exit()

finish = datetime.datetime.now()
print ( f"raw network: \t {finish-start}")


statisticsFunctions.PermutationTest(rawData, network, numPermutations = numPermutations, reBoot = True)
finish = datetime.datetime.now()
statisticsFunctions.printNetwork(network,"maiz_network_PermTest.csv")
print ( f"PERMUTATION test: \t {finish-start}")


statisticsFunctions.PermutationTest(rawData, network, bootstrap=True, numPermutations = numPermutations, reBoot = True)
finish = datetime.datetime.now()
statisticsFunctions.printNetwork(network,"maiz_network_complete.csv")
print ( f"BOOTSTRAP test: \t {finish-start}")
