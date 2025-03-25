import os
import sys
sys.path.append(r"D:\code\python\LHA\pyswcloader-main")
import pyswcloader
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pandas import ExcelWriter
from statistics import mean
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster,set_link_color_palette
import seaborn as sns
matplotlib.use('TkAgg')
from sklearn.cluster import AgglomerativeClustering


scorepath = r"D:\data\altas\morphologyfinal\LHA_dist.csv"
scoredata = pd.read_csv(scorepath)
indexlist = scoredata['Unnamed: 0'].tolist()
mydata = scoredata.drop(columns='Unnamed: 0')
mydata.index = indexlist
scores=mydata.copy(deep=True)

# clustering
""" info1 = pyswcloader.cluster.cluster(n_cluster=24,method='hierarchy',
                                   feature='precomputed',matrix=scores) """

# recompute
vec = squareform(scores,checks=False)
Z = linkage(vec,'ward')
f = fcluster(Z,t=24,criterion='maxclust')      # 24 subtypes


info = pd.DataFrame()
info['label']=f
info['filepath'] = list(scores.index)
info['neuron'] = [item for item in info.filepath]

fig1 = plt.figure(dpi=300,figsize=(20,4))
ax=plt.gca()
# set_link_color_palette(['#FF0000','#228B22','#4682B4','#A020F0'])
dendrogram(Z,color_threshold=7000,truncate_mode='none')
plt.axhline(y=7000,color='gray',linestyle='--')
ax.axes.xaxis.set_visible(True)
ax.axes.yaxis.set_visible(True)
plt.tick_params(axis='x',labelsize=0.5)
plt.show()

with ExcelWriter(r"D:\data\altas\morphologyfinal\info1112.xlsx") as writer1:
    info.to_excel(writer1)

order = pd.read_excel(r"D:\data\altas\morphologyfinal\order.xlsx",header=None)
orderlist = order.iloc[:,0].tolist()

index_resort = []
for i in orderlist:
    index_resort.append(indexlist[i])

neuronorder = pd.DataFrame(index_resort)
with ExcelWriter(r"D:\data\altas\morphologyfinal\neuronorder.xlsx") as writer2:
    neuronorder.to_excel(writer2)
