import sys
sys.path.append("D:\\code\\python\\LHA\\neuron-vis-master\\neuronVis")
import BrainRegion as BR
import IONData
import copy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pandas import ExcelWriter
import pickle
import matplotlib
matplotlib.use('TkAgg')
# %matplotlib inline
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, set_link_color_palette
from sklearn.metrics import silhouette_samples,silhouette_score
from scipy.cluster.hierarchy import fcluster
import matplotlib.colors as mcl
import Scene
import json
from jsonpath import jsonpath
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import pickle


cluster_path = r"D:\data\altas\morphologyfinal\neuronorder.xlsx"
cluster_data = pd.read_excel(cluster_path)
info_ip = '10.10.48.110'
iondata = IONData.IONData()


lhasamplelist = cluster_data['neuron'].tolist()
sclist = []
swcpath = "D:\\data\\altas\\swc"
for i in range(0,len(lhasamplelist)):
    name1,name2 = lhasamplelist[i].split('_')
    sampleid = name1
    neuronname = name2+'.swc'
    neuronpath = swcpath+"\\"+sampleid+"\\"+neuronname
    axondata = pd.read_csv(neuronpath,sep=' ',header=None)
    if axondata.shape[1]==7:
        axonlength = axondata.dropna(axis=0,how='any')
        xdata = axonlength.iloc[:,2]
        max_x = max(xdata)
        if max_x>=13200:
            is_sc = 1
        else:
            is_sc = 0
    else:
        is_sc = 0
    sclist.append(is_sc)



br = BR.BrainRegion()
br.praseJson()
def get_children_region(regionlist):
    regionchildrenlist=[]
    for region in regionlist:
        subbr= br.getRegion(region)
        for i in subbr.children:
            regionchildrenlist.append(i.name)
    return regionchildrenlist

def get_all_subregion(region):
    list=[]
    regionlist=br.getRegionList(region)
    for i in regionlist[0:]:
        list.append(i[0])
    return list


with open(r"D:\code\python\LHA\neuron-vis-master\resource\1.json") as f:
    myfile = json.load(f)

structuredata = myfile['msg'][0]['children']
bigregion = [i['acronym'] for i in structuredata]
# grey, fiber tracts, VS, grv, retina
greydata = structuredata[0]
result = jsonpath(greydata,'$..children[?(@.st_level==8)]')
for i in range(0,len(result)):
    del result[i]['children']

region8 = [x['acronym'] for x in result]



lhasamplelist = cluster_data['neuron'].tolist()
lhanamelist = []
lhaprojection = []
for i in range(0,len(lhasamplelist)):
    name1,name2 = lhasamplelist[i].split('_')
    sampleid = name1
    name = name2+'.swc'
    lhanamelist.append(lhasamplelist[i])
    neuronproperty = iondata.getNeuronPropertyByID(sampleid,name)
    projectlength = neuronproperty['projectregion']
    lhaprojection.append(projectlength)     # 916*663

info = pd.DataFrame(lhaprojection)
info.index=lhanamelist
info=info.fillna(0)
info1=info.copy(deep=True)


for i in region8:
    if i not in info1.columns.tolist():
        sub = get_children_region([i])
        if len(sub)!=0:
            info1[i]=0
            for j in get_children_region([i]):
                if j in info1.columns.tolist():
                    info1[i]+=info1[j]
                    info1.pop(j)
# 916*539
info1[info1<500]=0
info1.loc['sum']=info1.apply(lambda x:x.sum())
for i in info1.columns:
    if info1[i]['sum']==0:
        info1.pop(i)
# 917*443
info1.drop('sum')
info2=(np.log2(info1[0:916]/1000.0+1))
# 916*443
greysubregionlist=get_all_subregion('grey')[1:]
greysubregionlist = [x for x in greysubregionlist if x != 'unknow']
greysubregionlist = ['CUL4,5' if x == 'CUL4, 5' else x for x in greysubregionlist]
greysubregionlist.remove('HY')
greysubregionlist.remove('LHA')

info3=info2.copy(deep=True)
existing_columns = [col for col in greysubregionlist if col in info3.columns]
info3 = info3[existing_columns]    # 916*369
info4 = info3.T

fig1=plt.figure(figsize=(8,5),dpi=300)
ax=plt.gca()
g=sns.heatmap(info4,cmap='magma',yticklabels=True,xticklabels=True,
            cbar_kws={'shrink':0.5})
plt.xticks(fontsize=1)
plt.yticks(fontsize=1)
plt.show()


# spinal cord
mydata=info1.copy(deep=True)
df = pd.read_excel("D:\\data\\altas\\brainlist.xlsx")
regionlist = df['region'].to_list()
grouplist1 = df['group'].to_list()
all_columns = mydata.columns.to_list()
dividing = [x for x in regionlist if x not in all_columns]


for i in dividing:
    subregions = get_all_subregion(i)
    subregions = [x for x in subregions if x != 'unknow']
    newlist=[0]*mydata.shape[0]
    for j in subregions:
        if j in all_columns:
            divided_length = mydata[j].to_list()
            newlist=list(map(lambda x,y:x+y,newlist,divided_length))
    mydata[i]=newlist 
# 917*447
mydata = mydata[regionlist]
mydata1 = (np.log2(mydata[0:916]/3000.0+1))    # 916*46


mydata2 = mydata1.T
mydata2.index = regionlist

with ExcelWriter(r"D:\data\altas\LHAclassfy\somafig\cluster\heatmapdata.xlsx") as writer:
    mydata2.to_excel(writer)

dataarray = np.array(mydata2)

scnew = []
for i in range(0,len(sclist)):
    if sclist[i]==1:
        scnew.append(2)
    else:
        scnew.append(0)

scarray = np.array(scnew)
dataarraynew = np.r_[dataarray,[scarray]]

mydata3 = pd.DataFrame(dataarraynew)
regionlist.append('SC_C1')
mydata3.index = regionlist
mydata4 = mydata3.assign(group=regionlist)



colorlist = ['darkgreen']*7+['steelblue']*7+['darkred']*2+['darkmagenta']*13+['peru']*7+['y']*6+['grey']*5



colors = [(0,'#FFFFFF'),(0.4,'#FF4040'),(1,'#FF0000')]
newcmap = LinearSegmentedColormap.from_list('custom_cmap',colors)
xlist = [86,148,263,479,683,742]
ylist = [7,14,16,29,36,42]
fig2 = plt.figure(figsize=(12,8),dpi=300)
ax = plt.gca()
plt.rcParams.update({'font.family':'Arial'})
sns.heatmap(mydata4.iloc[:,0:916],yticklabels=True,
            cmap=newcmap,xticklabels=False,cbar_kws={'shrink':0.5})
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(True)
for i in range(0,45):
    ax.get_yticklabels()[i].set_color(colorlist[i])

for j in ylist:
    plt.axhline(y=j,color='grey',linestyle='--')
for h in xlist:
    plt.axvline(x=h,color='red',linestyle='--')
plt.show()

























savepath = "D:\\data\\altas\\20240714_LHA\\morphology_scene"
clusternumber=[29,57,4,58,44,34,37,65,50,101,37,23,
               40,34,70,6,7,29,17,16,43,47,35,33]

num=0
for i in range(0,24):
    scenelist = []
    for h in range(num,num+clusternumber[i]):
        name1,name2 = lhasamplelist[h].split('_')
        sampleid = name1
        name = name2+'.swc'
        scenedict = {'mirror':0,'sampleid':sampleid,'name':name,'soma':'LHA'}
        scenelist.append(scenedict)
    num=num+clusternumber[i]
    Scene.createScene(scenelist,savepath+"\\cluster{name}.nv".format(name=i+1))
