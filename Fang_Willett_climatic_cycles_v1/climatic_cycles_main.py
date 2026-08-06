# This py.file is the main codes to calculate Ecos-qm
# Author @X. Fang Aug. 2026
# #######################################################
# needed packages
import climatic_cycles_functions as ccf
import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio import features
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# #######################################################
# Section One: extend CHELSA TraCE21k precipitation from 21 kyr to ~100 kyr based on fits to insolation or benthic proxies
# 1 read insolation and benthic data
# insolation from Laskar 2004
insolation = pd.read_excel('proxy_insolation_Laskar2004_10Myr.xls')
insolation_data = np.array(insolation.iloc[-131:, 1])
insolation_t=np.array(insolation.iloc[-131:, 0])*(-1)
# benthic data from Lisiecki 2005
benthic=pd.read_excel('proxy_benthic_Lisiecki2005.xls')
benthic_data=np.array(benthic.iloc[:101, 1])
benthic_t=np.array(benthic.iloc[:101, 0])
benthic_data=benthic_data[::-1]
benthic_t=benthic_t[::-1]
# Linearly interpolate the insolation and benthic data to match the resolution of CHELSA traCE21k precipitation
insolation_data_h=ccf.interpolate_array(insolation_data, 9)
insolation_t_h=ccf.interpolate_array(insolation_t, 9)
benthic_data_h=ccf.interpolate_array(benthic_data, 9)
benthic_t_h=ccf.interpolate_array(benthic_t, 9)

# 2 Paleo precipitation from CHELSA traCE21k
# this precipitation is from 21ka to present
filename = 'CHELSA_TraCE21k_bio12_-200_V1.0891022337.tif'
P_data, P_bands, P_proj, firstX, firstY, nx, ny, intervalX, intervalY = ccf.getTiffdata(filename)
# put all precipitation data from 21ka to present into Ptotal
Ptotal=np.zeros((P_data.shape[0],P_data.shape[1],221))
for num in range(-200, 21, 1):
    print(num)
    filename = 'CHELSA_TraCE21k_bio12_' + str(num) + '_V1.0891022337' + '.tif'
    P_data, P_bands, P_proj, firstX, firstY, nx, ny, intervalX, intervalY = ccf.getTiffdata(filename)
    Ptotal[:,:,num+200]=P_data
# age array, from past to present, unit kyr
age = np.array([round(-0.1 * i + 2, 1) for i in range(-200, 21, 1)])

# 3 read Minindex, which is from climatic_cycles_preparation.py
filename='CHELSA_TraCE21k_891022337_Minindex_precipitation.tif'
with rasterio.open(filename) as src:
    my_transform = src.transform  # Affine transform
chelsaMinindex, im_bands, im_proj, firstX, firstY, nx, ny, intervalX, intervalY = ccf.getTiffdata(filename)
# Determine whether the precipitation cycles at each pixel resemble insolation or benthic curves
# If Minindex < 16, the precipitation pattern resembles insolation and the phase needs to be adjusted
ins_or_ben_check=np.zeros(chelsaMinindex.shape)
ins_or_ben_check[chelsaMinindex<16]=1
phase_change=(20-chelsaMinindex)*ins_or_ben_check

# 4 fitting results between CHELSA traCE21K precipiation and insolation/benthic data
bestA=np.zeros(chelsaMinindex.shape) # best fitting slope for every pixel
bestB=np.zeros(chelsaMinindex.shape) # best fitting intercept for every pixel
Plength=len(age)
for temy in range(chelsaMinindex.shape[0]):
    print(temy,'/',chelsaMinindex.shape[0])
    for temx in range(chelsaMinindex.shape[1]):
        temP=Ptotal[temy,temx,:] # from old time to present time
        if ins_or_ben_check[temy][temx]==0: # benthic
            temproxy = benthic_data_h[-Plength:]
        else: # insolation
            temphasechange=int(phase_change[temy][temx]*10)
            temproxy=insolation_data_h[-Plength-temphasechange:-temphasechange]
        A, B = np.polyfit(temproxy, temP, 1)
        bestA[temy][temx] = A
        bestB[temy][temx] = B

# 5 calculate and write fitted Pa,Pm and Pa/Pm
fitPm=np.zeros(chelsaMinindex.shape)
fitPa=np.zeros(chelsaMinindex.shape)
for temy in range(chelsaMinindex.shape[0]):
    for temx in range(chelsaMinindex.shape[1]):
        if ins_or_ben_check[temy][temx] == 0:  # benthic
            temPline=bestA[temy][temx]*benthic_data_h[:]+bestB[temy][temx]
            temPm=np.mean(temPline)
            temPa=(np.max(temPline)-np.min(temPline))/2
        else:# insolation
            temphasechange = int(phase_change[temy][temx] * 10)
            temPin=insolation_data_h[-Plength-temphasechange:-temphasechange]
            temPline=bestA[temy][temx] * temPin + bestB[temy][temx]
            temPm = np.mean(temPline)
            temPa = (np.max(temPline) - np.min(temPline)) / 2
        fitPm[temy][temx]=temPm
        fitPa[temy][temx]=temPa
# write fitted Pa,Pm and Pa/Pm
# Pm
im_width, im_height = fitPm.shape[1], fitPm.shape[0]
part_tVect = (firstX-intervalX/2, intervalX, 0.0, firstY-intervalY/2, 0.0, intervalY)
path='CHELSA_TraCE21k_891022337_Pm_fit_precipitation.tif'
ccf.writeTiff(fitPm,im_width, im_height,im_bands,part_tVect,im_proj,path)
ccf.generate_tfw(path,'tfw')
# Pa
path='CHELSA_TraCE21k_891022337_Pa_fit_precipitation.tif'
ccf.writeTiff(fitPa,im_width, im_height,im_bands,part_tVect,im_proj,path)
ccf.generate_tfw(path,'tfw')
# Pa/Pm
fitPa_Pm=fitPa/fitPm
path='CHELSA_TraCE21k_891022337_Pa_Pm_fit_precipitation.tif'
ccf.writeTiff(fitPa_Pm,im_width, im_height,im_bands,part_tVect,im_proj,path)
ccf.generate_tfw(path,'tfw')

# #######################################################
# Section Two: basin average precipitation history
# read basin shape files
shp_path ="Granite_basin.shp"
gdf_tem = gpd.read_file(shp_path)
gdf = gdf_tem.sort_values(by="Id").reset_index(drop=True) # make sure your shapefiles have Id and it is sorted by Id

# fitting precipitation history for every basin
withinPolyList=[]
for idx, row in gdf.iterrows():
    geom = row.geometry
    # Create a mask where True indicates pixels inside the basin polygon
    mask = features.geometry_mask([geom], transform=my_transform, invert=True, out_shape=chelsaMinindex.shape)
    # Get the row and column indices of all pixels within the polygon
    indices = np.argwhere(mask)
    withinPolyList.append(indices)

Pfitseries=np.zeros((len(gdf),len(benthic_data_h)))
for i in range(len(withinPolyList)):
    tempointsN=len(withinPolyList[i])
    tempoints=withinPolyList[i]
    temPfit=np.zeros(len(benthic_data_h))
    for j in range(tempointsN):
        temy=tempoints[j][0]
        temx=tempoints[j][1]
        if ins_or_ben_check[temy][temx] == 0:  # benthic
            temPm = bestA[temy][temx] * benthic_data_h[:] + bestB[temy][temx]
        else:  # insolation
            temphasechange = int(phase_change[temy][temx] * 10)
            temPin = insolation_data_h[-len(benthic_data_h) - temphasechange:-temphasechange]
            temPm = bestA[temy][temx] * temPin + bestB[temy][temx]
        temPfit+=temPm
    temPfit=temPfit/tempointsN
    Pfitseries[i][:]=temPfit

# #######################################################
# Section Three: basin information and the calculation of Ecos-qm
# read basin information
# including 10Be concentration, ksn,ksnq, ksn-qm using topotoolbox, average production rates, Pn, Pms, and Pmf
sample_table = pd.read_excel('Granite_basin.xls')
sample_Be10 = np.array(sample_table.iloc[:, 5])
sample_ero=np.array(sample_table.iloc[:, 7]) # Ecos, unit is mm/yr
sample_ksn = np.array(sample_table.iloc[:, 9])
sample_ksnq = np.array(sample_table.iloc[:, 10])
sample_ksnqm = np.array(sample_table.iloc[:, 11])
sample_Pn = np.array(sample_table.iloc[:, 12])
sample_Pms = np.array(sample_table.iloc[:, 13])
sample_Pmf = np.array(sample_table.iloc[:, 14])

rock_density=2.65 # g/cm3
Half_10Be=1.387e6 # y
L_10Be=np.log(2)/Half_10Be # namuda
attL_spal_10Be=160 # g/cm2 Spallation attenuation length (Braucher et al., 2011) [g/cm2]
attL_ms_10Be = 1500  # Slow muons attenuation length (Braucher et al., 2011) [g/cm2]
attL_mf_10Be = 4320 # Fast muons attenuation length (Braucher et al., 2011) [g/cm2]

# the key part, calculating Ecos-qm
# Ecosmodi_tList is the Ecos-qm for every basin, the unit is m/yr, and bestc_tList is the opitimal c for every basin
Ecosmodi_tList=np.zeros(sample_Be10.size)
bestc_tList=np.zeros(sample_Be10.size)
for tem_index in range(len(sample_Be10)):
    # get fitted precipitation for every basin
    if tem_index<20:
        temP=Pfitseries[tem_index][-921:]
        temage=np.array([100*i for i in range(921)])
    else:
        temP = Pfitseries[tem_index][-1001:]
        temage = np.array([100 * i for i in range(1001)])

    temtemEin = np.power(temP/1e3, 0.5) / 1e3
    temC10Be=sample_Be10[tem_index]
    P_spal_SLHL, P_ms_SLHL,P_mf_SLHL=sample_Pn[tem_index],sample_Pms[tem_index],sample_Pmf[tem_index]
    # temtemEin * best_c is the Ecos-qm
    # this error_func and minimize_scalar are to find the best c
    def error_func(c):
        temEin = temtemEin * c
        C10Be_last = ccf.nuclide_calculation(temage, temEin, len(temEin) - 2,
                                         P_spal_SLHL, P_mf_SLHL, P_ms_SLHL,
                                         attL_spal_10Be, attL_mf_10Be, attL_ms_10Be,
                                         rock_density, L_10Be)
        return abs(C10Be_last[-1] - temC10Be)
    res = minimize_scalar(error_func, bounds=(1e-4, 10), method='bounded')
    best_c = res.x
    best_error = res.fun
    temEcosmodi = np.mean(temtemEin * best_c)
    Ecosmodi_tList[tem_index] = temEcosmodi
    bestc_tList[tem_index] = best_c

# #######################################################
# Section Four: figure
need_index1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
need_index2=[20,21,22,23,24]
ksndiffer=np.abs((sample_ksnq-sample_ksnqm)/sample_ksnqm)
erodiffer=np.abs((sample_ero-Ecosmodi_tList*1e3)/(Ecosmodi_tList*1e3))

plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
fig = plt.figure(figsize=(8,3))
ax1=fig.add_subplot(1,2,1)
ax1.scatter(sample_ksnq[need_index1],sample_ero[need_index1],c='purple',marker='o',edgecolors='k')
ax1.scatter(sample_ksnq[need_index2],sample_ero[need_index2],c='purple',marker='d',edgecolors='k')
ax1.scatter(sample_ksnqm[need_index1],Ecosmodi_tList[need_index1]*1e3,c='y',marker='o',edgecolors='k')
ax1.scatter(sample_ksnqm[need_index2],Ecosmodi_tList[need_index2]*1e3,c='y',marker='d',edgecolors='k')
ax1.set_yticks([0.1,0.3,0.5,0.7])
ax1.set_xticks([200,400,600,800])
ax1.set_ylabel('Erosion rate (mm/yr)')
ax1.set_xlabel('Channel steepness')
ax1.xaxis.set_major_formatter(ScalarFormatter())
ax1.yaxis.set_major_formatter(ScalarFormatter())

ax2=fig.add_subplot(1,2,2)
ax2.scatter(ksndiffer[need_index1]*100,erodiffer[need_index1]*100,c='grey',marker='o',edgecolors='k')
ax2.scatter(ksndiffer[need_index2]*100,erodiffer[need_index2]*100,c='grey',marker='d',edgecolors='k')
ax2.set_xlim(-1,10)
ax2.set_ylim(-1,22)
ax2.set_xticks([0,5,10])
ax2.set_yticks([0,5,10,15,20])
ax2.set_ylabel('Erosion rate deviation %')
ax2.set_xlabel('Channel steepness deviation (%)')
plt.savefig('channel steepness vs erosion rates.jpg',dpi=400,bbox_inches='tight')
plt.show()


