# This py.file is preparation for calculating Ecos-qm
# Author @X. Fang Aug. 2026
# #######################################################
# needed packages
import climatic_cycles_functions as ccf
import numpy as np
from scipy.signal import savgol_filter

# #######################################################
# Section One: needed input files
# Paleo precipitation from CHELSA traCE21k
# this precipitation is from 21ka to present
# please download files from https://www.chelsa-climate.org/models/chelsa-trace21k and cut them
# here we only need the range from 89-102E, 23-37N
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

# #######################################################
# Section Two: calculate the required varibles chelsaMinindex,chelsaMaxindex
# find min precipitation and coresponding index, max precipitation and corresponding index during 21ka to present
chelsaMin, chelsaMinindex, chelsaMax, chelsaMaxindex = np.ones(P_data.shape), np.zeros(P_data.shape), np.zeros(P_data.shape), np.zeros(P_data.shape)
latTotal=P_data.shape[0]
for latn in range(P_data.shape[0]):
    print(latn,'/',latTotal)
    for lonn in range(P_data.shape[1]):
        temP=Ptotal[latn,lonn,:]
        # smooth
        # window_length: This parameter defines the size of the sliding window used for fitting. It must be an odd number, typically chosen between 5 and 21.
        # poly_order: This parameter specifies the order of the polynomial used for fitting. Lower orders are generally used for smoothing purposes, while higher orders make the filter more complex. It is usually selected between 0 and 3.
        temPsmooth=savgol_filter(temP, 21, 1, mode= 'nearest')

        temminindex = np.argmin(temPsmooth)
        temmin = np.min(temPsmooth)
        temmaxindex=np.argmax(temPsmooth)
        temmax=np.max(temPsmooth)
        temminindex=age[int(temminindex)]
        temmaxindex=age[int(temmaxindex)]

        chelsaMin[latn][lonn]=temmin
        chelsaMinindex[latn][lonn]=temminindex
        chelsaMax[latn][lonn] = temmax
        chelsaMaxindex[latn][lonn] = temmaxindex

# #######################################################
# Section Three: write the tiff files
# write tiff file minindex, maxindex
P_width, P_height = chelsaMaxindex.shape[1], chelsaMaxindex.shape[0]
part_tVect = (firstX-intervalX/2, intervalX, 0.0, firstY-intervalY/2, 0.0, intervalY)
path='CHELSA_TraCE21k_891022337_Maxindex_precipitation.tif'
ccf.writeTiff(chelsaMaxindex,P_width, P_height,P_bands,part_tVect,P_proj,path)
ccf.generate_tfw(path,'tfw')
path='CHELSA_TraCE21k_891022337_Minindex_precipitation.tif'
ccf.writeTiff(chelsaMinindex,P_width, P_height,P_bands,part_tVect,P_proj,path)
ccf.generate_tfw(path,'tfw')
# #######################################################