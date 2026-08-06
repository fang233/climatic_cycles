# This py.file is functions for preparation and main.py files
# Author @X. Fang Aug. 2026
import numpy as np
from osgeo import gdal
import osgeo.osr as osr
import os

# Function getTiffdata is to read tiff files
def getTiffdata(filename):
    gdal.AllRegister()
    dataset = gdal.Open(filename, gdal.GA_ReadOnly)
    nx = dataset.RasterXSize
    ny = dataset.RasterYSize
    im_bands = dataset.RasterCount  # band number
    im_proj = dataset.GetProjection()  # projection information
    im_data = dataset.ReadAsArray(0, 0, nx, ny)  # data
    tVect = dataset.GetGeoTransform()
    leftUpX,leftUpY = tVect[0],tVect[3]
    firstX,firstY=leftUpX+tVect[1]/2,leftUpY+tVect[5]/2
    return im_data,im_bands,im_proj,firstX,firstY,nx,ny,tVect[1],tVect[5]

# Function writeTiff is to write tiff files
def writeTiff(im_data,im_width, im_height,im_bands,im_geotrans,im_proj,path):
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32

    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    elif len(im_data.shape) == 2:
        im_data = np.array([im_data])
    else:
        im_bands, (im_height, im_width) = 1,im_data.shape
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_height, im_bands, datatype)
    if(dataset!= None):
        dataset.SetGeoTransform(im_geotrans)
        dataset.SetProjection(im_proj)
    for i in range(im_bands):
        dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset

# Function generate_tfw is to generate coresponding tfw files
def generate_tfw(infile, gen_prj):
    src = gdal.Open(infile)
    xform = src.GetGeoTransform()

    if gen_prj == 'prj':
        src_srs = osr.SpatialReference()
        src_srs.ImportFromWkt(src.GetProjection())
        src_srs.MorphToESRI()
        src_wkt = src_srs.ExportToWkt()

        prj = open(os.path.splitext(infile)[0] + '.prj', 'wt')
        prj.write(src_wkt)
        prj.close()

    src = None
    edit1=xform[0]+xform[1]/2
    edit2=xform[3]+xform[5]/2

    tfw = open(os.path.splitext(infile)[0] + '.tfw', 'wt')
    tfw.write("%0.8f\n" % xform[1])
    tfw.write("%0.8f\n" % xform[2])
    tfw.write("%0.8f\n" % xform[4])
    tfw.write("%0.8f\n" % xform[5])
    tfw.write("%0.8f\n" % edit1)
    tfw.write("%0.8f\n" % edit2)
    tfw.close()

# Function interpolate_array is to insert n_points linearly interpolated points between each pair of adjacent elements in the array arr.
def interpolate_array(arr, n_points):
    result = []
    for i in range(len(arr) - 1):
        # For each segment, insert n+2 points (including both endpoints), and then remove the last point to avoid duplication.
        interpolated = np.linspace(arr[i], arr[i + 1], n_points + 2)[:-1]
        result.extend(interpolated)
    result.append(arr[-1])
    return np.array(result)

# Function nuclide_calculation computes the present-day surface 10Be concentration for a given erosion rate history
def nuclide_calculation(tc,erosion_rates,start_mear_index,P1,P2,P3,attL1,attL2,attL3,density,L_10Be):
    timeStep = tc[1] - tc[0] # unit: yr
    # Spallation, slow muons, fast muons respectively
    L1, L2, L3 = attL1 / density, attL2 / density, attL3 / density

    # surface concentration at time t for Spallation, slow muons, fast muons respectively
    Css_final1, Css_final2, Css_final3 = np.zeros(tc.shape), np.zeros(tc.shape), np.zeros(tc.shape)

    # concentration start from 0
    # erosion depth for every time step
    ero_h = erosion_rates * timeStep*100  # cm
    # erosion depth for the whole history
    ero_h_acc = np.cumsum(ero_h)  # cm

    for i_index in range(start_mear_index, len(tc) - 1, 1):
        # print(i_index,'/',len(tc)-1)
        ero_h_acc_tem = ero_h_acc[:i_index + 1]
        ero_h_acc_tem = np.insert(ero_h_acc_tem, 0, 0)
        # erosion depth
        # depth
        ero_h_acc_tem = ero_h_acc_tem[-1] - ero_h_acc_tem
        ero_h_acc_tem = ero_h_acc_tem[:-1]
        #
        tem_coff = L_10Be * timeStep - 1
        powers_of_temcoff = tem_coff ** np.arange(i_index + 1)
        signs = (-1) ** np.arange(i_index + 1)

        temP1 = P1 * np.exp(-ero_h_acc_tem / L1) * timeStep
        temP2 = P2 * np.exp(-ero_h_acc_tem / L2) * timeStep
        temP3 = P3 * np.exp(-ero_h_acc_tem / L3) * timeStep

        temP1_rever = np.flip(temP1)
        temP2_rever = np.flip(temP2)
        temP3_rever = np.flip(temP3)

        Css_final1[i_index + 1] = np.dot(signs * powers_of_temcoff, temP1_rever)
        Css_final2[i_index + 1] = np.dot(signs * powers_of_temcoff, temP2_rever)
        Css_final3[i_index + 1] = np.dot(signs * powers_of_temcoff, temP3_rever)
    Css_total = Css_final1 + Css_final2 + Css_final3
    return Css_total
