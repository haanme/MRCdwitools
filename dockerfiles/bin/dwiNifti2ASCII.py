#!/usr/bin/env python

####################################################################
# Python 2.7 script for executing FA, MD calculations for one case #
####################################################################

# Directory where result data are located
experiment_dir = ''
subjectname = ''
avg_suffix = ''
avg_suffix_noavg = ''

from argparse import ArgumentParser
import sys
import os
import nibabel as nib
import bfitASCII_IO_B as bfitASCII_IO
from glob import glob
import numpy as np


def load_nifti(name, filename):
    print('Loading ' + name + ':' + filename)
    img = nib.load(filename)
    affine = img.get_affine()
    data = img.get_data()
    data_sform = img.get_header().get_sform()
    voxelsize = [abs(data_sform[0, 0]), abs(data_sform[1, 1]), abs(data_sform[2, 2])]
    return data, affine, voxelsize


def NIFTI2ASCII_MASK(niiname, bset, bindexes, cutoff, casename, outdir):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outdir = outdir + os.sep + os.sep + casename + '_ASCII.txt'
    print('Output directory will be ' + outdir)
    if os.path.isfile(outdir):
        os.remove(outdir)

    print('Reading:' + niiname)
    b_data = nib.load(niiname).get_data()
    zdim = b_data.shape[2]
    if len(b_data.shape) < 4:
        return 'ERROR: DWI does not have 4th dimension'
    tdim = b_data.shape[3]
    if not tdim == len(bset):
        return 'ERROR: DWI b-values ' + str(tdim) + ' while requested b-values were ' + str(len(bset))

    print(b_data.shape)
    mask = np.zeros_like(np.squeeze(b_data[:,:,:,0]))
    mask[mask.shape[0]*cutoff:mask.shape[0]*(1-cutoff),mask.shape[1]*cutoff:mask.shape[1]*(1-cutoff),:] = 1

    # Save in data in order z,y,x
    SIs = []
    for z_i in range(b_data.shape[2]):
        for y_i in range(b_data.shape[1]):
            for x_i in range(b_data.shape[0]):
                if mask[x_i, y_i, z_i] == 0:
                    continue
                SI = []
                for t_i in bindexes:
                    SI.append(b_data[x_i, y_i, z_i, t_i])
                SIs.append(SI)
        if (z_i % 10 == 0) or (z_i == b_data.shape[2]-1):
            print(str(z_i+1) + os.sep + str(zdim))
    print("total SIs:" + str(len(SIs)))
    subwindow = [0, 0, 0, 0]
    ROI_No = 0
    ROIslice = [0, b_data.shape[2]-1]
    name = subjectname
    data = { 'subwindow': subwindow, 'ROI_No': ROI_No, 'bset': bset, 'ROIslice': ROIslice, 'name': name, 'SIs': SIs, 'number': 0 }
    ASCIIio = bfitASCII_IO.bfitASCII_IO()
    print("Writing " + outdir)
    ASCIIio.Write3D(outdir, data)

    return outdir


def check_exists(name,filename,files_missing):
    if not os.path.exists(filename):
        print(name + ' does not exist [' + filename + ']')
        files_missing = files_missing + 1
    return files_missing


if __name__ == "__main__":
    # Parse input arguments into args structure
    parser = ArgumentParser()
    parser.add_argument("--basedir", dest="basedir", help="base subject directory", required=True)
    parser.add_argument("--outdir", dest="outdir", help="output directory", required=True)
    parser.add_argument("--case", dest="case", help="case", default='', required=False)    
    parser.add_argument("--checkexisting", dest="checkexisting", help="checkexisting", default=0, required=False)    
    args = parser.parse_args()

    bset = [0, 100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2000]
    b_indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    folders = glob(args.basedir + os.sep + '*')
    for folder in folders:
        basename = os.path.basename(folder)
        prefix = basename.split('_')[0]
        if len(args.case) > 0 and not args.case == basename and not args.case == prefix:
            continue
        if not os.path.exists(folder + os.sep + 'DWI.nii'):
            print('Missing ' + folder + os.sep + 'DWI.nii')
            continue
        ret = ''
        if int(float(args.checkexisting)) == 1:
            if not os.path.exists(folder + os.sep + 'DWI12b2000_ADCm.nii'):
                print('Missing ' + folder + os.sep + 'DWI12b2000_ADCm.nii')
                ret = NIFTI2ASCII_MASK(folder + os.sep + 'DWI.nii', bset, b_indexes, 0.25, prefix, args.outdir)
        else:
            ret = NIFTI2ASCII_MASK(folder + os.sep + 'DWI.nii', bset, b_indexes, 0.25, prefix + '_' + os.path.basename(args.basedir).split('_')[-1], args.outdir)
        if len(ret) > 0:
            print(ret)


