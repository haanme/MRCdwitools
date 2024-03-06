#!/usr/bin/env python

#
# Radiomics for Medical Imaging - Nifti to ASCII conversion.
#
# Copyright (C) 2019-2022 Harri Merisaari haanme@utu.fi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Created by Harri Merisaari on 2.6.2022.
#

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
from bvecIO_E import *
from glob import glob
import numpy as np


def find_cropping(img,bg_threshold):
    import numpy as np
    
    x_lo = 0
    for x_i in range(img.shape[0]):
        max_val = np.amax(img[x_i,:,:])
        if max_val >= bg_threshold:
            x_lo = x_i
            break
    x_hi = img.shape[0]-1
    for x_i in range(img.shape[0]-1,-1,-1):
        max_val = np.amax(img[x_i,:,:])
        if max_val >= bg_threshold:
            x_hi = x_i
            break
    y_lo = 0
    for y_i in range(img.shape[1]):
        max_val = np.amax(img[:,y_i,:])
        if max_val >= bg_threshold:
            y_lo = y_i
            break
    y_hi = img.shape[1]-1
    for y_i in range(img.shape[1]-1,-1,-1):
        max_val = np.amax(img[:,y_i,:])
        if max_val >= bg_threshold:
            y_hi = y_i
            break
    z_lo = 0
    for z_i in range(img.shape[2]):
        max_val = np.amax(img[:,:,z_i])
        if max_val >= bg_threshold:
            z_lo = z_i
            break
    z_hi = img.shape[2]-1
    for z_i in range(img.shape[2]-1,-1,-1):
        max_val = np.amax(img[:,:,z_i])
        if max_val >= bg_threshold:
            z_hi = z_i
            break
    return [x_lo, x_hi, y_lo, y_hi, z_lo, z_hi]


def load_nifti(name, filename):
    print('Loading ' + name + ':' + filename)
    img = nib.load(filename)
    affine = img.affine
    data = np.asanyarray(img.dataobj)
    data_sform = img.header.get_sform()
    voxelsize = [abs(data_sform[0, 0]), abs(data_sform[1, 1]), abs(data_sform[2, 2])]
    return data, affine, voxelsize


def NIFTI2ASCII_MASK(DWIfile, maskfile, bvalfile, subwindow, bvalselection):

    outfile = os.path.splitext(DWIfile)[0] + '_ASCII.txt'
    
    data, affine, voxelsize = load_nifti('DWI file', DWIfile)
    bvals, bval_string = read_bvals(bvalfile)
    if len(data.shape) < 4:
        return 'ERROR: DWI does not have 4th dimension'
    tdim = data.shape[3]
    if not tdim == len(bvals):
        return 'ERROR: DWI b-values ' + str(tdim) + ' while expected b-values based on b-value file were ' + str(len(bvals))
    if len(bvalselection) > 0:
        tdim = len(bvalselection)
        bvals_selected = []
        for bval_i in range(len(bvalselection)):
            bvals_selected.append(bvals[bvalselection[bval_i]])
        bvals = bvals_selected
    print('Selected b-values [' + str(len(bvals)) + ']:', str(bvals))

    mdata, maffine, mvoxelsize = load_nifti('mask file', maskfile)

    ROI_No = 0
    if len(subwindow) == 0:
        subwindow = [0, 0, 0, 0]
        ROIslice = [0, data.shape[2]-1]
    else:
        ROIslice = [subwindow[4], subwindow[5]]
        subwindow = [subwindow[0], subwindow[1], subwindow[2], subwindow[3]]
    name = os.path.splitext(os.path.basename(DWIfile))[0]

    SIs = np.zeros([len(mdata[mdata>0]), tdim])
    for t in range(tdim):
        if len(bvalselection) > 0:
            vol = data[:, :, :, bvalselection[t]]
        else:
            vol = data[:, :, :, t]
        SIs[:, t] = vol[mdata>0]

    data = { 'subwindow': subwindow, 'ROI_No': ROI_No, 'bset': bvals, 'ROIslice': ROIslice, 'name': name, 'SIs': SIs, 'number': 0 }
    ASCIIio = bfitASCII_IO.bfitASCII_IO()
    print("Writing " + outfile)
    ASCIIio.Write3D(outfile, data)

    return outfile


def check_exists(name,filename,files_missing):
    if not os.path.exists(filename):
        print(name + ' does not exist [' + filename + ']')
        files_missing = files_missing + 1
    return files_missing


if __name__ == "__main__":
    # Parse input arguments into args structure
    parser = ArgumentParser()
    parser.add_argument("--DWIfile", dest="DWIfile", help="DWI Nifti file", required=True)
    parser.add_argument("--maskfile", dest="maskfile", help="Nifti mask file", required=False, default='')
    parser.add_argument("--centercrop", dest="centercrop", help="Crop center NxN region, default [-1 no crop]", required=False, default='-1')
    parser.add_argument("--th", dest="th", help="Threshold value to create Nifti mask file, only if centercrop not given", required=False, default='100')
    parser.add_argument("--bvalfile", dest="bvalfile", help="ASCII file with b-values in one line space separated", required=True)
    parser.add_argument("--bvalselection", dest="bvalselection", help="0-based indexes of selected b-values in format [i1,i2,..,iN], default [empty]=all", required=False, default='')
    args = parser.parse_args()
    DWIfile = args.DWIfile
    maskfile = args.maskfile
    bvalfile = args.bvalfile
    centercrop = int(float(args.centercrop))
    if len(args.bvalselection) == 0:
        bvalselection = []
    else:
        bvalselection = [int(x) for x in args.bvalselection.replace(']','').replace('[','').split(',')]
    th = int(float(args.th))

    files_missing = 0
    files_missing = check_exists('DWI file', DWIfile, files_missing)
    files_missing = check_exists('bval file', bvalfile, files_missing)
    if files_missing > 0:
        sys.exit(1)
    
    # Create mask file if it does not exist
    subwindow = []
    files_missing = 0
    files_missing = check_exists('Mask file', maskfile, files_missing)
    if files_missing > 0:
        data, affine, voxelsize = load_nifti('DWI file', DWIfile)
        if centercrop == -1:
            print('Mask not found, using threshold ' + str(th) + ' for b0 (1st index)')
            x_lo, x_hi, y_lo, y_hi, z_lo, z_hi = find_cropping(data[:, :, :, 0], th)
        else:
            print('Mask not found, using center crop ' + str(centercrop))
            center_x = int(np.floor(data.shape[0] / 2))
            center_y = int(np.floor(data.shape[1] / 2))
            x_lo = center_x - int(np.floor(centercrop / 2))
            y_lo = center_x - int(np.floor(centercrop / 2))
            z_lo = 0
            x_hi = x_lo + centercrop
            y_hi = y_lo + centercrop
            z_hi = z_lo + data.shape[2]
        print('Using cropping ' + str((x_lo, x_hi, y_lo, y_hi, z_lo, z_hi)))
        mask = np.zeros([data.shape[0], data.shape[1], data.shape[2]])
        subwindow = [x_lo, x_hi, y_lo, y_hi, z_lo, z_hi]
        print((subwindow))
        mask[x_lo:x_hi, y_lo:y_hi, z_lo:z_hi] = 1
        print(np.amin(mask))
        print(np.amax(mask))
        img = nib.Nifti1Image(mask, affine)
        maskfile = os.path.splitext(DWIfile)[0] + '_mask.nii.gz'
        print('Saving mask ' + maskfile)
        nib.save(img, maskfile)
    
    NIFTI2ASCII_MASK(DWIfile, maskfile, bvalfile, subwindow, bvalselection)
