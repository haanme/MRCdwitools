#!/usr/bin/env python

#
# Radiomics for Medical Imaging - Nifti reconstruction from ASCII.
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

from argparse import ArgumentParser
import sys
import os
import nibabel as nib
import bfitASCII_IO_B as bfitASCII_IO
import numpy as np
from glob import glob

###############
# Main script #
###############
def load_nifti(name, filename):
    print('Loading ' + name + ':' + filename)
    img = nib.load(filename)
    affine = img.affine
    data = np.asanyarray(img.dataobj)
    data_sform = img.header.get_sform()
    voxelsize = [abs(data_sform[0, 0]), abs(data_sform[1, 1]), abs(data_sform[2, 2])]
    return data, affine, voxelsize


def fun_flipud_DATA(img):
    print('fun_flipud')
    imgout = np.zeros_like(img)
    for z in range(img.shape[2]):
        imgout[:,:,z] = np.flipud(img[:,:,z])
    return imgout


def ASCII2NIFTI_MASK(DWIfile, maskfile, ASCIIfile):

    ASCIIio = bfitASCII_IO.bfitASCII_IO()
    print('Reading:' + ASCIIfile)
    suffix = os.path.basename(ASCIIfile).split('_')[-2]
    print(suffix)
    data = ASCIIio.Read(ASCIIfile, False)
    pnames = data['parameters']
    print(pnames)
    data_data = data['data']
    print(data_data.shape)

    data, affine, voxelsize = load_nifti('DWI file', DWIfile)
    mdata, maffine, mvoxelsize = load_nifti('Mask file', maskfile)
    op = []
    filenames = []
    for i in range(len(pnames)):
        pname = pnames[i].strip('[').strip(']').replace('\'','').strip(',')
        filename = os.path.splitext(ASCIIfile)[0] + '_' + pname + '.nii.gz'
        print("Saving:" + filename)
        imgdata = np.zeros([data.shape[0], data.shape[1], data.shape[2]])
        imgdata[mdata > 0] = data_data[:, i]
        if 'flipud' in op:
            imgdata = fun_flipud_DATA(imgdata)
        img = nib.Nifti1Image(imgdata, affine)
        nib.save(img, filename)
        filenames.append(filename)

    return filenames

def check_exists(name,filename,files_missing):
    if not os.path.exists(filename):
        print(name + ' does not exist [' + filename + ']')
        files_missing = files_missing + 1
    return files_missing


if __name__ == "__main__":
    # Parse input arguments into args structure
    parser = ArgumentParser()
    parser.add_argument("--DWIfile", dest="DWIfile", help="DWI Nifti file", required=True)
    parser.add_argument("--maskfile", dest="maskfile", help="Nifti mask file", required=True)
    parser.add_argument("--ASCIIfile", dest="ASCIIfile", help="ASCII file with fitting results", required=True)
    args = parser.parse_args()
    DWIfile = args.DWIfile
    maskfile = args.maskfile
    ASCIIfile = args.ASCIIfile

    files_missing = 0
    files_missing = check_exists('DWI file', DWIfile, files_missing)
    files_missing = check_exists('Mask file', maskfile, files_missing)
    files_missing = check_exists('ASCII file', ASCIIfile, files_missing)
    if files_missing > 0:
        sys.exit(1)
    
    ASCII2NIFTI_MASK(DWIfile, maskfile, ASCIIfile)
