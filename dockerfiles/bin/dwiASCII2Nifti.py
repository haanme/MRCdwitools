#!/usr/bin/env python

####################################################################
# Python 2.7 script for executing FA, MD calculations for one case #
####################################################################

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
    affine = img.get_affine()
    data = img.get_data()
    data_sform = img.get_header().get_sform()
    voxelsize = [abs(data_sform[0,0]),abs(data_sform[1,1]),abs(data_sform[2,2])]
    return data, affine, voxelsize


def fun_flipud_DATA(img):
    print('fun_flipud')
    imgout = np.zeros_like(img)
    for z in range(img.shape[2]):
        imgout[:,:,z] = np.flipud(img[:,:,z])
    return imgout


def ASCII2NIFTI_MASK(folder, ASCIIname, bset, cutoff, casename, op, model, modalityprefix):

    ASCIIio = bfitASCII_IO.bfitASCII_IO()
    print('Reading:' + ASCIIname)
    suffix = os.path.basename(ASCIIname).split('_')[-2]
    print(suffix)
    data = ASCIIio.Read(ASCIIname, False)
    pnames = data['parameters']
    print(pnames)

    dwifilename = folder + os.sep + 'DWI.nii'
    print('Reading:' + dwifilename)
    ref_data = nib.load(dwifilename)
    ref_img = ref_data.get_data()
    ref_affine = ref_data.get_affine()
    print(ref_data.shape)

    mask = np.zeros_like(np.squeeze(ref_img[:,:,:,0]))
    mask[int(ref_img.shape[0]*cutoff):int(ref_img.shape[0]*(1-cutoff)),int(ref_img.shape[1]*cutoff):int(ref_img.shape[1]*(1-cutoff)),:] = 1

    filenames = []
    for i in range(0,len(pnames)):
        pname = pnames[i].strip('[').strip(']').replace('\'','').strip(',')
        if pname == 'C' or pname == 'RMSE':
            if suffix == 'Kurt':
                pname += 'k'
            if suffix == 'Mono':
                pname += 'm'
        if not model == 'all' and not suffix == model:
            continue
        print((suffix, model))
        print(pname)
        filename = folder + os.sep + modalityprefix + '_' + pname
        print("Saving:" + filename)
        imgdata = np.zeros([ref_img.shape[0], ref_img.shape[1], ref_img.shape[2]])
        SI_i = 0
        for z_i in range(ref_img.shape[2]):
            for y_i in range(ref_img.shape[1]):
                for x_i in range(ref_img.shape[0]):
                    if mask[x_i, y_i, z_i] == 0:
                        continue
                    imgdata[x_i, y_i, z_i] = data['data'][SI_i, i]
                    SI_i = SI_i + 1
            if (z_i % 10 == 0) or (z_i == ref_img.shape[2]-1):
                print(str(z_i+1) + os.sep + str(ref_img.shape[2]))
        if 'flipud' in op:
            imgdata = fun_flipud_DATA(imgdata)
        img = nib.Nifti1Image(imgdata, ref_affine)
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
    parser.add_argument("--ASCIIfile", dest="ASCIIfile", help="ASCII file with fitting results", required=True)
    args = parser.parse_args()

    bset = [0, 100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2000]
    filenames = glob(args.ASCIIdir + os.sep + '*results.txt')
    print(str(len(filenames)) + ' filenames names from:' + args.ASCIIdir + os.sep + '*results.txt')
    head,tail = os.path.split(args.basedir)
    print((head,tail,os.path.basename(head)))
    if len(args.prefix) == 0:
        modalityprefix = os.path.basename(head).split('_')[2]
    else:
        modalityprefix = args.prefix
    print(modalityprefix)

    folders = glob(args.basedir + os.sep + '*')
    print(str(len(folders)) + ' folder names from:' + args.basedir)
    start = False
    failed_conversions = []
    for filename in filenames:
        basename = os.path.basename(filename)
        prefix = basename.split('_')[0]
        if len(args.case) > 0 and not args.case == basename and not args.case == prefix:
            continue
        for folder in folders:
            basefolder = os.path.basename(folder)
            folderprefix = basefolder.split('_')[0]            
            if not prefix == folderprefix:
                # print(prefix + os.sep + ' does not match ' + os.sep + folderprefix)
                continue
            if not os.path.exists(folder + os.sep + 'DWI.nii'):
                # print(folder + os.sep + 'DWI.nii does not exist')
                continue
            ASCII2NIFTI_MASK(folder, filename, bset, 0.25, prefix, args.op, args.model, modalityprefix)
