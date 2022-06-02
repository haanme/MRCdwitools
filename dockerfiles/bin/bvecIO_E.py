#!/usr/bin/env python
#-*- coding: UTF-8 -*-

#
# Radiomics for Medical Imaging - b-vector file I/O.
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

import os
import numpy as np

# Reads fsl rotation matrix
def read_mat(mat_file):
    f = open(mat_file, 'r')
    mat = np.zeros([4, 4])
    line_i = 0
    for line in f:
        l = repr(line)
        splitted = line.split(' ')
        values = []
        for s_i in range(len(splitted)):
            if len(splitted[s_i].strip()) > 0:
                values.append(float(splitted[s_i].strip()))
        mat[line_i,:] = values
        line_i = line_i + 1
    f.close()
    return mat

def read_bvals(bval_file_orig):

    root, ext2 = os.path.splitext(bval_file_orig)

    f = open(bval_file_orig, 'r')
    bval_string = []
    line_i = 0
    for line in f:
        l = repr(line)
        splitted = line.split(' ')
        bvals = []
        for s_i in range(len(splitted)):
            if len(splitted[s_i].strip()) == 0:
                continue
            bvals.append(float(splitted[s_i]))
            bval_string = splitted[s_i]
        line_i = line_i + 1
        break
    f.close()

    return bvals, bval_string

def read_bvecs(bvec_file_orig):

    root, ext2 = os.path.splitext(bvec_file_orig)

    f = open(bvec_file_orig, 'r')
    bvec_strings = [[],[],[]]
    line_i = 0
    for line in f:
        l = repr(line)
        splitted = line.split(' ')
        values = np.zeros([len(splitted)])
        for s_i in range(len(splitted)):
            if len(splitted[s_i].strip()) == 0:
                continue
            values[s_i] = float(splitted[s_i])
            bvec_strings[line_i].append(splitted[s_i])
        line_i = line_i + 1
    f.close()

    bvecs = np.zeros([3,len(bvec_strings[0])])
    for b_i in range(bvecs.shape[1]):
        bvecs[:,b_i] = [float(bvec_strings[0][b_i]), float(bvec_strings[1][b_i]), float(bvec_strings[2][b_i])]
    bvecs = bvecs.T

    return bvecs, bvec_strings

def read_bvecs_Transpose(bvec_file_orig):

    root, ext2 = os.path.splitext(bvec_file_orig)

    f = open(bvec_file_orig, 'r')
    bvec_strings = []
    line_i = 0
    for line in f:
        l = repr(line)
        splitted = line.split(' ')
        for s_i in range(len(splitted)):
            if len(splitted[s_i].strip()) == 0:
                continue
            bvec_strings.append(splitted[s_i])
        line_i = line_i + 1
    f.close()

    bvecs = np.zeros([3, line_i])
    for b_i in range(line_i):
        bvecs[:, b_i] = [float(bvec_strings[b_i*3]), float(bvec_strings[b_i*3+1]), float(bvec_strings[b_i*3+2])]
    bvecs = bvecs.T

    return bvecs, bvec_strings

def write_bvecs(bvec_file, bvecs):

    f_bvec_0 = open(bvec_file, 'w')
    bvals = bvecs.shape[1]

    for bvec_i in range(bvals):
        f_bvec_0.write(' %s' % str(bvecs[0][bvec_i]))
    f_bvec_0.write('\n')
    for bvec_i in range(bvals):
        f_bvec_0.write(' %s' % str(bvecs[1][bvec_i]))
    f_bvec_0.write('\n')
    for bvec_i in range(bvals):
        f_bvec_0.write(' %s' % str(bvecs[2][bvec_i]))
    f_bvec_0.write('\n')
    f_bvec_0.close()

    return bvec_file

def write_bvecs_Transpose(bvec_file, bvecs):

    import nibabel as nib

    f_bvec = open(bvec_file, 'w')
    bvals = bvecs.shape[0]

    for bvec_i in range(bvals):
        f_bvec.write('%s %s %s\n' % (str(bvecs[bvec_i][0]), str(bvecs[bvec_i][1]), str(bvecs[bvec_i][2])))
    f_bvec.close()

    return bvec_file

def split_bvecs_0_64(bvec_file_orig):

    import nibabel as nib

    root, ext2 = os.path.splitext(bvec_file_orig)

    bval_file = experiment_dir + os.sep + subjectname + os.sep + subjectname + '_fixed.bval'
    bvec_file = experiment_dir + os.sep + subjectname + os.sep + subjectname + '_fixed.bvec'

    f = open(bvec_file_orig, 'r')
    read_state = 0
    gradient_i = 1
    included_gradients = []
    included_gradients.append([0, 0])
    bvecs = np.zeros([3, 64])
    bvec_strings = [[],[],[]]
    line_i = 0
    for line in f:
        l = repr(line)
        splitted = line.split(' ')
        values = np.zeros([len(splitted)-1])
        for s_i in range(len(splitted)-1):
            values[s_i] = float(splitted[s_i])
            bvec_strings[line_i].append(splitted[s_i])
        bvecs[line_i, :] = values
        line_i = line_i + 1
    f.close()
    bvecs = bvecs.T

    bval_file_0 = root + '_DWI_0' + '.bval'
    bvec_file_0 = root + '_DWI_0' + ext2

    f_bval_0 = open(bval_file_0, 'w')
    f_bvec_0 = open(bvec_file_0, 'w')

    f_bval_0.write('0')
    for bvec_i in range(64):
        if bvec_i > 0:
            f_bval_0.write(' %f' % 1000.0)
        f_bvec_0.write(' %s' % bvec_strings[0][bvec_i])
    f_bval_0.write('\n')
    f_bvec_0.write('\n')
    f_bval_0.close()

    for bvec_i in range(64):
        f_bvec_0.write(' %s' % bvec_strings[1][bvec_i])
    f_bvec_0.write('\n')

    for bvec_i in range(64):
        f_bvec_0.write(' %s' % bvec_strings[2][bvec_i])
    f_bvec_0.write('\n')
    f_bvec_0.close()

    return bval_file_0, bvec_file_0
