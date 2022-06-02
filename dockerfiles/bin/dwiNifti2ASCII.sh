#!/bin/sh

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

dwi=$1
mask=$2
th=$3
if [[ -z "$dwi" ]]
then
    echo "DWI Nifti filename was not given"
fi
if [[ -z "$mask" ]]
then
    echo "Mask Nifti filename was not given"
fi
if [[ -z "$th" ]]
then
    python /usr/local/bin/dwiNifti2ASCII.py --DWIfile $dwi --maskfile $mask
else
    python /usr/local/bin/dwiNifti2ASCII.py --DWIfile $dwi --maskfile $mask --th $th
fi
