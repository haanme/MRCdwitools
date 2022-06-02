#!/bin/sh

#
# Radiomics for Medical Imaging - script for running fitting to all 
# text files with 'ASCII'-files in the current directory.
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
#  Created by Harri Merisaari on 3/10/14.
#

# Binary to execute
binary="dwifit"
# Model name, see dwifit executable
model=$1
if [[ -z "$model" ]]
then
    echo "Model name not given see dwifit help for details"
fi
# Settings filename
settings_file="settings.ini"
# Max CPU usage for starting new fittings automatically
start_cpu=95.0
idle_lim=1
# Execution paths (current directory by default)
spaths="."

echo "Starting "$0" with PID "$(sh -c 'echo $PPID')
echo "Starting "$0" with PID "$(sh -c 'echo $PPID') >> run_dwifit_for_all_ASCII.txt
# Go trough paths
for spath in $spaths
do
echo "Next subfolder:"$spath
# Execute new files until all have corresponding results file with 'results'-suffix
while true
do
    clear
    last | head -5
    echo "---------------------------------------------------------------------"
    top -b -n 2 | grep 'Cpu(s)'
    echo "---------------------------------------------------------------------"
    w
    echo "---------------------------------------------------------------------"
    ps aux | grep dwifit | grep -v 'grep'
    echo "---------------------------------------------------------------------"
    # Get number of ongoing processes and other users
    no_cpu=$(top -b -n 2 | grep 'Cpu(s)' | tail -1 | awk -F' ' '{print $2}' | sed 's/%us,//')
    no_processes=$(ps aux | grep dwifit | grep -v 'grep' | wc -l)
    # Get current hour
    current_H=$(date | tr -s " " | cut -d" " -f4 | tr -s ":" | cut -d":" -f1)
    echo "Current hour is "$current_H
    current_D=$(date | tr -s " " | cut -d" " -f3)
    echo "Current day is "$current_D
    # Print start-up conditions
    echo "---------------------------------------------------------------------"
    echo $no_processes" processes still running"
    echo "Next starting subfolder:"$spath
    cpu_condition=$(echo $no_cpu'<'$start_cpu | bc -l)
    if [ $cpu_condition -gt 0 ];
    then
        echo "Cpu load "$no_cpu"% is less than "$start_cpu
    else
        echo "Cpu load "$no_cpu"% is not less than "$start_cpu
    fi
    # Run only if all conditions match
	if [ $cpu_condition -gt 0 ];
	then
		echo $(date)" Starting in subfolder "$spath"/"
		echo $(date)" Starting in subfolder "$spath"/" >> run_con_CPU_log.txt
		execution=1
		for file in $(ls *ASCII.txt)
		do
        		resfile=$(echo ${file%.*}"_"$model"_results.txt" | sed -e 's/\.//g' | sed -e 's/resultstxt/results.txt/g')
	        	if [ -f $resfile ];
	        	then
                		echo "File "$resfile" exist already"
        		        continue
	        	fi
        		no_processing=$(ps aux | grep $file | grep $model | grep -v 'grep' | wc -l)
	        	if [ $no_processing -gt 0 ];
	        	then
                		echo "File "$resfile" is being processed already"
        		        continue
		        fi
        		execution=0
	        	echo $binary $file $model BFGS $settings_file
	        	$binary $file $model BFGS $settings_file
        		break
		done
		echo "Execution status:"$execution
		if [ $execution -ne 0 ];
		then
			break
		fi
	fi
	sleep 10
done
done
