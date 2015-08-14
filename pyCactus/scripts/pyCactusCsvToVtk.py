#!/usr/bin/env python
# pyCactusCsvToVtk.py
# =====================
# Converts CACTUS instantaneous wake node and wake grid data files from CSV format to VTK using
# the pyCactus and pyCactusWake modules.
#
# The wake node data is moved to '/WakeNodeVTK' and the grid data in '/WakeGridVTK' in the 
# Expects that the wake node data and wake grid data match the patterns:
#
#    [case_name]_Wake_*.csv
#    [case_name]_WakeDef_*.csv"""
#
# Usage: pyCactusCsvToVtk.py [case_path] [case_name] [output_path]

import os
import time
import glob
import argparse

import pyCactus
import pyCactusWake

if __name__ == '__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser(description="""Convert CACTUS wake node and grid data files from CSV to VTK using pyCactus and pyCactusWake modules.

The wake node data is moved to '/WakeNodeVTK' and the grid data in '/WakeGridVTK' in the 
Expects that the wake node data and wake grid data match the patterns:

   [case_name]_Wake_*.csv
   [case_name]_WakeDef_*.csv""")

    parser = argparse.ArgumentParser()
    parser.add_argument("case_path", help="path to CSV files to be converted.", type=str)
    parser.add_argument("case_name", help="prefix of output files, e.g. [case_name]_WakeData_1.csv", type=str)
    parser.add_argument("output_path", help="desired output path for VTK files.", type=str)
    parser.add_argument("--delete", help="delete original CSV files after converting.", action="store_true", default=False)
    
    args = parser.parse_args()

    # read in command-line arguments
    case_path   = args.case_path
    case_name   = args.case_name
    output_path = args.output_path

    # guess the geom filename
    geom_fname, = glob.glob('*.geom')

    # load CactusRun instance
    run = pyCactus.CactusRun(case_path, case_name, geom_fname)
    
    # make output directories if they don't exist
    wake_node_out_path = os.path.abspath(output_path + '/WakeNodeVTK')
    wake_grid_out_path = os.path.abspath(output_path + '/WakeGridVTK')

    # make directories
    if not os.path.exists(wake_grid_out_path):
        try:
            os.makedirs(wake_grid_out_path)
        except:
            print 'There was an error creating the directory' + wake_grid_out_path

    if not os.path.exists(wake_node_out_path):
        try:
            os.makedirs(wake_node_out_path)
        except:
            print 'There was an error creating the directory' + wake_node_out_path

    # write VTK series
    try:
        run.wakegrid.write_vtk_series(wake_grid_out_path, name=case_name + 'WakeGrid')
    except:
        print 'Error: There was an error writing the wake grid data out (possibly because it doesn\'t exist.)'
    else:
        # delete the CSV files if the --delete flag is specified
        if args.delete:
            for f in run.wakegrid_filenames:
                os.remove(f)
                print f + ' was deleted.'

    try:            
        run.wakeelems.write_vtk_series(wake_node_out_path, name=case_name + 'WakeNode')
    except:
        print 'Error: There was an error writing the wake node data out (possibly because it doesn\'t exist.)'
    else:
        if args.delete:
        # delete the CSV files if the --delete flag is specified
            for f in run.wake_filenames:
                os.remove(f)
                print f + ' was deleted.'
