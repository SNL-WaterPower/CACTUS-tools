#!/usr/bin/env python
# pyCactusCsvToVtk.py
# =====================
# Converts CACTUS instantaneous wake node and wake grid data files from CSV format to VTK using
#   the pyCactus and pyCactusWake modules.
#
# The wake node data is moved to '/WakeNodeVTK' and the grid data in '/WakeGridVTK' in the 
# Expects that the wake node data and wake grid data match the patterns:
#
#   [case_name]_Wake_*.csv
#   [case_name]_WakeDef_*.csv
# 
#   Usage: pyCactusCsvToVtk.py [case_dir] [case_name] [output_dir]

import os
import sys
import time
import glob

import pyCactus
import pyCactusWake

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Not enough command line arguments.'
        print '    Usage:    pyCactusCsvToVtk.py [case_path] [case_name] [output_dir]'
        sys.exit()

    # read in command-line arguments
    case_path   = sys.argv[1]
    case_name   = sys.argv[2]
    output_path = sys.argv[3]

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
        
    try:            
        run.wakeelems.write_vtk_series(wake_node_out_path, name=case_name + 'WakeNode')
    except:
        print 'Error: There was an error writing the wake node data out (possibly because it doesn\'t exist.)'
