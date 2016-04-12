#!/usr/bin/env python

import os
import time
import glob
import argparse

import sys

if __name__ == '__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser(description="""
Convert CACTUS wake node and grid data files from CSV to VTK using pyCactus and
pyCactusWake modules.

The wake node data is moved to '/WakeElemVTK' and the field data to '/FieldVTK'.
Recursive glob is used to find wake element data and field data matching the
default patterns:

   *_WakeData_*.csv
   *_FieldData_*.csv

However, alternate search patterns may also be specified using:

    --field_fnames_pattern 
    --wakeelem_fnames_pattern
""",formatter_class=argparse.RawDescriptionHelpFormatter)
                                    

    parser.add_argument("case_path",
                        help="path to CSV files to be converted.",
                        type=str)
    parser.add_argument("case_name",
                        help="prefix of output files, \
                        e.g. [case_name]_WakeData_1.csv",
                        type=str)
    parser.add_argument("output_path",
                        help="desired output path for VTK files.",
                        type=str)
    parser.add_argument("--delete",
                        help="delete original CSV files after converting.",
                        action="store_true",
                        default=False)
    parser.add_argument("--wakeelem_fnames_pattern",
                        help="glob pattern for wake element data files.",
                        type=str,
                        default=None)
    parser.add_argument("--field_fnames_pattern",
                        help="glob pattern for wake field data files.",
                        type=str,
                        default=None)
    parser.add_argument("--path_to_pyCactus",
                        help="path to pyCactus module directory (e.g.,\
                        '~/Repos/CACTUS-tools/'). This must be specified if \
                        pyCactus is not installed (is not available in the \
                        Python path)",
                        type=str,
                        default='')

    args = parser.parse_args()

    # read in command-line arguments
    case_path               = args.case_path
    case_name               = args.case_name
    output_path             = args.output_path
    wakeelem_fnames_pattern = args.wakeelem_fnames_pattern
    field_fnames_pattern    = args.field_fnames_pattern
    path_to_pyCactus        = args.path_to_pyCactus

    # append path
    sys.path.append(os.path.abspath(path_to_pyCactus))
    
    # import pyCactus
    import pyCactus.CactusRun

    # assemble a keyword args dictionary for the glob filename patterns
    # if they were entered
    kwargs = {}
    if wakeelem_fnames_pattern is not None:
        kwargs['wakeelem_fnames_pattern'] = wakeelem_fnames_pattern
    if field_fnames_pattern is not None:
        kwargs['field_fnames_pattern'] = field_fnames_pattern

    # load CactusRun instance
    run = pyCactus.CactusRun(case_path, case_name, **kwargs)
    
    # make output directories if they don't exist
    wakeelem_out_path = os.path.abspath(output_path + '/wakeelementVTK')
    field_out_path     = os.path.abspath(output_path + '/fieldVTK')

    # make directories
    if not os.path.exists(field_out_path):
        try:
            os.makedirs(field_out_path)
        except:
            print 'There was an error creating the directory' + field_out_path

    if not os.path.exists(wakeelem_out_path):
        try:
            os.makedirs(wakeelem_out_path)
        except:
            print 'There was an error creating the directory' + wakeelem_out_path

    # write VTK series
    try:
        run.field.write_vtk_series(field_out_path, name=case_name + 'Field')
    except:
        print 'Error: There was an error writing the field data out (possibly because it doesn\'t exist.)'
    else:
        # delete the CSV files if the --delete flag is specified
        if args.delete:
            for f in run.wakegrid_filenames:
                os.remove(f)
                print f + ' was deleted.'

    try:            
        run.wakeelems.write_vtk_series(wakeelem_out_path, name=case_name + 'WakeElem')
    except:
        print 'Error: There was an error writing the wake element data out (possibly because it doesn\'t exist.)'
    else:
        if args.delete:
        # delete the CSV files if the --delete flag is specified
            for f in run.wake_filenames:
                os.remove(f)
                print 'Deleted file: ' + f
