# split_file_by_time.py
# =====================
# Splits a CACTUS time-series data file into multiple files, with each file containing a single
#   timestep. Assumes that the time data is in the first column of the file.
#
#   Names the files as datafileX.csv, where X is the timestep number. The timestep number is
#   calculated assuming a constant timestep. The optional parameters ts_start and ts_interval
#   specify the starting timestep and number of timesteps between data writes.
#
#   Usage: split_file_by_time.py [datafile.csv] [output_dir] [ts_start] [ts_interval]
#

import os
import sys
import time as pytime

def split_file_by_time(data_filename, output_dir, ts_start, ts_interval):
    # split the filename
    path, fname = os.path.split(os.path.abspath(data_filename))
    basename, ext = fname.split('.')

    # initialize values
    lines     = ''
    time_prev = ''
    ts_count = 0

    # make output directory if necessary
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # loop through the file
    with open(path + '/' + basename + '.' + ext) as f:
        # get current time
        etime_start = pytime.time()
        etime_file = etime_start
        for i, line in enumerate(f):
            if i == 0:
                header = line

            else:
                time = line.split(',')[0]
                if time_prev == '' or time == time_prev:
                    # add the current line to the current timestep string
                    lines = lines + line

                elif time != time_prev:
                    # at a new timestep, write the header and lines to a file
                    outfilename = output_dir + '/' + basename + '_' + str(ts_start + ts_count*ts_interval) + '.' + ext
                    out = open(outfilename, 'w')
                    out.write(header)
                    out.write(lines)
                    out.close()

                    print 'Finished writing file: ' + outfilename + ' in %2.2f seconds' % (pytime.time() - etime_file)
                    etime_file = pytime.time()

                    # make the current line the starting line
                    lines = line

                    # increment the timestep counter
                    ts_count = ts_count + 1

                time_prev = time

        # write the final timestep
        outfilename = output_dir + '/' + basename + '_' + str(ts_start + ts_count*ts_interval) + '.' + ext
        out = open(outfilename, 'w')
        out.write(header)
        out.write(lines)
        out.close()
        print 'Finished writing file: ' + outfilename + ' in %2.2f seconds' % (pytime.time() - etime_file)
        print 'Completed in %2.2f seconds' % (pytime.time() - etime_start)

if __name__ == '__main__':
    # check for command line args
    if len(sys.argv) < 3:
        print 'Error -- need to specify input files.'
        print '\tUsage: split_file_by_time.py [datafile.csv] [output_dir] [ts_start] [ts_interval]'
        sys.exit()

    # read in command line arguments
    data_filename       = sys.argv[1]
    output_dir          = sys.argv[2]

    # read in optional parameters; if they don't exist, use defaults as shown
    try: ts_start       = sys.argv[3]
    except: ts_start    = 1

    try: ts_interval    = sys.argv[4]
    except: ts_interval = 1

    split_file_by_time(data_filename, output_dir, ts_start, ts_interval)