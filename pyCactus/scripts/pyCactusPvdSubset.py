#!/usr/bin/env python

import argparse
import xml.etree.cElementTree as ET


def pvd_time_subset(pvd_filename,pvd_filename_new,t_start=0,t_end=None):
    """Extracts a time subset from a ParaView collection file (*.pvd) in the
       range t_start to t_end.

       The timesteps are expected to be stored with the "timestep" attribute.
       The new collection file is written to the specified filename."""
    
    # load the XML file into a tree
    e = ET.parse(pvd_filename).getroot()
    
    # the collectiontree is the first node of the root tree
    collectiontree = e[0]

    # initialize empty list
    dataset_removelist = []

    # loop through collections
    for dataset in collectiontree:
        t = dataset.attrib['timestep']
        fname = dataset.attrib['file']

        # find which datasets meet t criterion
        if float(t) < t_start:
            dataset_removelist.append(dataset)
        
        if t_end != None:
            if float(t) > t_end:
                dataset_removelist.append(dataset)

    # remove the datasets
    for dataset in dataset_removelist:
        collectiontree.remove(dataset)
    
    # write
    e = ET.ElementTree(e)
    e.write(pvd_filename_new, xml_declaration=True)


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description="""
Extracts a time subset from a ParaView collection file (.pvd) given a range of
timesteps and writes it to a new .pvd file.""",formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("pvd_filename", help="path to *.pvd", type=str)
    parser.add_argument("pvd_filename_new", help="path to new *.pvd file", type=str)
    parser.add_argument("--t_start", help="time before which to clip *.pvd file.", type=float)
    parser.add_argument("--t_end", help="time after which to clip *.pvd file.", type=float)
    
    args = parser.parse_args()

    # read in command-line arguments
    pvd_filename     = args.pvd_filename
    pvd_filename_new = args.pvd_filename_new
    t_start          = args.t_start
    t_end            = args.t_end


    # call the converter function
    pvd_time_subset(pvd_filename,pvd_filename_new,t_start,t_end)

    return parser

if __name__ == "__main__":
    main()
