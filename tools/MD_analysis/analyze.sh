#!/bin/bash

box_x = 10
box_y = 10
box_z = 30

dt = 1000

n_frames = 1000
n_frames_rdf = 50

julia /home/salberti/MPMDP/tools/MD_analysis/MSD_args.jl "$n_frames" "$dt"
echo "MSD started."
date
julia /home/salberti/MPMDP/tools/MD_analysis/VACF_args.jl "$n_frames" "$dt"
echo "VACF started."
date
julia /home/salberti/MPMDP/tools/MD_analysis/RDF_args.jl "$n_frames_rdf" "$box_x" "$box_y" "$box_z"
echo "RDF started."
date
