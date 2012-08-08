#!/bin/bash

#This script should loop through a set of .asc files and generate z-f
#profile files. It should then append an identifier to the input file
#when specifying the filename of the output file This script assumes
#that the tabulated data is ordered t v i, otherwise the algorithm
#will break!

for file in *.asc
do
    ls "$file"  # Lists all files in $PWD (current directory).
    # generate z-f profile for file and send the output to folder z-f
    ./zfgenerator -f "$file"
done
