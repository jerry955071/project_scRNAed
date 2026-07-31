#!/bin/bash
# Usage:
# ./bgzip-tabix.sh {input} {output.gz} {output.tbi}

# compress with bgzip
bgzip -o $2 $1

# index with tabix
tabix -s 1 -b 2 -e 2 -c 'R' $2

# check if the index file was created
if [ -f $3 ]; then
    echo "Index file created successfully: $3"
else
    echo "Failed to create index file: $3"
fi