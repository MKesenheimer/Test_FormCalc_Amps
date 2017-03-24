#!/bin/bash

# does the script run on Linux or Mac?
UNAME=$(uname)
if [[ $UNAME == Darwin ]]; then
    # name of sed(gsed for Mac OSX)
    SED=gsed
else
    # name of sed (sed for Linux)
    SED=sed
fi

filename=$(basename "$1")
extension="${filename##*.}"
filename="${filename%.*}"

$SED -i -e "s/\t/      /g" $1
$SED -i -e "s/=        /=/g" $1
$SED -i -e "s/(        /(/g" $1
#mv $filename.f $filename.F
