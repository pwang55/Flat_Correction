#!/bin/bash

path_file=$1


usage='
Usage:

	In script folder:
	$ ./make_twilight_flat.sh path_to_file/flatlist.txt

	In flat files folder:
	$ path_to_script/./make_twilight_flat.sh flatlist.txt

flatlist.txt should contain the raw flat files names. Different filters and month can mix together.
This script run the four python script and generates the twilight flat flat_asu1_feb.fits (example name)

'


# If no argument was given at all, show docs and end the script.
if [ -z "$path_file" ]
then
        echo "$usage"
        exit 1
fi


# Get the absolute directory path of this script so that it can find python scripts
script_dir=$(cd `dirname $0` && pwd)


python $script_dir/flat_sub_bias.py $path_file
python $script_dir/compare_biassub_flat2nights.py $path_file
python $script_dir/correct_flat_with_compare.py $path_file
python $script_dir/weighted_combine_twilightflat.py $path_file





