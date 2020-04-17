
To create twilight flats, you need to create night flats first, before running the scripts here.

Usage:
	In script folder:

	`$ ./make_twilight_flat.sh path_to_file/flatlist.txt`

	In flat files folder:

	`$ path_to_script/./make_twilight_flat.sh flatlist.txt`

flatlist.txt should contain the raw flat files names. Different filters and month can mix together.
This script run the four python script and generates the twilight flat flat_asu1_feb.fits (example name).

Each python script can be run individually in case you just need to run some part of the processes.

To view usage of each python script, simply execute the python scripts without arguments.


