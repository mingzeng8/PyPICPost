# PyPICPost
python3 scripts for particle-in-cell (PIC) code post-processing
Written in python3. Using "property" which python2 does not support.
Supporting OSIRIS, HiPACE (not OpenPMD) and QuickPIC outputs.
"outfile.py" is the main script which one can use it to get some plots directly.
"W_vs_t.py" is the script for ploting the time evolution of the laser size.
"FrameMovie.py" is for plotting a series of snapshots with driver and plasma densities.

Optional packages:
parse, for detecting the available file numbers of the PIC output. To install with pip:
  pip install parse
or on a computer without root access:
  pip install --user parse

openpmd-viewer, for code_name being 'fbpic'. To install with anaconda:
  conda install -c conda-forge openpmd-viewer
