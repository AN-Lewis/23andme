## 23andMe raw data to PLINK file converter

This program takes as input one or more 23andMe raw data files and outputs:
* A [PLINK PED](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map) /
  [MERLIN PED](https://csg.sph.umich.edu/abecasis/Merlin/tour/input_files.html#relate) file
* A [PLINK MAP file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map).
* A [MERLIN MAP file](https://csg.sph.umich.edu/abecasis/Merlin/tour/input_files.html#mapfile)
* A [MERLIN DAT file](https://csg.sph.umich.edu/abecasis/Merlin/tour/input_files.html#pedfile)
These files can in turn be used with
[PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/),
[MERLIN](https://csg.sph.umich.edu/abecasis/Merlin/), and other bioinformatics
tools.

To use this program:

1. Execute `./get-hapmap.sh` to download HapMap Phase II data from NCBI. If you
   skip this step, the generated map file will not have genetic distances.
2. Download raw data files from 23andMe, unzip them, and put them in the cases,
   controls, and unknowns directories. Each file *must* have a `.txt` extension.
3. Execute `./convert.py`.

### Options
* `--cases`
    * The directory that contains the raw data files of the cases.
    * Defaults to `./cases`.
* `--controls`
    * The directory that contains the raw data files of the controls.
    * Defaults to `./controls`.
* `--unknowns`
    * The directory that contains the raw data files of the unknowns.
    * Defaults to `./unknowns`.
* `-r`, `--recursive`
    * Search the case and control directories recursively.
* `--family`
    * The family ID to put in the generated PED file, and the base filename of
      the generated PED and MAP files.
    * Defaults to `FAM001`.
* `--no-parents`
    * Set all parents to 0 in the PED file; do not attempt to infer them.
* `--no-sexes`
    * Set all sexes to 0 in the PED file; do not attempt to infer them.
    * Implies `--no-parents`.
* `--spacing`
    * The minimum genetic distance, in centimorgans, between markers in the
      output. Use this option to thin the data if it is too big for your linkage
      analysis program.
    * Defaults to `0`.
* `--out`
    * The directory to create the PED and MAP files in.
    * Defaults to the current directory.

### License

Copyright (c) 2016 Alex Henrie

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
