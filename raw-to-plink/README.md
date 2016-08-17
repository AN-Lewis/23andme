## 23andMe raw data to MAP file converter

This program takes as input one or more 23andMe raw data files and outputs a
single [PED file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map) and
a single [MAP file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map).
These files can in turn be used with
[PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) or other bioinformatics
tools.

To use this program:

1. Download raw data files from 23andMe, unzip them, and put them in the cases,
   controls, and unknowns directories.
2. Execute `./raw-to-plink`.

## License

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
