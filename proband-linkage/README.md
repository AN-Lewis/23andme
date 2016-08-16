## 23andMe Proband Linkage Analyzer

This program compares shared segments reported by 23andMe to find segments
shared by a case group but not a control group or vice-versa.

To use this program:

1. Go to https://you.23andme.com/tools/relatives/ and request to share ancestry
   reports with your case relatives and control relatives.
2. Go to https://you.23andme.com/tools/relatives/dna/ and compare yourself to
   one other relative. This will only work if the relative has shared their
   ancestry report with you.
3. Copy and paste the tabular data into your favorite spreadsheet program.
4. Save the spreadsheet as a CSV in either the cases directory or the controls
   directory.
5. Repeat steps 1-4 until you have at least 10 cases and 10 controls. The more
   the better.
6. Install [SciPy](https://www.scipy.org/) by running `sudo pip install scipy`.
7. Execute `./proband-linkage.py`. The program will ask you if you, the proband,
   are a case or a control.

The program will output a tab-separated spreadsheet. If you are a case, the
output should contain segments that members of the case group share with you
significantly more often than members of the control group. If you are a
control, the output should contain segments that members of the case group share
with you significantly less often than members of the control group.

While this program can tell you where a clinically significant variant is
located in the genome, it cannot tell you what exactly the variant is.

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
