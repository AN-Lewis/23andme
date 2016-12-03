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
   directory. The file *must* have a `.csv` extension.
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

### Options

* `--cases`
    * The directory that contains the linkage files of the cases.
    * Defaults to `./cases`.
* `--controls`
    * The directory that contains the linkage files of the controls.
    * Defaults to `./controls`.
* `-r`, `--recursive`
    * Search the case and control directories recursively.
* `--proband`
    * Set this to `case` if the proband is a case, `control` if the proband is a
      control, and `unknown` otherwise.
    * The script will prompt for this value if it is not provided on the command
      line.
* `-a`, `--alpha`
    * Exclude results with p-values less than this value.
    * Defaults to 0.05.
    * Pass `--alpha=1 --no-bonferroni` to see all results.
* `--no-bonferroni`
    * Do not divide alpha by 23 (the number of human chromosomes and therefore
      the number of independent/unlinked comparisons) to compensate for the
      [multiple comparison problem](https://en.wikipedia.org/wiki/Multiple_comparisons_problem)
      (a [Bonferroni correction](https://en.wikipedia.org/wiki/Bonferroni_correction)).
* `--method`
    * Set this to `chi` to compute p-values by
      [chi-squared test](https://en.wikipedia.org/wiki/Chi-squared_test),
      `fisher` for
      [Fisher exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test),
      `g` for the [G-test](https://en.wikipedia.org/wiki/G-test),
      and `auto` to use Fisher as long as it does not reduce the power. `auto`
      will choose Fisher if at least 3 cells in the 4-cell contingency table
      have expected values greater than or equal to 5 and no cell has an
      expected value less than 1, otherwise it will choose chi-squared.
    * Defaults to `auto`.
* `--no-yates`
    * If using the chi-squared test, do not apply the
      [Yates continuity correction](https://en.wikipedia.org/wiki/Yates%27s_correction_for_continuity).
      The Yates continuity correction increases the p value, and in theory does
      not reduce the power.
    * Pass `--method=chi --no-yates` for a simple z test of two proportions.
* `--misfits`
    * For each segment in the results, print the base filename of each file that
      did not fit its "case" or "control" category.
* `--randomize`
    * Ignore the folder structure and randomly assign each file to the case or
      control group. If the proband's affection is not specified with the
      `--proband=` option, the proband is also randomly assigned to the case or
      control group. This will give you an empirical impression of the false
      positive rate.
* `--help`
    * Print a synopsis of the available options.

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
