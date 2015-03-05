

cDNA yield after C1 runs
========================

Summary
-------

This page shows where to retrieve the [data files](#datasets) containing the
measurement of DNA concentration in the 96-well plates after collection of the
cDNAs from the capture array.  It then presents some [commands](#functions) to
import the data in `R`, and uses them to produce a consolidated
[table](#output-format) for the five runs, that is [saved](#save-data) for
integrated analysis later.  Finally, some [quality controls](#qc) are ran,
showing strong variation between runs for the DNA yield after the on-chip PCR.
An inspection of the [standard curves](#standard-curves) rules out that it
could be a simple problem of measurement.


Method
------

After a C1 run, the products are recovered from the outlets and transferred to a
96-well plate.  Following the standard procedure (PN 100-5950 A1, page 25),
2 μL are quantified on a 384-well fluorometer using the PicoGreen dye.  See
also _PicoGreen Quantitation of DNA: Effective Evaluation of Samples Pre-or
Post-PCR._ Susan J. Ahn, José Costa and Janet Rettig Emanuel, _Nucl. Acids
Res._ (1996) [24(13):2623-2625](http://dx.doi.org/10.1093/nar/24.13.2623).

The measured fluorescence intensities are transferred in an Excel sheet
provided by Fluidigm (PN 100-6160), which converts them to DNA concentrations
by linear regression to a standard curve.

The `R` commands here load the Excel sheets and produce one consolidated table
for all runs, on which a few quality controls are run.


<a name='datasets'>Datasets</a>
-------------------------------

Each file is named following the serial ID of the C1 chip from which DNA was collected.

The files are available from a [temporary location](BASEURL=https://briefcase.riken.jp/proself/publicweb/publicweb.go/RCkIQA6ZAohAcycBOlFL0_ZSt5-IUtPtwtvH_w-vUo79) and will be deposited in a proper place
later.  The files need to be downloaded before running this knitr script.


```bash
shasum -c << __CHECKSUMS__
6d58020277bd2eb42f5702eabe17194671d9e87f  1772-062-248.picogreen.xlsx
83844dd8e077d54a619e6bff396992b23d8a26e8  1772-062-249.picogreen.xlsx
656fb4bd93466cd4f005881a046f3f17c17f4a78  1772-064-103.picogreen.xlsx
f4e2064d7c09826a8ad14ed2a95a66255018da39  1772-067-038.picogreen.xlsx
c9a1ece68f037589d40b382f59d202dbb90425de  1772-067-039.picogreen.xlsx
__CHECKSUMS__
```

```
## shasum: 1772-062-248.picogreen.xlsx: 
## 1772-062-248.picogreen.xlsx: FAILED open or read
## shasum: 1772-062-249.picogreen.xlsx: No such file or directory
## 1772-062-249.picogreen.xlsx: FAILED open or read
## shasum: 1772-064-103.picogreen.xlsx: No such file or directory
## 1772-064-103.picogreen.xlsx: FAILED open or read
## shasum: 1772-067-038.picogreen.xlsx: No such file or directory
## 1772-067-038.picogreen.xlsx: FAILED open or read
## shasum: 1772-067-039.picogreen.xlsx: No such file or directory
## 1772-067-039.picogreen.xlsx: FAILED open or read
## shasum: WARNING: 5 listed files could not be read
```


<a name='functions'>Functions to load the data in `R`</a>
---------------------------------------------------------


```r
library(gdata)
```

```
## Error in library(gdata): there is no package called 'gdata'
```

```r
library(ggplot2)
library(reshape)
```

```
## Error in library(reshape): there is no package called 'reshape'
```

These files in _Excel_ format contain the final result in their third sheet,
from line 42 to 49.

First, extract the concentrations as a 8 × 12 table.  An example is shown for
run `1772-062-248`.


```r
readConcentrationTable <- function (FILE) {
  picogreen <- read.xls( FILE
                       , sheet =  3
                       , skip  = 41
                       , nrow  =  8
                       , head  = FALSE
                       , blank.lines.skip = FALSE )[,1:13]
  colnames(picogreen) <- c('row', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12')
  picogreen
}
readConcentrationTable('1772-062-248.picogreen.xlsx')
```

```
## Error in readConcentrationTable("1772-062-248.picogreen.xlsx"): could not find function "read.xls"
```

Then, transform this table to have one measurement by line, with the run ID,
and the row and column name.  This is a typical format when plotting data with
_[ggplot2](http://ggplot2.org/)_.


```r
meltConcentrationTable <- function (RUN, TABLE) {
  picogreen <- melt(TABLE, id.vars='row')
  colnames(picogreen) <- c('row', 'column', 'concentration')
  picogreen[,"run"] <- RUN
  picogreen[,"well"] <- paste(picogreen$row, picogreen$column, sep='')
  picogreen <- picogreen[, c('run', 'well', 'row', 'column', 'concentration')]
  picogreen
}
```

The function below outputs the data for one run, provided that a properly named
file (`run id` plus `.picogreen.xlsx`) is available in the same directory.


```r
read_pg <- function(RUN) {
  FILE <- paste(RUN, 'picogreen.xlsx', sep='.')
  picogreen <- readConcentrationTable(FILE)
  picogreen <- meltConcentrationTable(RUN, picogreen)
  picogreen
}
head(read_pg('1772-062-248'))
```

```
## Error in readConcentrationTable(FILE): could not find function "read.xls"
```

The standard curve is also in sheet 3, from row 2 to 12.  Here is a function to
get the standard curve from one run (with background correction already
applied).  The output is self-explanatory.


```r
read_sc <- function(RUN) {
  FILE <- paste(RUN, "picogreen.xlsx", sep = ".")
  sc <- read.xls(FILE, sheet=3, skip=2, nrow=10, header=FALSE)[,c(2,5)]
  sc <- cbind(RUN, sc)
  colnames(sc) <- c('run', 'dna', 'fluorescence')
  sc$fluorescence <- as.numeric(as.character(sc$fluorescence))
  return(sc)
}
read_sc('1772-062-248')
```

```
## Error in read_sc("1772-062-248"): could not find function "read.xls"
```

<a name="output-format">Consolidated file (format)</a>
------------------------------------------------------

The file `cDNA_concentration.csv` is made from the files above using `R`, and has
the following columns.

### `run`

The serial ID of the C1 chip for a given run.  Example: `1772-062-248`.

### `well`

The coordinates in the 96-well plate where the cDNAs have been transferred at
the end of the C1 run. Examples: `A01`, `F08`, `C12`, etc.  Combined with the
run ID, this uniquely identifies a cell.

### `row`

The row coordinates in the 96-well plate where the cDNAs have been transferred
at the end of the run.  Possible values: `A`, `B`, `C`, `D`, `E`, `F`, `G` and
`H`.

### `column`

The column coordinates in the 96-well plate where the cDNAs have been
transferred at the end of the run.  Possible values: `01`, `02`, `03`, `04`,
`05`, `06`, `07`, `08`, `09`, `10`, `11` and `12`.

### `concentration`

The DNA concentration, in ng/μL.


<a name="save-data">Consolidated file (preparation)</a>
-------------------------------------------------------

Each Excel sheet is loaded in `R`, the DNA concentrations are extracted, and
added to a table saved under the name `cDNA_concentration.csv`.

The `cDNA_concentration.csv` file contains the data for the following runs.


```r
RUNS <- c('1772-062-248', '1772-062-249', '1772-064-103', '1772-067-038', '1772-067-039')
```

Create a `picogreen` table for the first run being processed, append the other
runs, and save the file.


```r
for (RUN in RUNS) {
  if (! exists('picogreen'))
    {picogreen <- read_pg(RUN)}
  else
    {picogreen <- rbind(picogreen, read_pg(RUN))}
}
```

```
## Error in readConcentrationTable(FILE): could not find function "read.xls"
```

```r
summary(picogreen)
```

```
## Error in summary(picogreen): object 'picogreen' not found
```

```r
write.csv(file='cDNA_concentration.csv', picogreen, row.names=FALSE)
```

```
## Error in is.data.frame(x): object 'picogreen' not found
```

<a name='qc'>Quality control</a>
--------------------------------

There is a strong variation between runs.


```r
qplot( data=picogreen
     , run
     , concentration
     , geom="boxplot",
     , colour=run) + coord_flip()
```

```
## Error in ggplot(data, aesthetics, environment = env): object 'picogreen' not found
```

Still, for most runs except _1772-064-103_, it is possible to detect
low-concentration outliers, were there probably was no cell in the chamber.
Note that the scale of each histogram is different.


```r
qplot( data=picogreen
     , concentration
     , geom="histogram"
     , colour=run) + facet_wrap( ~run, scales='free')
```

```
## Error in ggplot(data, aesthetics, environment = env): object 'picogreen' not found
```

<a name='standard-curves'>Comparison between standard curves</a>
----------------------------------------------------------------

Load the data from all runs.


```r
for (RUN in RUNS) {
    if (!exists("sc")) {
        sc <- read_sc(RUN)
    } else {
        sc <- rbind(sc, read_sc(RUN))
    }
}
```

```
## Error in read_sc(RUN): could not find function "read.xls"
```

While the DNA concentrations in `1772-062-248` might have been overestimated
as suggested by the shifted standard curve, overall the calibration of the
fluorometer is stable and does not explain the strong variations of the DNA
yield. 


```r
ggplot(
  sc,
  aes(
    x=fluorescence,
    y=dna,
    colour=run)
) + geom_point() +
    geom_line() + 
    scale_x_log10('Average fluorescence (background subtracted)') +
    scale_y_log10('DNA concentration (ng/μL)')
```

```
## Error in ggplot(sc, aes(x = fluorescence, y = dna, colour = run)): object 'sc' not found
```

Session Info
------------


```r
sessionInfo()
```

```
## R version 2.15.3 (2013-03-01)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=C                 LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_0.9.3.1 knitr_1.9      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-2   dichromat_2.0-0    digest_0.6.3       evaluate_0.5.5     formatR_1.0        grid_2.15.3       
##  [7] gtable_0.1.2       labeling_0.2       MASS_7.3-23        munsell_0.4.2      plyr_1.8           proto_0.3-10      
## [13] RColorBrewer_1.0-5 reshape2_1.2.2     scales_0.2.3       stringr_0.6.2      tools_2.15.3
```
