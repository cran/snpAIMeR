# snpAIMeR

<!-- badges: start -->

<!-- badges: end -->

The goal of snpAIMeR is to assess the diagnostic power of SNP combinations using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within the [adegenet](https://github.com/thibautjombart/adegenet) R package.

Its value is in (1) identifying ancestry informative markers (AIMs) and (2) evaluating how well different marker combinations and can predict an unknown sample's population of origin.

The user provides candidate markers, SNP genotypes from individuals of known origin, a range of panel sizes, and a threshold value for an acceptable rate of correct sample identification.

snpAIMeR tests every marker combination within the specified minimum and maximum panel sizes. For each cross-validation replicate, individuals are randomly divided with 80% for the DAPC and 20% withheld as test samples. Results from the DAPC are used to predict the population of origin for each test individual, which is then compared with the known population label from the input file.

Because of the number of possible combinations, we recommend testing no more than 15 markers. For example, testing 15 markers in panel sizes of 1 to 15 (32,767 total combinations) with 1,000 cross-validation replicates on a system with 48 processor cores took about 5 hours and 20 GB RAM. To mitigate run time, snpAIMeR automatically uses n - 1 the number of available processor cores. Reducing the number of cross-validation replicates also reduces run time, however, we recommend no less than 100 replicates.

## Requirements

STRUCTURE formatted genotype file (.str or .stru). Individuals must have population assignments.

## Installation

You can install the released version of snpAIMeR like so:

```         
install.packages("snpAIMeR")
```

## Usage

```         
library(snpAIMeR)
```

### Run interactively

```         
snpAIMeR("interactive", verbose=TRUE)
```

Upon execution, the user is prompted with the following (do not quote paths):

```         
Enter path to working directory: 
Enter path to STRUCTURE file:
```

Then, the user is prompted (by adegenet) for information about the SNP genotype file:

```         
How many genotypes are there? 
How many markers are there? 
Which column contains labels for genotypes ('0' if absent)? 
Which column contains the population factor ('0' if absent)? 
Which other optional columns should be read (press 'return' when done)? 
Which row contains the marker names ('0' if absent)? 
Are genotypes coded by a single row (y/n)? 
```

Finally, after a few messages about the data (again from adegenet), the user is prompted for the following:

```         
Minimum number of markers in combination:
Maximum number of markers in combination:
Assignment rate threshold (minimum rate of successful assignments):
Number of cross-validation replicates:
```

### Run without interaction

```         
snpAIMeR("non-interactive", "/path/to/yaml/file", verbose=TRUE)
```

Non-interactive mode requires a config file in YAML format. Example [here](https://github.com/OksanaVe/snpAIMeR/blob/main/snpAIMeR_config.yml) and in the help documentation.

## Output

-   "Single_marker_assign_rate.pdf" has the mean correct assignment rate for each marker.
-   "All_combinations_assign_rate.csv" has the mean correct assignment rate for each combination tested (average of all cross-validation replicates).
-   "Panel_size_assign_rate.csv" and "Panel_size_assign_rate.pdf" have the mean correct assignment rate for each panel size tested (average of all combinations).
-   "Above_threshold_assign_rate.csv" lists the combinations with a mean correct assignment rate above the user-specified threshold.
