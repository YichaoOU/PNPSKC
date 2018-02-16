# PNPSKC
### a Greedy algorithm for solving the Positive Negative Partial Set K-Cover problem

## Getting Started

The Greedy_PNPSKC algorithm can be used to do feature selection in a sense that multiple features are expected to co-occur in the positive samples rather than in the negative samples.

```
Git clone https://github.com/unfashionable/PNPSKC.git

python Greedy_PNPSKC.py HRGP_dataset/HRGP_raw_both_strand.csv 3
```

You should be able to see a result file `Greedy.3.result`.

### Parameters

The first parameter is the input file --- csv format. The first column is used for index. Using the HRGP dataset as an example, the first column is a list of sequence IDs, the header row is a list of motif names. If a motif j occurs in sequence i, then cell (i,j) = 1; otherwise, it is 0. The last column is the label; 1 indicates positive, other numbers indicates negative, e.g. 0 or -1.

The second parameter is the `K`. K is the number of times that a sample is expected to be covered, which is denoted as k-covered. In other words, K is the number of features that are expected to co-occur in a positive sample. 

### Prerequisites
The following python libs are required:

`Pandas`

`sklearn`

`numpy`

## Authors

Yichao Li, Allan Showalter, and Lonnie Welch

A Novel Discriminative Set Multi-Cover Model for Discovering DNA Motifs and Motif Pairs

