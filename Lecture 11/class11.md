---
title: "Bioinformatics Class 11"
author: "Peng Fei Huang"
date: "5/9/2018"
output: 
  html_document: 
    keep_md: yes
---




## PDB Statistics 

Import out PDB statistics CSV file and calculate precent strucutres by experimental method 

```r
p <- read.csv("Data Export Summary.csv", row.names = 1)
```



```r
percent <- (p$Total/sum(p$Total)) *100
names(percent) <-  row.names(p)
percent
```

```
##               X-Ray                 NMR Electron Microscopy 
##         89.51673340          8.71321614          1.51239392 
##               Other        Multi Method 
##          0.16986775          0.08778879
```


## Using Bio3D

Load the bio3d package



```r
library(bio3d)
```

Read in our HIV-Pr structure 


```r
pdb <- read.pdb("1hsg")
```

```
##   Note: Accessing on-line PDB file
```

```r
pdb
```

```
## 
##  Call:  read.pdb(file = "1hsg")
## 
##    Total Models#: 1
##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 172  (residues: 128)
##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, xyz, seqres, helix, sheet,
##         calpha, remark, call
```


```r
attributes(pdb)
```

```
## $names
## [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  
## 
## $class
## [1] "pdb" "sse"
```


```r
head(pdb$atom)
```

```
##   type eleno elety  alt resid chain resno insert      x      y     z o
## 1 ATOM     1     N <NA>   PRO     A     1   <NA> 29.361 39.686 5.862 1
## 2 ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
## 3 ATOM     3     C <NA>   PRO     A     1   <NA> 29.760 38.071 4.022 1
## 4 ATOM     4     O <NA>   PRO     A     1   <NA> 28.600 38.302 3.676 1
## 5 ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508 37.541 6.342 1
## 6 ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296 37.591 7.162 1
##       b segid elesy charge
## 1 38.10  <NA>     N   <NA>
## 2 40.62  <NA>     C   <NA>
## 3 42.64  <NA>     C   <NA>
## 4 43.40  <NA>     O   <NA>
## 5 37.87  <NA>     C   <NA>
## 6 38.40  <NA>     C   <NA>
```

# Print a subset of $atom data for the first two atoms

```r
pdb$atom[1:2, c("eleno", "elety", "x","y","z")]
```

```
##   eleno elety      x      y     z
## 1     1     N 29.361 39.686 5.862
## 2     2    CA 30.307 38.663 5.319
```

# Note that individual $atom records can also be accessed like this

```r
pdb$atom$elety[1:2]
```

```
## [1] "N"  "CA"
```

# Which allows us to do the following

```r
#plot.bio3d(pdb$atom$b[pdb$calpha], sse=pdb, typ="l", ylab=“B-factor”)
```

# Print a summary of the coordinate data in $xyz

```r
pdb$xyz
```

```
## 
##    Total Frames#: 1
##    Total XYZs#:   5058,  (Atoms#:  1686)
## 
##     [1]  29.361  39.686  5.862  <...>  30.112  17.912  -4.791  [5058] 
## 
## + attr: Matrix DIM = 1 x 5058
```

# Examine the row and column dimensions

```r
dim(pdb$xyz)
```

```
## [1]    1 5058
```

# Print coordinates for the first two atom

```r
pdb$xyz[ 1, atom2xyz(1:2) ]
```

```
## [1] 29.361 39.686  5.862 30.307 38.663  5.319
```

# Select all C-alpha atoms (return their indices)

```r
ca.inds <- atom.select(pdb, "calpha")
ca.inds
```

```
## 
##  Call:  atom.select.pdb(pdb = pdb, string = "calpha")
## 
##    Atom Indices#: 198  ($atom)
##    XYZ  Indices#: 594  ($xyz)
## 
## + attr: atom, xyz, call
```


# Print details of the first few selected atoms

```r
head( pdb$atom[ca.inds$atom, ] )
```

```
##    type eleno elety  alt resid chain resno insert      x      y     z o
## 2  ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
## 9  ATOM     9    CA <NA>   GLN     A     2   <NA> 30.158 36.492 2.199 1
## 18 ATOM    18    CA <NA>   ILE     A     3   <NA> 29.123 33.098 3.397 1
## 26 ATOM    26    CA <NA>   THR     A     4   <NA> 29.774 30.143 1.062 1
## 33 ATOM    33    CA <NA>   LEU     A     5   <NA> 27.644 27.003 1.144 1
## 41 ATOM    41    CA <NA>   TRP     A     6   <NA> 30.177 24.150 1.279 1
##        b segid elesy charge
## 2  40.62  <NA>     C   <NA>
## 9  41.30  <NA>     C   <NA>
## 18 34.13  <NA>     C   <NA>
## 26 30.14  <NA>     C   <NA>
## 33 30.12  <NA>     C   <NA>
## 41 30.82  <NA>     C   <NA>
```


# And selected xyz coordinates

```r
head( pdb$xyz[, ca.inds$xyz] )
```

```
## [1] 30.307 38.663  5.319 30.158 36.492  2.199
```


# Select chain A

```r
a.inds <- atom.select(pdb, chain="A")
```

# Select C-alphas of chain A

```r
ca.inds <- atom.select(pdb, "calpha", chain="A")
```

# We can combine multiple selection criteria to return their intersection

```r
cab.inds <- atom.select(pdb, elety=c("CA","CB"), chain="A", resno=10:20)
```




#Q8. Use the Bio3D write.pdb() function to write out a protein only PDB file for viewing in VMD. Also write out a second separate PDB file for the ligand with residue name MK1


```r
inds.ligand <- atom.select(pdb, "ligand")
inds.protein <- atom.select(pdb, "protein")
inds.protein
```

```
## 
##  Call:  atom.select.pdb(pdb = pdb, string = "protein")
## 
##    Atom Indices#: 1514  ($atom)
##    XYZ  Indices#: 4542  ($xyz)
## 
## + attr: atom, xyz, call
```


Check we have what we want

```r
pdb$atom[inds.ligand$atom,]
```

```
##        type eleno elety  alt resid chain resno insert      x      y      z
## 1515 HETATM  1517    N1 <NA>   MK1     B   902   <NA>  9.280 23.763  3.004
## 1516 HETATM  1518    C1 <NA>   MK1     B   902   <NA>  9.498 23.983  4.459
## 1517 HETATM  1519    C2 <NA>   MK1     B   902   <NA> 10.591 24.905  4.962
## 1518 HETATM  1520    C3 <NA>   MK1     B   902   <NA> 10.591 24.864  6.466
## 1519 HETATM  1521    O1 <NA>   MK1     B   902   <NA> 10.937 23.849  7.057
## 1520 HETATM  1522    N2 <NA>   MK1     B   902   <NA> 10.193 25.953  7.094
## 1521 HETATM  1523    C4 <NA>   MK1     B   902   <NA> 10.145 26.250  8.490
## 1522 HETATM  1524    C5 <NA>   MK1     B   902   <NA>  9.379 27.577  8.641
## 1523 HETATM  1525    C6 <NA>   MK1     B   902   <NA> 11.398 26.347  9.074
## 1524 HETATM  1526    C7 <NA>   MK1     B   902   <NA>  9.364 25.283  9.268
## 1525 HETATM  1527    N3 <NA>   MK1     B   902   <NA> 11.819 24.282  4.355
## 1526 HETATM  1528    C8 <NA>   MK1     B   902   <NA> 11.753 23.776  2.961
## 1527 HETATM  1529    C9 <NA>   MK1     B   902   <NA> 10.440 23.182  2.493
## 1528 HETATM  1530   C10 <NA>   MK1     B   902   <NA> 13.083 24.963  4.552
## 1529 HETATM  1531   C11 <NA>   MK1     B   902   <NA> 14.203 24.064  5.078
## 1530 HETATM  1532    O2 <NA>   MK1     B   902   <NA> 15.242 24.884  4.634
## 1531 HETATM  1533   C12 <NA>   MK1     B   902   <NA> 14.440 23.761  6.569
## 1532 HETATM  1534   C13 <NA>   MK1     B   902   <NA> 15.573 22.821  7.005
## 1533 HETATM  1535   C14 <NA>   MK1     B   902   <NA> 15.644 22.664  8.534
## 1534 HETATM  1536   C15 <NA>   MK1     B   902   <NA> 16.733 21.750  8.961
## 1535 HETATM  1537   C16 <NA>   MK1     B   902   <NA> 18.058 21.916  8.553
## 1536 HETATM  1538   C17 <NA>   MK1     B   902   <NA> 19.037 21.016  8.947
## 1537 HETATM  1539   C18 <NA>   MK1     B   902   <NA> 18.673 19.939  9.758
## 1538 HETATM  1540   C19 <NA>   MK1     B   902   <NA> 17.347 19.773 10.176
## 1539 HETATM  1541   C20 <NA>   MK1     B   902   <NA> 16.374 20.687  9.772
## 1540 HETATM  1542   C21 <NA>   MK1     B   902   <NA> 15.447 21.440  6.373
## 1541 HETATM  1543    O3 <NA>   MK1     B   902   <NA> 14.367 20.831  6.397
## 1542 HETATM  1544    N4 <NA>   MK1     B   902   <NA> 16.583 20.913  5.924
## 1543 HETATM  1545   C22 <NA>   MK1     B   902   <NA> 16.692 19.500  5.604
## 1544 HETATM  1546   C23 <NA>   MK1     B   902   <NA> 18.067 18.945  5.936
## 1545 HETATM  1547    O4 <NA>   MK1     B   902   <NA> 19.061 19.938  5.729
## 1546 HETATM  1548   C24 <NA>   MK1     B   902   <NA> 18.226 17.726  5.057
## 1547 HETATM  1549   C25 <NA>   MK1     B   902   <NA> 17.476 17.904  3.760
## 1548 HETATM  1550   C26 <NA>   MK1     B   902   <NA> 17.500 17.363  2.496
## 1549 HETATM  1551   C27 <NA>   MK1     B   902   <NA> 16.613 17.872  1.541
## 1550 HETATM  1552   C28 <NA>   MK1     B   902   <NA> 15.722 18.906  1.865
## 1551 HETATM  1553   C29 <NA>   MK1     B   902   <NA> 15.683 19.479  3.129
## 1552 HETATM  1554   C30 <NA>   MK1     B   902   <NA> 16.504 19.061  4.128
## 1553 HETATM  1555   C31 <NA>   MK1     B   902   <NA>  8.033 23.100  2.604
## 1554 HETATM  1556   C32 <NA>   MK1     B   902   <NA>  6.666 23.739  2.876
## 1555 HETATM  1557   C33 <NA>   MK1     B   902   <NA>  6.158 24.808  2.124
## 1556 HETATM  1558    N5 <NA>   MK1     B   902   <NA>  4.911 25.430  2.300
## 1557 HETATM  1559   C34 <NA>   MK1     B   902   <NA>  4.207 24.839  3.348
## 1558 HETATM  1560   C35 <NA>   MK1     B   902   <NA>  4.654 23.774  4.136
## 1559 HETATM  1561   C36 <NA>   MK1     B   902   <NA>  5.905 23.211  3.897
##      o     b segid elesy charge
## 1515 1 28.25  <NA>     N   <NA>
## 1516 1 30.30  <NA>     C   <NA>
## 1517 1 27.27  <NA>     C   <NA>
## 1518 1 28.85  <NA>     C   <NA>
## 1519 1 29.59  <NA>     O   <NA>
## 1520 1 22.29  <NA>     N   <NA>
## 1521 1 23.47  <NA>     C   <NA>
## 1522 1 27.66  <NA>     C   <NA>
## 1523 1 21.71  <NA>     C   <NA>
## 1524 1 22.75  <NA>     C   <NA>
## 1525 1 28.91  <NA>     N   <NA>
## 1526 1 26.24  <NA>     C   <NA>
## 1527 1 27.47  <NA>     C   <NA>
## 1528 1 20.86  <NA>     C   <NA>
## 1529 1 21.68  <NA>     C   <NA>
## 1530 1 15.87  <NA>     O   <NA>
## 1531 1 21.49  <NA>     C   <NA>
## 1532 1 26.89  <NA>     C   <NA>
## 1533 1 28.67  <NA>     C   <NA>
## 1534 1 26.89  <NA>     C   <NA>
## 1535 1 29.22  <NA>     C   <NA>
## 1536 1 29.22  <NA>     C   <NA>
## 1537 1 30.97  <NA>     C   <NA>
## 1538 1 29.25  <NA>     C   <NA>
## 1539 1 29.96  <NA>     C   <NA>
## 1540 1 29.35  <NA>     C   <NA>
## 1541 1 32.66  <NA>     O   <NA>
## 1542 1 31.19  <NA>     N   <NA>
## 1543 1 29.22  <NA>     C   <NA>
## 1544 1 28.82  <NA>     C   <NA>
## 1545 1 28.32  <NA>     O   <NA>
## 1546 1 32.05  <NA>     C   <NA>
## 1547 1 31.29  <NA>     C   <NA>
## 1548 1 32.00  <NA>     C   <NA>
## 1549 1 28.00  <NA>     C   <NA>
## 1550 1 29.01  <NA>     C   <NA>
## 1551 1 27.70  <NA>     C   <NA>
## 1552 1 31.86  <NA>     C   <NA>
## 1553 1 36.25  <NA>     C   <NA>
## 1554 1 42.75  <NA>     C   <NA>
## 1555 1 47.41  <NA>     C   <NA>
## 1556 1 51.38  <NA>     N   <NA>
## 1557 1 50.60  <NA>     C   <NA>
## 1558 1 49.34  <NA>     C   <NA>
## 1559 1 44.71  <NA>     C   <NA>
```


```r
head(pdb$atom[inds.protein$atom,])
```

```
##   type eleno elety  alt resid chain resno insert      x      y     z o
## 1 ATOM     1     N <NA>   PRO     A     1   <NA> 29.361 39.686 5.862 1
## 2 ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
## 3 ATOM     3     C <NA>   PRO     A     1   <NA> 29.760 38.071 4.022 1
## 4 ATOM     4     O <NA>   PRO     A     1   <NA> 28.600 38.302 3.676 1
## 5 ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508 37.541 6.342 1
## 6 ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296 37.591 7.162 1
##       b segid elesy charge
## 1 38.10  <NA>     N   <NA>
## 2 40.62  <NA>     C   <NA>
## 3 42.64  <NA>     C   <NA>
## 4 43.40  <NA>     O   <NA>
## 5 37.87  <NA>     C   <NA>
## 6 38.40  <NA>     C   <NA>
```


```r
pdb.ligand <- trim.pdb(pdb, inds=inds.ligand)
pdb.ligand
```

```
## 
##  Call:  trim.pdb(pdb = pdb, inds = inds.ligand)
## 
##    Total Models#: 1
##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
## 
##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 45  (residues: 1)
##      Non-protein/nucleic resid values: [ MK1 (1) ]
## 
## + attr: atom, helix, sheet, seqres, xyz,
##         calpha, call
```

```r
pdb.protein <-  trim.pdb(pdb,inds=inds.protein)
```



```r
write.pdb(pdb.ligand, file="1hsg_ligand.pdb")
write.pdb(pdb.protein, file="1hsg_protein.pdb")
```

















