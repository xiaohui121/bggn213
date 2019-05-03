class06 HW
================
Xiaohui Lyu
2019/4/19

## Analyzes protein drug interactions

### **Step.1** read in data

``` r
library(bio3d)
```

### **Step.2** set up function

The function called *pdint* has 3 arguments: the protein name as **x**,
**the chain name** , which is *A* by default, the **elety**(a character
vector of atom names), which is *CA* by default. It mainly input the
name and the region of the protein and output the plot which describe
the protein drug interaction

``` r
pdint <- function(x,chain="A",elety="CA"){
#reading data from database
  s <- read.pdb(x)
#select the chain of the protein and create a new pdb object
  s.chain <- trim.pdb(s, chain=chain, elety=elety)
#select the region to be plotted
  s.b <- s.chain$atom$b
  plotb3(s.b, sse=s.chain, typ="l", ylab = "Bfactor")
} 
```

### **Step.3** try with 4AKE

``` r
pdint("4AKE")
```

    ##   Note: Accessing on-line PDB file

![](class6-bio3d_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
