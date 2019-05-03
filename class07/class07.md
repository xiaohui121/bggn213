Class 7: R functions and packages
================
Xiaohui Lyu
2019/4/24

## More on function writting

First we weill revisit our function from last day

``` r
source("http://tinyurl.com/rescale-R")
```

Test the **rescale()**
    function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#rescale(c(1:10,"string"))
```

``` r
x <- c(1:10,"string")
!is.numeric(x)
```

    ## [1] TRUE

``` r
##rescale2(x)
```

## Function practice

Write a function to identify NA elements in two vectors

Start with a simple example input where I know what the answer should
be.

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

I am looking for the position where it is TRUE in both vectors.

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum to find how many

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

This is my working snippet of code\!\!

``` r
both_na <- function(x,y){
  sum(is.na(x)& is.na(y))
}
```

``` r
both_na(x,y)
```

    ## [1] 1

``` r
both_na(c(NA,NA,NA),c(NA,NA,1))
```

    ## [1] 2

``` r
both_na(c(NA,NA,NA),c(NA,NA,1,NA))
```

    ## Warning in is.na(x) & is.na(y): 长的对象长度不是短的对象长度的整倍数

    ## [1] 3

Check the length of inputs are equal

``` r
x <- c(NA,NA,NA)
y <- c(NA,NA,1,NA,NA,NA)
length(x)
```

    ## [1] 3

``` r
length(y)
```

    ## [1] 6

``` r
length(x) == length(y)
```

    ## [1] FALSE

``` r
length(x) != length(y)
```

    ## [1] TRUE

Add the checking into the function

``` r
both_na2 <- function(x, y) {
  if(length(x) != length(y)) {
    stop("Input x and y should be the same length")
  }
  sum( is.na(x) & is.na(y) )
}
```

Try both\_na3 with extra features

``` r
both_na3 <- function(x, y) {
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length")
  }
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number <- sum(na.in.both)
  na.which <- which(na.in.both)
  message("Found ", na.number, " NA's at position(s):",
  paste(na.which, collapse=", ") )
  return( list(number=na.number, which=na.which) )
}
```

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

## Example grade function

``` r
grade <- function(x){
  (sum(x,na.rm = TRUE)-min(x,na.rm=TRUE))/(length(x)-1)
}
```

``` r
grade(c(100,100,100,100,100,100,100,90))
```

    ## [1] 100

``` r
grade(c(100,NA,90,90,90,90,97,80))
```

    ## [1] 79.57143

``` r
students <- read.csv("student_homework.csv",row.names = 1)
ans <- apply(students,1,grade)
sort(ans,decreasing = TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

## One last function

Find the intersection of two sets

``` r
source("http://tinyurl.com/rescale-R")
x <- df1$IDs
y <- df2$IDs

intersect(x,y)
```

    ## [1] "gene2" "gene3"

learn intersect from google and then find *%in%* in see also which is
better

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
y %in% x
```

    ## [1]  TRUE FALSE  TRUE FALSE

select the sentences and choose from **code -\> extract function**

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y],
        y[y %in% x])
}
```

``` r
merge(df1,df2,by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1
