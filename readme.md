# GEOmirror

A package for researchers in mainland China to downlaod the GEO dataset,  which is just a **replacement** of getGEO function from GEOquery package.

### just one function: geoChina

Install the **development version** from Github:

```r
library(devtools)
install_github("jmzeng1314/GEOmirror")
library(GEOmirror)
```

Then use it to download GEO dataset, as below :

```
geoChina('gse1009') 
geoChina('GSE27533') 
geoChina('GSE95166') 
```

### warnings

- 1.only the expression profiling by array datasets will be offered by our package.
- 2.it's free for all of us, so we cann't make sure the internet connection will always be better than GEOquery. 

