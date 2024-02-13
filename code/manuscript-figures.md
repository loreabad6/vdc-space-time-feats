

# Install packages

Run this chunk to get the necessary packages. Note that, to reproduce,
some packages require the development version fetched from GitHub.

``` r
pkgs = c(
  "dplyr","patchwork",
  "purrr","readr", "stars",
  "stringr","units","zen4R"
)
install.packages(pkgs)
remotes::install_github("tidyverse/ggplot2")
remotes::install_github("huizezhang-sherry/cubble")
```

# Libraries

``` r
library(cubble)
library(dplyr)
library(ggplot2)
library(here)
library(patchwork)
library(purrr)
library(readr)
library(stars)
library(stringr)
library(units)
library(zen4R)
```

# Figure 1

Created with sketchpad.io and draw.io See source file in
code/vdc-schema.drawio

# Figure 2

Code snippet with the result of the following script The script creates
a VDC for Fagradalsfjall It is in array format, using the `stars`
package

Data comes from Zenodo so the script shows demo code with steps to 1.
download to temporary directory 2. pre-process 3. coerce to a vdc

## Download

``` r
dir = tempdir()
download_zenodo(
  doi = "10.5281/zenodo.7866738",
  path = dir,
  files = list(
    "outlines_pedersen_etal2022_v12.zip"
  ),
  overwrite = FALSE,
  timeout = 600
)
```

    [zen4R][INFO] ZenodoRecord - Download in sequential mode 
    [zen4R][INFO] ZenodoRecord - Will download 1 file from record '7866738' (doi: '10.5281/zenodo.7866738') - total size: 2.3 MiB 
    [zen4R][INFO] Downloading file 'outlines_pedersen_etal2022_v12.zip' - size: 2.3 MiB
    [zen4R][INFO] File downloaded at 'C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs'.
    [zen4R][INFO] ZenodoRecord - Verifying file integrity... 
    [zen4R][INFO] File 'outlines_pedersen_etal2022_v12.zip': integrity verified (md5sum: fc7a74d235274b9707bc5acb296b45ea)
    [zen4R][INFO] ZenodoRecord - End of download 

``` r
# Unzip
files = list.files(here(dir), full.names = TRUE)
lapply(files, unzip, exdir = here(dir, "unzipped"))
```

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    [[1]]
    NULL

    [[2]]
    NULL

    [[3]]
    NULL

    [[4]]
     [1] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210320_1240_A6D_Pedersen_etal2022_v12.gpkg"
     [2] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210321_1130_HEL_Pedersen_etal2022_v12.gpkg"
     [3] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210331_1210_A6D_Pedersen_etal2022_v12.gpkg"
     [4] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210405_1010_A6D_Pedersen_etal2022_v12.gpkg"
     [5] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210405_1416_A6D_Pedersen_etal2022_v12.gpkg"
     [6] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210406_1338_A6D_Pedersen_etal2022_v12.gpkg"
     [7] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210408_1325_A6D_Pedersen_etal2022_v12.gpkg"
     [8] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210412_1210_A6D_Pedersen_etal2022_v12.gpkg"
     [9] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210320_0745_HEL_Pedersen_etal2022_v12.gpkg"
    [10] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210322_1322_PLE_Pedersen_etal2022_v12.gpkg"
    [11] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210323_1005_A6D_Pedersen_etal2022_v12.gpkg"
    [12] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210326_1252_PLE_Pedersen_etal2022_v12.gpkg"
    [13] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210329_1319_PLE_Pedersen_etal2022_v12.gpkg"
    [14] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210330_1311_PLE_Pedersen_etal2022_v12.gpkg"
    [15] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210418_1230_A6D_Pedersen_etal2022_v12.gpkg"
    [16] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210421_1330_A6D_Pedersen_etal2022_v12.gpkg"
    [17] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210426_1515_A6D_Pedersen_etal2022_v12.gpkg"
    [18] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210428_1249_PLE_Pedersen_etal2022_v12.gpkg"
    [19] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210503_1545_A6D_Pedersen_etal2022_v12.gpkg"
    [20] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210510_1242_A6D_Pedersen_etal2022_v12.gpkg"
    [21] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210518_1730_A6D_Pedersen_etal2022_v12.gpkg"
    [22] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210602_1522_A6D_Pedersen_etal2022_v12.gpkg"
    [23] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210611_1250_A6D_Pedersen_etal2022_v12.gpkg"
    [24] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210626_1330_A6D_Pedersen_etal2022_v12.gpkg"
    [25] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210702_1249_PLE_Pedersen_etal2022_v12.gpkg"
    [26] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210727_1000_A6D_Pedersen_etal2022_v12.gpkg"
    [27] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210808_1717_A6D_Pedersen_etal2022_v12.gpkg"
    [28] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210909_1600_A6D_Pedersen_etal2022_v12.gpkg"
    [29] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210917_1330_A6D_Pedersen_etal2022_v12.gpkg"
    [30] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210930_1420_A6D_Pedersen_etal2022_v12.gpkg"

``` r
# Find geopackage files
fn_gpkg = list.files(
  here(dir, "unzipped"),
  pattern = "Outline.*gpkg$",
  full.names = TRUE,
  recursive = TRUE
)
```

## Pre-process

``` r
# Create function to read in each file and to extract 
# the date from the filename
read_fun = function(x){
  read_sf(x) |> 
    mutate(
      fn_gpkg = tail(str_split(x, '/')[[1]], n=1),
      datetime = as.POSIXct(
        paste(str_split(fn_gpkg, "_")[[1]][2:3], collapse = ""),
        format = "%Y%m%d%H%M"
      )
    ) |> 
    st_set_crs(3057)
}

# Call map2 to read in the files
# and bind them into one single sf object
outlines = map(fn_gpkg, read_fun) |> bind_rows()
```

    Warning: st_crs<- : replacing crs does not reproject data; use st_transform for
    that

``` r
# Combine polygons from single date into multipolygons
# Make geometry valid
outlines = outlines |> 
  group_by(datetime) |> 
  summarise(geom = st_combine(geom)) |> 
  ungroup() |> 
  st_make_valid()
```

## VDC coercion

``` r
# Add gid identifier
outlines = outlines |> 
  mutate(gid = 1L)

# Create array
a = array(
  data = outlines$geom, 
  dim = c(
    length(unique(outlines$gid)),
    length(unique(outlines$gid)),
    length(unique(outlines$datetime))
  ),
  dimnames = list(
    geom_sum = unique(outlines$gid),
    gid = unique(outlines$gid),
    datetime = unique(outlines$datetime)
  )
)

# Create dimensions object
# Summary geometry is the centroid of the union of all geometries
# The point parameter indicates if the value refers to a point (location)
# or to a pixel (area) value
dim_cent = st_dimensions(
  geom_sum = st_centroid(st_make_valid(st_union(outlines$geom))), # approach centroid
  gid = unique(outlines$gid),
  datetime = unique(outlines$datetime),
  point = c(TRUE, TRUE, FALSE)
)

# Coerce to cube
# The output of this code constitutes Figure 2
(cube_arr = st_as_stars(
  list(geometry = a), 
  dimensions = dim_cent)
)
```

    stars object with 3 dimensions and 1 attribute
    attribute(s):
             geometry  
     MULTIPOLYGON : 2  
     POLYGON      :28  
     epsg:3057    : 0  
     +proj=lcc ...: 0  
    dimension(s):
             from to               refsys point
    geom_sum    1  1 ISN93 / Lambert 1993  TRUE
    gid         1  1                   NA  TRUE
    datetime    1 30              POSIXct FALSE
                                                  values
    geom_sum                       POINT (339860 380008)
    gid                                                1
    datetime 2021-03-20 07:45:00,...,2021-09-30 14:20:00

# Figure 3

Base R plot of lava flow outlines a) shows all the lava flows for the
whole time period b) shows the result of filtering the VDC between two
dates Plots are combine with the patchwork library, hence the use of
`wrap_elements()`

``` r
oldpar = par(no.readonly = TRUE)
par(mar = c(1.8,1,1.1,0.3), bg = "transparent")
wrap_elements(
  ~cube_arr |>
    # use pseudo-function to extract geometries from vdc attributes
    (\(x) plot(x$geometry, axes = TRUE))(),
  clip = FALSE
) /
wrap_elements(
  ~cube_arr |>
    filter(datetime > "2021-03-18", datetime < "2021-03-25") |> 
    # use pseudo-function to extract geometries from vdc attributes
    (\(x) plot(x$geometry, axes = TRUE))(),
  clip = FALSE
) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
```

![](manuscript-figures_files/figure-commonmark/unnamed-chunk-6-1.png)

``` r
par(oldpar)
```

# Figure 4

Code snippet with the result of the following script The script creates
a VDC for Butangbunasi It is in tabular format, using the `cubble`
package

Data comes from Zenodo so the script shows demo code with steps to 1.
download to temporary directory 2. pre-process 3. coerce to a vdc

## Download

``` r
# Fetch also the CSV file with addiontional info
dir = tempdir()
download_zenodo(
  doi = "10.5281/zenodo.10635102",
  path = dir,
  files = list(
    "outlines.zip",
    "Butangbunasi_OBIA_statistics.csv"
  ),
  overwrite = FALSE,
  timeout = 100
)
```

    [zen4R][INFO] ZenodoRecord - Download in sequential mode 
    [zen4R][INFO] ZenodoRecord - Will download 2 files from record '10635102' (doi: '10.5281/zenodo.10635102') - total size: 144 KiB 
    [zen4R][INFO] Downloading file 'outlines.zip' - size: 142.7 KiB
    [zen4R][INFO] Downloading file 'Butangbunasi_OBIA_statistics.csv' - size: 1.2 KiB
    [zen4R][INFO] Files downloaded at 'C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs'.
    [zen4R][INFO] ZenodoRecord - Verifying file integrity... 
    [zen4R][INFO] File 'outlines.zip': integrity verified (md5sum: d0034de915b5cae20d9be02899550e9a)
    [zen4R][INFO] File 'Butangbunasi_OBIA_statistics.csv': integrity verified (md5sum: 8ebe12e999df3e41f862db5dc31a57e7)
    [zen4R][INFO] ZenodoRecord - End of download 

``` r
# Unzip
files = list.files(here(dir), full.names = TRUE)
lapply(files, unzip, exdir = here(dir, "unzipped"))
```

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    Warning in FUN(X[[i]], ...): error 1 in extracting from zip file

    [[1]]
    NULL

    [[2]]
    NULL

    [[3]]
    NULL

    [[4]]
    NULL

    [[5]]
    NULL

    [[6]]
    NULL

    [[7]]
    NULL

    [[8]]
    NULL

    [[9]]
     [1] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1984.gpkg"      
     [2] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1989.gpkg"      
     [3] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1990.gpkg"      
     [4] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1992.gpkg"      
     [5] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1994.gpkg"      
     [6] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1996.gpkg"      
     [7] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_1998.gpkg"      
     [8] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2000.gpkg"      
     [9] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2001.gpkg"      
    [10] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2004.gpkg"      
    [11] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2005.gpkg"      
    [12] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2005_10_03.gpkg"
    [13] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2008.gpkg"      
    [14] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2008_03_17.gpkg"
    [15] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2009.gpkg"      
    [16] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2010.gpkg"      
    [17] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2013.gpkg"      
    [18] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2015.gpkg"      
    [19] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2016.gpkg"      
    [20] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2018.gpkg"      
    [21] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines/Butangbunasi_2021.gpkg"      

    [[10]]
     [1] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210320_1240_A6D_Pedersen_etal2022_v12.gpkg"
     [2] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210321_1130_HEL_Pedersen_etal2022_v12.gpkg"
     [3] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210331_1210_A6D_Pedersen_etal2022_v12.gpkg"
     [4] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210405_1010_A6D_Pedersen_etal2022_v12.gpkg"
     [5] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210405_1416_A6D_Pedersen_etal2022_v12.gpkg"
     [6] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210406_1338_A6D_Pedersen_etal2022_v12.gpkg"
     [7] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210408_1325_A6D_Pedersen_etal2022_v12.gpkg"
     [8] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210412_1210_A6D_Pedersen_etal2022_v12.gpkg"
     [9] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210320_0745_HEL_Pedersen_etal2022_v12.gpkg"
    [10] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210322_1322_PLE_Pedersen_etal2022_v12.gpkg"
    [11] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210323_1005_A6D_Pedersen_etal2022_v12.gpkg"
    [12] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210326_1252_PLE_Pedersen_etal2022_v12.gpkg"
    [13] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210329_1319_PLE_Pedersen_etal2022_v12.gpkg"
    [14] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210330_1311_PLE_Pedersen_etal2022_v12.gpkg"
    [15] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210418_1230_A6D_Pedersen_etal2022_v12.gpkg"
    [16] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210421_1330_A6D_Pedersen_etal2022_v12.gpkg"
    [17] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210426_1515_A6D_Pedersen_etal2022_v12.gpkg"
    [18] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210428_1249_PLE_Pedersen_etal2022_v12.gpkg"
    [19] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210503_1545_A6D_Pedersen_etal2022_v12.gpkg"
    [20] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210510_1242_A6D_Pedersen_etal2022_v12.gpkg"
    [21] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210518_1730_A6D_Pedersen_etal2022_v12.gpkg"
    [22] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210602_1522_A6D_Pedersen_etal2022_v12.gpkg"
    [23] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210611_1250_A6D_Pedersen_etal2022_v12.gpkg"
    [24] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210626_1330_A6D_Pedersen_etal2022_v12.gpkg"
    [25] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210702_1249_PLE_Pedersen_etal2022_v12.gpkg"
    [26] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210727_1000_A6D_Pedersen_etal2022_v12.gpkg"
    [27] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210808_1717_A6D_Pedersen_etal2022_v12.gpkg"
    [28] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210909_1600_A6D_Pedersen_etal2022_v12.gpkg"
    [29] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210917_1330_A6D_Pedersen_etal2022_v12.gpkg"
    [30] "C:/Users/b1066081/AppData/Local/Temp/Rtmp6RYJQs/unzipped/outlines_pedersen_etal2022_v12/Outline_20210930_1420_A6D_Pedersen_etal2022_v12.gpkg"

    [[11]]
    NULL

``` r
# Find geopackage files
mapping_ls = list.files(
  here(dir, "unzipped", "outlines"), 
  pattern = ".gpkg",
  full.names = TRUE
)

# Read in CSV file
stats = read_csv(
  here(dir, "Butangbunasi_OBIA_statistics.csv")
)
```

    Rows: 21 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr  (2): satellite_sensor, file
    dbl  (2): landslide_area_ha, lake_area_ha
    date (1): date

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Set-up a read function that fetches filename
read_fun = function(x){
  read_sf(x) |> 
    mutate(
      fn_gpkg = tail(str_split(x, '/')[[1]], n=1),
    )
}

# Read files and combine into single sf object
# Join with CSV file
mapping = lapply(mapping_ls, read_fun) |> 
  bind_rows() |> 
  left_join(stats, by = c("fn_gpkg" = "file"))
```

## Pre-process

``` r
# All pre-processing steps are chained but have 
# documentation in each line
landslides = mapping |> 
  # Extract only landslide class
  filter(Class == "landslide") |> 
  # Coerce date column to Date class
  mutate(date = as.Date(date)) |> 
  # Remove outline for 2018
  # this is a reference outline and does not follow a
  # Typhoon event according to Hoelbling et al., 2020
  filter(date != as.Date("2018-11-08")) |>
  # Group by date which in combination with summarise
  # will union all single polygons into a multipolygon
  group_by(date) |> 
  # Summarise sensor and area information for landslide
  summarise(
    sensor = first(satellite_sensor),
    area = set_units(first(landslide_area_ha), 'ha')
  ) |> 
  # Add unique identifier
  mutate(gid = 1L)

# Compute geom_sum as the centroid of the union of the geometries
landslidescent = landslides |> 
  group_by(gid) |> 
  mutate(geom_sum = st_centroid(st_union(geom))) |> 
  relocate(geom_sum, .after = gid) |> 
  st_as_sf(sf_column_name = 'geom_sum')
```

## VDC coercion

``` r
# Call as cubble to coerce from sf to cubble
# key is the gid column while index is date
cube_tab = as_cubble(
  landslidescent, key = gid, index = date
)

# The output of this code corresponds to Figure 4a
cube_tab |> 
  face_spatial()
```

    ℹ The cubble is already in the nested form

    # cubble:   key: gid [1], index: date, nested form, [sf]
    # spatial:  [271664.917737363, 2568861.92906261, 271664.917737363,
    #   2568861.92906261], WGS 84 / UTM zone 51N
    # temporal: date [date], sensor [chr], area [[ha]], geom [MULTIPOLYGON [m]]
        gid       x        y           geom_sum ts               
    * <int>   <dbl>    <dbl>        <POINT [m]> <list>           
    1     1 271665. 2568862. (271664.9 2568862) <tibble [20 × 4]>

``` r
# The truncated output of this code corresponds to Figure 4b
cube_tab |> 
  face_temporal()
```

    # cubble:   key: gid [1], index: date, long form
    # temporal: 1984-12-12 -- 2021-08-28 [8D], has gaps!
    # spatial:  x [dbl], y [dbl], geom_sum [POINT [m]]
         gid date       sensor     area                                         geom
       <int> <date>     <chr>      [ha]                           <MULTIPOLYGON [m]>
     1     1 1984-12-12 Landsat 5  66.2 (((271637.5 2568620, 271637.5 2568583, 2716…
     2     1 1989-10-23 Landsat 5  62.4 (((273712.5 2566845, 273562.5 2566845, 2735…
     3     1 1990-10-10 Landsat 5  78.2 (((273737.5 2566845, 273737.5 2566820, 2735…
     4     1 1992-10-31 Landsat 5 121.  (((273462.5 2566820, 273462.5 2566870, 2734…
     5     1 1994-09-03 Landsat 5  94.1 (((273500 2566820, 273462.5 2566820, 273437…
     6     1 1996-08-23 Landsat 5 118.  (((273462.5 2566820, 273462.5 2566845, 2734…
     7     1 1998-11-01 Landsat 5  96.3 (((273775 2566820, 273712.5 2566820, 273500…
     8     1 2000-09-27 Landsat 7 120.  (((273500 2566783, 273500 2566820, 273437.5…
     9     1 2001-09-14 Landsat 7 121.  (((273825 2566820, 273825 2566783, 273800 2…
    10     1 2004-07-12 Landsat 5 123.  (((273825 2566933, 273862.5 2566933, 273862…
    11     1 2005-09-17 Landsat 5 151.  (((273787.5 2566683, 273775 2566683, 273775…
    12     1 2005-10-03 Landsat 5 151.  (((273712.5 2566870, 273712.5 2566845, 2737…
    13     1 2008-03-17 Landsat 5 191.  (((273737.5 2566758, 273675 2566758, 273675…
    14     1 2008-08-24 Landsat 5 164.  (((273775 2566783, 273737.5 2566783, 273462…
    15     1 2009-09-12 Landsat 5 424.  (((272062.5 2568133, 272000 2568133, 272000…
    16     1 2010-12-20 Landsat 5 396.  (((273650 2566658, 273612.5 2566658, 273612…
    17     1 2013-06-03 Landsat 8 399.  (((273262.5 2566370, 273275 2566370, 273275…
    18     1 2015-11-16 Landsat 8 396.  (((272175 2568258, 272237.5 2568258, 272237…
    19     1 2016-12-04 Landsat 8 382.  (((272087.5 2568220, 272112.5 2568220, 2721…
    20     1 2021-08-28 Landsat 8 388.  (((272062.5 2568108, 272000 2568108, 272000…

# Figure 5

This figure shows in a) a multi-dimensional representation of the
Butangbunasi landslide and b) a time series of the landslide area To
create a) the geometry of the original landslide data (outside of the
cube) was distorted by multiplying it by a shear matrix and adding a
shift in the y axis to stack them on top of each other

``` r
# Create shear matrix
sm = matrix(c(2.5, 1.2, 0, 1), 2, 2)

# Apply shear matrix
ldsl_shear = landslides |>
  mutate(
    area = st_area(geom),
    geom = geom * sm,
    # sequence along date, i.e. 1 to 20
    shift_y = seq_along(date)
  ) |> 
  # Add lost crs
  st_set_crs(st_crs(landslides)) 

ldsl_shift = ldsl_shear |> 
  rowwise() |> 
  # add a shift to stack outlines on top of each other
  mutate(
    geom = geom + c(0, shift_y * 4000),
    # Marker to place date labels
    y_label = st_coordinates(st_centroid(geom))[,'Y']
  ) |> 
  ungroup() |> 
  st_as_sf()

# Plot a)
spaceplot = ggplot(ldsl_shift) +
  geom_sf(
    aes(fill = date), 
    color = "black",
    show.legend = FALSE
  ) +
  # Add date label, the value of x is added after visual inspection
  geom_text(aes(label = date, y = y_label), x = 3751500, size = 3.5) +
  scale_fill_viridis_c("Date", direction = 1, trans = "date", option = "D") +
  # Limits are expanded after visual inspection
  coord_sf(xlim = c(3745000, 3767000), ylim = c(2576500, 2646000), clip = "off") +
  theme_void()+
  theme(
    text = element_text(size = 18)
  )


# Plot b), directly from cubble object in temporal face
area = cube_tab |> 
  face_temporal() |> 
  rename(Area = area) |> 
  ggplot() +
  aes(
    x = date, y = Area,
    color = date,
    shape = sensor, 
    group = 1
  ) +
  geom_point(size = 2) + geom_line() +
  scale_y_units(unit = "km2") +
  scale_color_viridis_c("Date", direction = 1, trans = "date", option = "D") +
  scale_shape("Sensor") +
  guides(
    color = guide_colorbar(barwidth = 20, barheight = 0.8)
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    axis.title.x = element_blank(),
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.title.position = "top"
  )

# Set-up layout
layout = "
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAABBBBBB
AAAACCCCCC
AAAADDDDDD
AAAADDDDDD
AAAADDDDDD
AAAADDDDDD
"

# Combine plots with patchwork
spaceplot + area + plot_spacer() + guide_area() + 
  plot_layout(guides = "collect", design = layout) + 
  plot_annotation(tag_levels = "a", tag_suffix = ")")
```

![](manuscript-figures_files/figure-commonmark/unnamed-chunk-10-1.png)

# R Session Info

``` r
sessionInfo()
```

    R version 4.3.1 (2023-06-16 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)

    Matrix products: default


    locale:
    [1] LC_COLLATE=English_Austria.utf8  LC_CTYPE=English_Austria.utf8   
    [3] LC_MONETARY=English_Austria.utf8 LC_NUMERIC=C                    
    [5] LC_TIME=English_Austria.utf8    

    time zone: Europe/Vienna
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] zen4R_0.9          units_0.8-5.3      stringr_1.5.1      stars_0.6-5       
     [5] sf_1.0-15          abind_1.4-5        readr_2.1.4        purrr_1.0.2       
     [9] patchwork_1.2.0    here_1.0.1         ggplot2_3.4.4.9000 dplyr_1.1.3       
    [13] cubble_0.3.0      

    loaded via a namespace (and not attached):
     [1] gtable_0.3.4       anytime_0.3.9      cubelyr_1.0.2      xfun_0.40         
     [5] rdflib_0.2.8       tsibble_1.1.4      tzdb_0.4.0         vctrs_0.6.5       
     [9] tools_4.3.1        generics_0.1.3     curl_5.0.2         parallel_4.3.1    
    [13] tibble_3.2.1       proxy_0.4-27       fansi_1.0.6        pkgconfig_2.0.3   
    [17] KernSmooth_2.23-22 redland_1.0.17-17  assertthat_0.2.1   lifecycle_1.0.4   
    [21] farver_2.1.1       compiler_4.3.1     atom4R_0.3-3       textshaping_0.3.6 
    [25] munsell_0.5.0      keyring_1.3.1      ncdf4_1.21         htmltools_0.5.7   
    [29] class_7.3-22       yaml_2.3.8         crayon_1.5.2       pillar_1.9.0      
    [33] tidyr_1.3.0        ellipsis_0.3.2     classInt_0.4-10    tidyselect_1.2.0  
    [37] zip_2.3.0          digest_0.6.33      stringi_1.7.12     labeling_0.4.3    
    [41] rprojroot_2.0.3    fastmap_1.1.1      grid_4.3.1         colorspace_2.1-0  
    [45] cli_3.6.2          magrittr_2.0.3     XML_3.99-0.16      utf8_1.2.4        
    [49] e1071_1.7-13       withr_3.0.0        scales_1.3.0       bit64_4.0.5       
    [53] lubridate_1.9.2    timechange_0.2.0   roxygen2_7.2.3     rmarkdown_2.25    
    [57] httr_1.4.7         bit_4.0.5          ragg_1.2.5         hms_1.1.3         
    [61] evaluate_0.23      knitr_1.45         viridisLite_0.4.2  gridGraphics_0.5-1
    [65] rlang_1.1.3        Rcpp_1.0.11        glue_1.7.0         DBI_1.2.1         
    [69] xml2_1.3.5         vroom_1.6.3        rstudioapi_0.15.0  jsonlite_1.8.8    
    [73] R6_2.5.1           systemfonts_1.0.4 
