# SPH5105 Assessment 2: kriging case counts across Singapore

R code for SPH5105 Spatial Epidemiology, Individual Assessment 2 (Hein Thu).

`assignment2_ht.r` is the complete analysis:

1. Data preparation
2. Variogram fitting
3. Ordinary kriging on a 1 km grid
4. Green space
5. Outbreak profiling and intervention centres
6. Robustness

## Running it

Open R (4.1 or later; tested on 4.6.0) in an empty folder and run:

```r
source("assignment2_ht.r")
```

- **Data:** the 30-site case table from the brief is embedded in the script. On the first run it downloads:
  - the Singapore boundary (GADM 4.1);
  - MRT/LRT station exits and bus stops (Land Transport Authority, via data.gov.sg).
- **Packages:** any missing ones are installed.
- **Outputs:** all figures and tables go to `output/`. Every number quoted in the report is in `output/table8_key_statistics.csv`.
- **Run time:** about 20 minutes the first time, mostly the nested cross-validation. Later runs take about 2 minutes because that result is cached in `data/derived/`.

## Data licences

- GADM: free for academic use, redistribution not permitted, so it is downloaded rather than included.
- LTA datasets: Singapore Open Data Licence.
