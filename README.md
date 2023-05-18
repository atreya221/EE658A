# EE658A FUZZY PROJECT


## Generating Fuzzy Rules by Learning from Examples - Wang & Mendel


### Truck Backer-Upper Control Problem


For truck backer-upper control, run 


1. `truckBackerUpper_all.m`  for example 1 - with complete numerical data (no linguistic rules)
2. `truckBackerUpper_truncated.m`  for example 2 - with truncated data pairs (first 3 pairs from each sequence)
3. `truckBackerUpper_truncated_linguistic.m`  for example 2 - with truncated data pairs and added linguistic rules


### Mackey-Glass Chaotic Time Series Prediction Problem


For Mackey-Glass time series prediction, run 


1. `mackeyGlass_700.m`  for example 1 - training using 700 data points of `x(t)`
2. `mackeyGlass_200.m`  for example 2 - training using 200 data points of `x(t)`


For generating Mackey-Glass chaotic time series, I used the files `mackeyglass.m`, `mackeyglass_eq.m` and `mackeyglass_rk4.m` from the following source [https://www.mathworks.com/matlabcentral/fileexchange/24390-mackey-glass-time-series-generator](https://www.mathworks.com/matlabcentral/fileexchange/24390-mackey-glass-time-series-generator)


### Project directory structure

```bash
│   EE658A_Project_Presentation_pdf.pdf
│   EE658A_Project_Presentation_ppt.pptx
│   EE658A_Project_Report.pdf
│   README.md
│
├───MackeyGlassChaoticTimeSeries
│       mackeyglass.m
│       mackeyGlass_200.m
│       mackeyGlass_700.m
│       mackeyglass_eq.m
│       mackeyglass_rk4.m
│       mgchaotic.dat
│       timeseriesregions.m
│
└───TruckBackerUpper
    │   truckBackerUpper_all.m
    │   truckBackerUpper_truncated.m
    │   truckBackerUpper_truncated_linguistic.m
    │
    └───data
            table1.csv
            table10.csv
            table11.csv
            table12.csv
            table13.csv
            table14.csv
            table2.csv
            table3.csv
            table4.csv
            table5.csv
            table6.csv
            table7.csv
            table8.csv
            table9.csv
```
