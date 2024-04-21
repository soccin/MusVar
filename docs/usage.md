# MusVar: Usage

First must create a manifest file `Proj_No_sample_manifest.csv` with the following columns:

|sample|patient |type|
|------|--------|----|
|S-ID-1|P-ID-1  |Tumor/Normal|
|507   |C-WK3YU0|Normal|
|508   |C-WK3YU0|Tumor|

Also make sure to add the pooled normal to the mapping file and manifest if needed.


Then
```
Rscript MusVar/bin/makeSarekInputSomatic.R sample_mapping.txt
```
to make the input .csv and then

```
./MusVar/runSarek.sh sarek_input_somatic.csv &
```

