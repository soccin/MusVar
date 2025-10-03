# Mouse Variant Pipeline

## Branch proj/17454

Version to work with Mouse WGS samples. This version/branch has the
updated varDict parallel calling and a custom BED file for project
17454 the deletes a region that is hyper-amplified

```
cat $ASSETS/Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_CLEAN__1000pad.bed \
  | grep -v "^14[[:space:]]19414856" \
  >proj_17454_fixed.bed
```

And a custom post script to use it.

- Fix config to have LSF params be set different for different assays

## Version: 1.3.2

Pulled in updates for sarek 3.4.4-A (for real this time). Also fixed memory in MarkDups

## Install

- Remember to clone with submodules `--recurse-submodules` 
- Install nextflow in `bin`

