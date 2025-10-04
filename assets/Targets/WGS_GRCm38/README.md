# WGS Resources

## Intervals

Appears we need to mask out problematic regions for WGS samples. So use the Boyle blacklist to create a masked version of the full genome to use. Also get rid of chrM just in case.


### Boyle-Lab/Blacklist

[https://github.com/Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist)

[mm10-blacklist.v2](https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/mm10-blacklist.v2.bed.gz)


### Generated Files

- **GRCm38.subBlack.bed**: Main genome intervals with blacklist
  regions masked out and chrM excluded. Used for both
  INTERVAL_BED_FILE and VARDICT_BED_FILE in target.resources.sh.

The blacklist masking process subtracts the problematic regions
from the original GRCm38.main.bed file to create the filtered
interval file used in variant calling.
