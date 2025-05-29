read_as_metrics<-function(ff) {
    read_tsv(ff,comment="#",show_col_types=F,progress=F) %>%
        filter(CATEGORY=="PAIR") %>%
        select(SAMPLE,everything()) %>%
        select(-LIBRARY,-READ_GROUP,-PCT_HARDCLIP)
}

require(tidyverse)
asm=fs::dir_ls("post/metrics",recur=T,regex=".as.txt$") %>%
    map(read_as_metrics) %>%
    bind_rows
