#!/usr/bin/env bash

config=/media/stephen/Data/20120810_LewisS_AM_MMhmeDIP/hmeDip-Seq/config/local-config.yaml
snakemake --use-conda -kp --notemp --configfile "${config}" \
--conda-prefix /media/stephen/Data/HiC-subsample/analysis/.snakemake "${@}"
