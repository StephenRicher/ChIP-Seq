#!/usr/bin/env bash

config=/home/u/sr467/scratch/projects/20120810_LewisS_AM_MMhmeDIP/config/config.yaml
snakemake --profile slurm --use-conda -kp --notemp --configfile "${config}" \
--conda-prefix /home/u/sr467/scratch/projects/ChIP-Seq/.snakemake \
--cluster-config /home/u/sr467/scratch/projects/20120810_LewisS_AM_MMhmeDIP/config/cluster-config.yaml \
--wrapper-prefix file:/home/u/sr467/scratch/snakemake-wrappers/ "${@}"
