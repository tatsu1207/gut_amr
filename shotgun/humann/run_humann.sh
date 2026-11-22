#!/bin/bash
# Community-wide functional profiling of the 72 KJM metagenomes (MetaCyc pathways).
# MetaPhlAn 4.1 (vJun23) -> HUMAnN 3.9 (ChocoPhlAn v201901_v31 + UniRef90 EC-filtered), host-removed reads.
# Same settings as ~/analyses/MAG_1/tools/{3_taxonomy,4_humann}.sh. Finished samples are skipped on rerun.
# Work (temp) on local disk; final tables copied to the vault.
set -u
RES=/home/unnot/vault/kjm/hiseq/results                  # <sample>/<sample>_nhost_{1,2}.fastq.gz
WORK=/home/unnot/analyses/KJM/humann_work                # local scratch
OUT=/home/unnot/vault/kjm/hiseq/humann/per_sample
MPA_DB=/home/database/metaphlan_vJun23
MPA_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202307
CHOCO=/home/database/humann3/chocophlan
UNIREF=/home/database/humann3/uniref_ec_filtered
JOBS=4; CPU=12
mkdir -p $WORK $OUT

for r1 in $RES/*/*_nhost_1.fastq.gz; do
    s=$(basename $r1 _nhost_1.fastq.gz)
    [ -f $OUT/${s}_pathabundance.tsv ] && { echo "skip $s"; continue; }
    while [ $(jobs -rp | wc -l) -ge $JOBS ]; do wait -n; done
    (
        d=$WORK/$s; rm -rf $d; mkdir -p $d
        cat $r1 ${r1/_1.fastq.gz/_2.fastq.gz} > $d/$s.fq.gz
        if ! conda run -n humann3_mpa41 metaphlan $d/$s.fq.gz --input_type fastq --nproc $CPU --offline \
                --bowtie2db $MPA_DB --index $MPA_INDEX --bowtie2out $d/mpa.bt2.bz2 \
                -o $d/${s}_metaphlan.tsv > $d/metaphlan.log 2>&1; then
            echo "FAIL metaphlan $s"; exit 1; fi
        if ! conda run -n humann3_mpa41 humann --input $d/$s.fq.gz --output $d --threads $CPU \
                --nucleotide-database $CHOCO --protein-database $UNIREF \
                --taxonomic-profile $d/${s}_metaphlan.tsv --remove-temp-output \
                --o-log $d/humann_detail.log > $d/humann.log 2>&1; then
            echo "FAIL humann $s"; exit 1; fi
        cp $d/${s}_metaphlan.tsv $d/${s}_genefamilies.tsv $d/${s}_pathcoverage.tsv $OUT/
        grep -E "reads|Unaligned|Total|Initial" $d/humann_detail.log | head -30 > $OUT/${s}_humann_stats.txt
        cp $d/${s}_pathabundance.tsv $OUT/          # copied last: marks the sample as finished
        rm -rf $d
        echo "done $s $(date +%T)"
    ) &
done
wait
echo "all finished"
