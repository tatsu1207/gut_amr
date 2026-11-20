#!/bin/bash
# Search predicted proteins of all 72 HiSeq metagenomes for butyrate-pathway genes.
# hmmsearch: 36 KOfam profiles (butyrate pathways + single-copy markers)
# diamond:   butyryl-CoA:acetate CoA-transferase (but) positive/negative reference set
source ~/miniforge3/etc/profile.d/conda.sh; conda activate base
B=/home/unnot/vault/kjm/hiseq/butyrate_analysis; R=/home/unnot/vault/kjm/hiseq/results
run_one() {
  s=$1; B=/home/unnot/vault/kjm/hiseq/butyrate_analysis; R=/home/unnot/vault/kjm/hiseq/results
  [ -s $B/hmm_out/$s.tbl ] || hmmsearch --cpu 4 --noali -o /dev/null --tblout $B/hmm_out/$s.tbl.tmp --domtblout $B/hmm_out/$s.domtbl.tmp \
      $B/ref/butyrate_kos.hmm $R/$s/${s}_proteins.faa && mv $B/hmm_out/$s.tbl.tmp $B/hmm_out/$s.tbl && mv $B/hmm_out/$s.domtbl.tmp $B/hmm_out/$s.domtbl
  [ -s $B/diamond_out/$s.tsv ] || diamond blastp -q $R/$s/${s}_proteins.faa -d $B/ref/but_ref.dmnd -o $B/diamond_out/$s.tsv.tmp \
      --threads 4 --max-target-seqs 1 --evalue 1e-10 --quiet --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
      && mv $B/diamond_out/$s.tsv.tmp $B/diamond_out/$s.tsv
  echo "done $s"
}
export -f run_one
ls $R | grep '^ARG_' | xargs -P 36 -I{} bash -c 'run_one {}'
echo ALL_DONE
