# gut_amr

Analysis code for a mouse study of how a high-fat diet (HFD) shapes the gut microbiome, resistome and
butyrate-synthesis capacity before, during and after treatment with amoxicillin, ciprofloxacin or
erythromycin. Mice were fed a normal diet or HFD for 13 weeks, given one antibiotic for 3 days and
followed for 14 days by 16S rRNA gene amplicon sequencing (8 timepoints), shotgun metagenomics
(3 timepoints, 72 samples) and gas chromatography of cecal short-chain fatty acids.

Raw reads: NCBI SRA, BioProject **PRJNA1501316**.

## Layout

| Folder | Scripts | Produces |
|---|---|---|
| `16S/` | `create_metadata.py` | sample metadata from ASV-table sample names |
| | `run_diversity.R` | alpha/beta diversity, PERMANOVA (Supplementary Table 1) |
| | `run_sseries_stats.R`, `run_pairwise_permdisp.R` | repeated-measures models, PERMDISP (Supplementary Data 2) |
| | `differential_abundance.R`, `diet_comparisons.R` | ALDEx2 differential abundance |
| | `genus_stats.py` | every genus percentage quoted in the Results (451 samples, >= 5,000 reads) |
| | `d0_baseline.R` | pre-diet (diet start, samples `ARG_D0_*`) vs T0 diet comparison (Supplementary Data 2) |
| | `taxonomy_plots.py` | Fig. 3 panels (genus % of genus-assigned reads) |
| | `plot_fig1_physiology.py`, `plot_fig2_diversity.py`, `plot_supplementary_figures.py` | Figs 1, 2 and Supplementary Figs 1, 2 |
| `shotgun/preprocessing_resistome/` | `run_resistome.sh` | fastp QC, host removal (Bowtie2, GRCm39), Bowtie2 mapping to CARD, `samtools idxstats` |
| | `run_arg_quantify.sh`, `calculate_arg_abundance.py` | 16S (SILVA) read counts and ARG copies per 16S rRNA gene copy |
| `shotgun/resistome/` | `reanalyze_arg_model_types.py` | per-gene abundance under three CARD model filters (all / rRNA variants excluded / protein homolog models only) |
| | `readlevel_stats.py` | **every read-level statistic quoted in the manuscript**, incl. the three-way ANOVA (Supplementary Data 4); main analysis = `homolog_only` |
| | `resistome_figures.py` | Fig. 4A-E |
| `shotgun/assembly_binning/` | `run_assembly*.sh` | MEGAHIT, Prodigal, RGI (CARD, Strict + Perfect), geNomad, MOB-typer |
| | `run_binning*.sh`, `run_gtdbtk.sh`, `run_mag_analysis.sh`, `run_analysis.sh` | MetaBAT2, MaxBin2, DAS Tool, CheckM2, GTDB-Tk, ARG-MAG linkage |
| | `run_arg_reanalysis.sh` | AMRFinderPlus cross-check (also an ARGs-OAP run not used in the paper) |
| | `plot_mag_figures.py` | Fig. 5 panels |
| `shotgun/mobilome/` | `run_mobilome_analysis.sh`, `analyze_mobilome.py`, `*_plasmid_summaries.sh` | ARG mobility (plasmid/phage/conjugation) and Supplementary Data 11 |
| `shotgun/butyrate/` | `run_search.sh`, `analyze_butyrate.py`, `plot_fig6_butyrate.py` | butyrate-pathway genes (KOfam/HMMER, DIAMOND for But), Fig. 6, Supplementary Data 12 |
| `shotgun/lmrD/` | `classify_mobility.py` | mobility rule applied to the lmrD homologues (Supplementary Data 8) |
| `manuscript/` | `compose_multipanel_figures.py`, `build_*`, `generate_manuscript_pdf.py`, `export_panels.py` | figure assembly, one PNG per panel, and the submission documents |

## Notes

- Timepoints: data files and sample IDs use the original codes (T0h, T24h, T72h, S72h, SD7, SD10, SD14,
  matching SRA); figures and the manuscript show T0, T1, T3, R3, R7, R10, R14 (T = days since the first
  dose, R = days since the last dose). `manuscript/tplabels.py` converts the labels when a figure is saved.
- Scripts were run in place on the authors' server and contain absolute paths (`/home/unnot/vault/kjm/...`);
  edit these before running elsewhere.
- Read-level resistome results use **CARD protein homolog models only**. rRNA gene variant, protein
  variant, overexpression and knockout models describe wild-type genes whose resistance depends on
  mutations that read mapping cannot genotype; `readlevel_stats.py no_rrna` reproduces the sensitivity
  analysis that retains the protein-coding variant models.
- Never sum CARD "simple class" columns to obtain totals: a gene annotated to k drug classes appears in
  all k columns. Totals here are sums over genes.
- Not included: the DADA2 amplicon-processing script (parameters are given in the Methods), and the
  interactive steps behind Supplementary Data 9 (DIAMOND self-scores of CARD homolog models and the direct
  homology search of the 72 assemblies) and the lmrD direct search; their outputs are in the
  supplementary tables. `classify_mobility.py` re-implements the mobility rule described in the Methods.

## Main software

R 4.3 (DADA2, vegan, ALDEx2), SILVA 138.1, fastp, Bowtie2, samtools, CARD/RGI (DIAMOND aligner),
MEGAHIT, Prodigal, MetaBAT2, MaxBin2, DAS Tool, CheckM2, skani, GTDB-Tk v2.4 (GTDB release 226),
geNomad, MOB-suite, HMMER v3.4 with KOfam (KEGG release 2024-01), DIAMOND v2.1.8,
AMRFinderPlus (database 2026-03-24.1); Python 3 with pandas, SciPy and statsmodels.

## License

MIT (see `LICENSE`).
