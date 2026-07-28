# Bash Scripts for Bioinformatics Pipelines

Collection of shell scripts and sbatch job files for HPC-based bioinformatics analyses.

## Subdirectories

- **Blast** -- BLAST sequence similarity search scripts (blastp, blastn, makeblastdb), plus HMMER protein domain search and phylogenetic tree construction with FastTree.
- **DAP_SEQ** -- DNA Affinity Purification sequencing pipeline for transcription factor binding site identification, including peak calling, motif discovery, and HOMER annotation.
- **WGS** -- Whole Genome Sequencing variant calling and analysis pipelines using BWA, GATK, bcftools, and BSA (Bulked Segregant Analysis).
- **GWAS** -- Genome-Wide Association Study scripts for running GEMMA, building VCF-based genotype matrices, setting covariates, and generating Manhattan/Q-Q plots.
- **RNA_SEQ** -- RNA-seq analysis pipelines for read processing, Kallisto-based quantification, differential expression analysis, and GO enrichment.
- **evolution** -- Evolutionary analysis scripts for CDS cleaning, codon-aware alignment (MAFFT + PAL2NAL), trimming, and maximum-likelihood tree construction with IQ-TREE.
- **ParaAT2.0** -- Paralog and ortholog identification using ParaAT2.0, with batch processing and Ka/Ks ratio calculation.
