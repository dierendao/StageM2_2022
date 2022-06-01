#!/bin/sh

# create PLZF_RARA.fa
perl -ne 'printf "%s\n",substr($_,0,60)' PLZF_
RARA_Origin.txt

# Find intersection of 2 files
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' "file1" "file2" > out4
comm -12 <(sort "file1") <(sort "file2") > out3
grep -Fxf "file1" "file2" > out2
join <(sort "file1") <(sort "file2") > out1 

# make custom GTF for PLZF_RARA (gene.gtf)
echo -e 'PLZF_RARA\tunknown\texon\t1\t3993\t.\t+\t.\tgene_id "PLZF_RARA"; transcript_id "PLZF_RARA"; gene_name "PLZF_RARA"; gene_biotype "protein_coding";' > gene_fusion.gtf

# copy ref.fa
cp Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus.GRCm39.dna.primary_assembly_PLZF.fa
 
# add gene.fa to ref.fa
cat gene_fusion.fa >> Mus_musculus.GRCm39.dna.primary_assembly_PLZF.fa

# copy ref.gtf
cp Mus_musculus.GRCm39.105.filtered.gtf Mus_musculus.GRCm39.105.filtered_PLZF.gtf

# add gene.gtf to ref.gtf
cat gene_fusion.gtf >> Mus_musculus.GRCm39.105.filtered_PLZF.gtf
