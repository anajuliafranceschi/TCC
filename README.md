# TCC 🍇 💻 🧬

# Removing control sequences from the main archive
CV07 - CV08 - CV09 - CV10 - CV11 - CV12 - CV17 - CV18 - CV19 - CV20 - CV21 - CV22 - CV29 - CV30 - CV31 - CV32 - CV46 - CV33 - CV40 - CV41 - CV42 - CV43 - CV44 - CV45

# Mapping the sequences with Vitis labrusca (Var. Concord) using HISAT2
hisat2-build /media/ext5tb/anajulia/montagem2/fungi_reads/ncbi_dataset/data/GCA_011039315.1/GCA_011039315.1_ASM1103931v1_genomic.fna reference_index_vlabrusca

#!/bin/bash

ref=/home/sporisor/marcella/r570/reference/SofficinarumxspontaneumR570_771_v2.0.fa
index=/home/sporisor/marcella/r570/reference/hisat2_index/SofficinarumxspontaneumR570_771_v2.0

# Index reference genome
#hisat2-build "$ref" "$index"

# Map reads against genome
for i in /home/sporisor/marcella/r570/SP_inoc/ssci_free/*PE1.fastq.gz
do
        file2=`echo "$i" | sed "s/PE1/PE2/g"`
        rep=$(basename "$i" | sed -E "s/SP80-3280_inoc_rep([0-9]+)_*/\1/g")

        echo "Analyzing $i and $file2"

source activate bioinfo

        hisat2 -p 100 --rg-id SP80-3280_inoc --rg SM:SP80-3280_inoc \
        --summary-file ./summary_SP80-3280_inoc_rep_"$rep"_mapped_r570.txt \
        -x "$index" -1 "$i" -2 "$file2" \
        -S ./SP80-3280_inoc_rep_"$rep"_mapped_r570.sam

conda deactivate

module load Bio/SAMtools/1.20

        # Map convert to bam and index mapping
        samtools sort ./SP80-3280_inoc_rep_"$rep"_mapped_r570.sam > ./SP80-3280_inoc_rep_"$rep"_mapped_r570.bam
        samtools index ./SP80-3280_inoc_rep_"$rep"_mapped_r570.bam
        rm ./SP80-3280_inoc_rep_"$rep"_mapped_r570.sam

done










