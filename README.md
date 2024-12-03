# TCC 🍇 💻 🧬

# Removing control sequences from the main archive
CV07 - CV08 - CV09 - CV10 - CV11 - CV12 - CV17 - CV18 - CV19 - CV20 - CV21 - CV22 - CV29 - CV30 - CV31 - CV32 - CV46 - CV33 - CV40 - CV41 - CV42 - CV43 - CV44 - CV45

# Mapping the sequences with Vitis labrusca (Var. Concord) using HISAT2
        #!/bin/bash
        
        #Caminhos para o genoma de referência e o índice HISAT2
        ref=/media/ext5tb/anajulia/montagem2/fungi_reads/ncbi_dataset/data/GCA_011039315.1/GCA_011039315.1_ASM1103931v1_genomic.fna
        index=/media/ext5tb/anajulia/montagem2/fungi_reads/hisat2_index
        
        #Construir o índice do genoma de referência (descomente a linha abaixo se ainda não tiver o índice)
        hisat2-build "$ref" "$index"
        
        #Caminho para as leituras (substitua conforme necessário)
        reads_dir=/media/ext5tb/anajulia/montagem2/fungi_reads
        
        #Alinhar leituras contra o genoma de referência
        for i in "$reads_dir"/*_cut_PE1.fastq.gz
        do
            # Substituir _R1 por _R2 para os arquivos paired-end
            file2=$(echo "$i" | sed "s/_PE1/_PE2/g")
        
            # Extrair nome do arquivo base
            sample_name=$(basename "$i" | sed "s/_cut_PE1.fastq.gz//")
        
            echo "Analisando amostra: $sample_name"
        
            # Executar o HISAT2
            hisat2 -p 25 --rg-id "$sample_name" --rg SM:"$sample_name" \
            --summary-file ./summary_"$sample_name"_mapped.txt \
            -x "$index" -1 "$i" -2 "$file2" \
            -S ./"$sample_name"_mapped.sam
        
            # Converter SAM para BAM, ordenar e indexar
            samtools sort ./"$sample_name"_mapped.sam -o ./"$sample_name"_mapped.bam
            samtools index ./"$sample_name"_mapped.bam
        
            # Remover o arquivo SAM para economizar espaço
            rm ./"$sample_name"_mapped.sam
        
        done

# Removing the non mapped sequences (the ones that are probably fungi sequences)
        #!/bin/bash
        
        #Caminhos
        bam_dir="/media/ext5tb/anajulia/montagem2/fungi_reads/hisat2"  # Diretório onde os arquivos BAM estão localizados
        output_dir="/media/ext5tb/anajulia/montagem2/fungi_reads/unmapped_fastq"  # Diretório para os arquivos FASTQ de saída
        
        #Criar diretório de saída se não existir
        mkdir -p "$output_dir"
        
        #Loop para processar cada par de arquivos BAM
        for bam in "$bam_dir"/*_mapped.bam
        do
            # Obter o nome base do arquivo
            sample_name=$(basename "$bam" "_mapped.bam")
        
            echo "Processando $sample_name"
        
            # Extrair sequências não mapeadas usando Samtools
            samtools view -b -f 12 -F 256 "$bam" > "$output_dir"/"$sample_name"_unmapped.bam
        
            # Converter BAM de não mapeados para FASTQ usando Bedtools
            bedtools bamtofastq -i "$output_dir"/"$sample_name"_unmapped.bam \
            -fq "$output_dir"/"$sample_name"_unmapped_R1.fastq \
            -fq2 "$output_dir"/"$sample_name"_unmapped_R2.fastq
        
            # Remover arquivo BAM de não mapeados para economizar espaço
            rm "$output_dir"/"$sample_name"_unmapped.bam
        
            echo "$sample_name processado."
        done
        
        echo "Processamento concluído."

# Converting fastq files in fasta files
        #!/bin/bash
        
        #Diretório onde os arquivos FASTQ estão localizados
        fastq_dir="/media/ext5tb/anajulia/montagem2/fungi_reads/unmapped_fastq"
        #Diretório de saída para os arquivos FASTA
        output_dir="/media/ext5tb/anajulia/montagem2/fungi_reads/unmapped_fastq/fasta_files"
        
        #Criar diretório de saída, se não existir
        mkdir -p "$output_dir"
        
        #Loop sobre todos os arquivos FASTQ no diretório
        for fastq_file in "$fastq_dir"/*.fastq
        do
            # Extrair o nome base do arquivo (sem extensões .fastq)
            base_name=$(basename "$fastq_file" | sed 's/.fastq//')
        
            # Converter o arquivo FASTQ em FASTA usando seqtk
            echo "Convertendo $fastq_file para FASTA..."
            seqtk seq -A "$fastq_file" > "$output_dir"/"$base_name".fasta
        
            echo "$fastq_file convertido com sucesso para $output_dir/$base_name.fasta"
        done
        
        echo "Conversão concluída para todos os arquivos."


# Obtaining statistc information about the fasta files
seqkit stats CV39-NGS488_S35_L001_unmapped_R2.fasta (example)

# Checking sequencing bias using Infer Experiment (Galaxy)
Was used the bam file obtained from the aligment between sequencing reads and Vitis Vinifera genome (fasta file and gff file). The gff file was converted into bed file. Both of them were put in the Galaxy tool (Infer Experiment).

Results: This is PairEnd Data
Fraction of reads failed to determine: 0.0137
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0048
Fraction of reads explained by "1+-,1-+,2++,2--": 0.9815

The library is Reverse Forward 
![image](https://github.com/user-attachments/assets/a5c69b00-ba19-4fc0-afe8-5bcdc06eec68)
![image](https://github.com/user-attachments/assets/93f48502-a6ff-45e5-94c1-83329950b2ca)    

# Organizing files for Trinity
Trinity code has some specifics for the files, specially multiple files 
![image](https://github.com/user-attachments/assets/6b170023-2128-46bf-8ccf-0654f870b663) 

So I needed to put all the files in a table, and after, convert the table in a tab-delimited text file. I did this step for both of the datasets (interaction sequencing files and germinates spores sequencing files). 
For example: cond_A	cond_A_rep1	/media/ext5tb/anajulia/montagem2/fungi_reads/fungi_fasta/CV1-NGS488_S1_L001_cut_PE1.fasta	/media/ext5tb/anajulia/montagem2/fungi_reads/fungi_fasta/CV1-NGS488_S1_L001_cut_PE2.fasta
![image](https://github.com/user-attachments/assets/96fa3066-4506-41b8-bcd1-676479c23a6c) (line 1 in the table)

# Trinity Assembly
After organizing the files in the correct template for trinity, it is time for use it.

The assembly was conducted in some attempts, in order to test and verify some files
The first attempt used the wrong files, located in: /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_fasta (Estes arquivos não foram submetidos ao mapeamento pelo hisat)
THe code used was: docker run -v $(pwd):$(pwd) trinityrnaseq/trinityrnaseq Trinity \
--seqType fa --samples_file /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_fasta/trinitydata_interaction_renamed.txt \
--max_memory 150G --CPU 30 --output /media/ext5tb/anajulia/montagem2/trinity_output_renamed > trinity_renamed_run.log (run_trinity_renamed.sh)

The second attempt used the fastq files found in: /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/unmapped_fastq (Arquivos submetidos ao mapeamento para retirada de reads da videira)
The code used was: docker run --user $(id -u):$(id -g) -v $(pwd):$(pwd) trinityrnaseq/trinityrnaseq Trinity \
--seqType fq --samples_file /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/trinitydata_interaction_2.txt \
--max_memory 150G --CPU 40 --SS_lib_type RF --output /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/trinity_output_attempt5 > trinity_attempt5_run.log


CONDA - eval "$(/media/ext5tb/anajulia/miniconda3/bin/conda shell.bash hook)"

sudo: gustavoc
senha: genomics10,

# TrinityStats
################################
        ## Counts of transcripts, etc.
################################
Total trinity 'genes':	74185
Total trinity transcripts:	105628
Percent GC: 44.67

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 3196
	Contig N20: 2396
	Contig N30: 1920
	Contig N40: 1545
	Contig N50: 1215

	Median contig length: 400
	Average contig: 723.03
	Total assembled bases: 76372710

# BUSCO Analysis
        # BUSCO version is: 5.8.0 
        # The lineage dataset is: basidiomycota_odb10 (Creation date: 2020-09-10, number of genomes: 133, number of BUSCOs: 1764)
        # Summarized benchmarking in BUSCO notation for file /mnt/pulsar/files/staging/10326218/inputs/dataset_15dc61c8-4ee9-422d-9f9d-342227968d2e.dat
        # BUSCO was run in mode: euk_tran

	***** Results: *****

	C:85.0%[S:4.5%,D:80.4%],F:6.9%,M:8.2%,n:1764	   
	1499	Complete BUSCOs (C)			   
	80	Complete and single-copy BUSCOs (S)	   
	1419	Complete and duplicated BUSCOs (D)	   
	121	Fragmented BUSCOs (F)			   
	144	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched

 # Finding ORFs in transcripts
 The script biopython_orf_find.py was used, and it found 246618 orfs in the assembly file (output: orfs_montageminteracao.fasta)
 The script get_longestORF.py was used to select the longest orf in the orfs file, and it found 99361 orfs (output: longestorfs_montageminteracao.fasta)
 6267 transcripts without ORFs

 # Using the CDHits to colapse the assembly
cd-hit -i trinity_output_attempt5.Trinity.fasta -o montageminteracao_colapsed.fasta -T 20 -M 0 -c 0.9 -d 0

################################
	## Counts of transcripts, etc.
################################
Total trinity 'genes':	72600
Total trinity transcripts:	83533
Percent GC: 44.55

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 2794
	Contig N20: 2044
	Contig N30: 1591
	Contig N40: 1215
	Contig N50: 890

	Median contig length: 353
	Average contig: 605.20
	Total assembled bases: 50554586

	 # BUSCO version is: 5.8.0 
	# The lineage dataset is: basidiomycota_odb10 (Creation date: 2020-09-10, number of genomes: 133, number of BUSCOs: 1764)
	# Summarized benchmarking in BUSCO notation for file /mnt/pulsar/files/staging/10361100/inputs/dataset_796b3a49-62e3-4763-b764-50ffe9ffb893.dat
	# BUSCO was run in mode: euk_tran

	***** Results: *****

	C:84.8%[S:35.3%,D:49.5%],F:6.9%,M:8.4%,n:1764	   
	1495	Complete BUSCOs (C)			   
	622	Complete and single-copy BUSCOs (S)	   
	873	Complete and duplicated BUSCOs (D)	   
	121	Fragmented BUSCOs (F)			   
	148	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched











