# TCC 🍇 💻 🧬

# fastqc
 automatic script (run_fastqc.sh)
	 #!/bin/bash
	mkdir -p ./FastQC_raw
	for i in *fastq.gz
	do
	fastqc "$i" -o ./FastQC_raw -t 8
	done

# multiqc
	 multiqc /media/ext5tb/anajulia/montagem2/raw_reads_data/FastQC_raw

# cutadapt
First we needed to find the adapters used is the sequencing. This was the kit used TruSeq® Stranded mRNA Library Prep kit e o sequenciador NovaSeq 6000. So we entered the Illumina site and found the sequences used as adapters.
https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html
https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/CDIndexes.htm -> ILLUMINA SITE
TruSeq DNA and RNA CD Indexes -> SEQUENCING TYPE

The following sequences are used for adapter trimming.

Read 1

AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Read 2

AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

The script below was used (cutadapt.sh)
	
 	#!/bin/bash
	
	for i in /media/ext5tb/anajulia/montagem2/raw_reads_data/*_R1_001.fastq.gz
	do
	    file2=$(echo "$i" | sed "s/_R1/_R2/g")
	    output1=$(basename "$i" | sed "s/_R1_001.fastq.gz//g")
	
	    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	    -m 50 --max-n 2 -q 20 -j 6 \
	    -o "${output1}_cut_PE1.fastq.gz" \
	    -p "${output1}_cut_PE2.fastq.gz" \
	    "$i" "$file2" > "${output1}_cutadapt_summary.txt"
	
	done


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

# Trinity Assembly - Plant-pathogen interaction
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
Counts of transcripts, etc.
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
BUSCO version is: 5.8.0 
The lineage dataset is: basidiomycota_odb10 (Creation date: 2020-09-10, number of genomes: 133, number of BUSCOs: 1764)
Summarized benchmarking in BUSCO notation for file /mnt/pulsar/files/staging/10326218/inputs/dataset_15dc61c8-4ee9-422d-9f9d-342227968d2e.dat
BUSCO was run in mode: euk_tran

	***** Results: *****

	C:85.0%[S:4.5%,D:80.4%],F:6.9%,M:8.2%,n:1764	   
	1499	Complete BUSCOs (C)			   
	80	Complete and single-copy BUSCOs (S)	   
	1419	Complete and duplicated BUSCOs (D)	   
	121	Fragmented BUSCOs (F)			   
	144	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched

 # Assemblying with different kmers sizes (WRONG!!!)
 In order to obtain the best assembly, the assembly will be performed again, but this this time, using 2 different kmers sizes. The same code will be used, and only the parameter "--KMER_SIZE" will be add. After the assembly with Trinity, the CD-Hit tool will be used to collapse the two assemblies, removing redundancy and clustering similar sequences to generate a non-redundant dataset for downstream analyses.

Assembly 1: docker run --user $(id -u):$(id -g) -v $(pwd):$(pwd) trinityrnaseq/trinityrnaseq Trinity \
--seqType fq --samples_file /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/trinitydata_interaction_2.txt \
--max_memory 150G --CPU 40 --SS_lib_type RF --min_kmer_cov 25 --output /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/trinity_output_kmer25 > trinity_kmer25_run.log

Total trinity 'genes':	10505
Total trinity transcripts:	16502
Percent GC: 45.62

	Contig N10: 2117
	Contig N20: 1727
	Contig N30: 1449
	Contig N40: 1213
	Contig N50: 1028

	Median contig length: 521
	Average contig: 719.06
	Total assembled bases: 11865985

 C:38.0%[S:4.4%,D:33.6%],F:9.5%,M:52.5%,n:1764	   
	670	Complete BUSCOs (C)			   
	78	Complete and single-copy BUSCOs (S)	   
	592	Complete and duplicated BUSCOs (D)	   
	168	Fragmented BUSCOs (F)			   
	926	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched		   


Assembly 2: docker run --user $(id -u):$(id -g) -v $(pwd):$(pwd) trinityrnaseq/trinityrnaseq Trinity \
--seqType fq --samples_file /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/trinitydata_interaction_2.txt \
--max_memory 150G --CPU 40 --SS_lib_type RF --min_kmer_cov 31 --output /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/trinity_output_kmer31 > trinity_kmer31_run.log

Total trinity 'genes':	8654
Total trinity transcripts:	13785
Percent GC: 45.79

	Contig N10: 2069
	Contig N20: 1654
	Contig N30: 1393
	Contig N40: 1173
	Contig N50: 996

	Median contig length: 521
	Average contig: 708.11
	Total assembled bases: 9761362

 C:33.6%[S:3.5%,D:30.2%],F:8.3%,M:58.0%,n:1764	   
	593	Complete BUSCOs (C)			   
	61	Complete and single-copy BUSCOs (S)	   
	532	Complete and duplicated BUSCOs (D)	   
	147	Fragmented BUSCOs (F)			   
	1024	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched

 Then we need to combine these two assemblies:
 cat trinity_output_kmer25.Trinity.fasta trinity_output_kmer31.Trinity.fasta > combined_kmers_interaction.fasta

 And then, we used the cd-hits-est to remove redundance in the sequences:
 cd-hit-est -i combined_kmers_interaction.fasta -o collapsed_kmers_interaction.fasta -c 0.90 -n 8 -T 25 -M 100000

 
 # Finding ORFs in transcripts
 The script biopython_orf_find.py was used, and it found 246618 orfs in the assembly file (output: orfs_montageminteracao.fasta)
 The script get_longestORF.py was used to select the longest orf in the orfs file, and it found 99361 orfs (output: longestorfs_montageminteracao.fasta)
 6267 transcripts without ORFs

 # Using the CDHits to colapse the assembly
cd-hit -i trinity_output_attempt5.Trinity.fasta -o montageminteracao_colapsed.fasta -T 20 -M 0 -c 0.9 -d 0

TrinityStats
Counts of transcripts, etc.
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

BUSCO version is: 5.8.0 
The lineage dataset is: basidiomycota_odb10 (Creation date: 2020-09-10, number of genomes: 133, number of BUSCOs: 1764)
Summarized benchmarking in BUSCO notation for file /mnt/pulsar/files/staging/10361100/inputs/dataset_796b3a49-62e3-4763-b764-50ffe9ffb893.dat
BUSCO was run in mode: euk_tran

	***** Results: *****

	C:84.8%[S:35.3%,D:49.5%],F:6.9%,M:8.4%,n:1764	   
	1495	Complete BUSCOs (C)			   
	622	Complete and single-copy BUSCOs (S)	   
	873	Complete and duplicated BUSCOs (D)	   
	121	Fragmented BUSCOs (F)			   
	148	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched

 # Finding ORFs in the colapsed assembly
 The script biopython_orf_find.py was used, and it found 187800 orfs in the assembly file (output: orfs_montagemcolapsada.fasta)
 The script get_longestORF.py was used to select the longest orf in the orfs file, and it found 77737 orfs (output: longestorfs_montagemcolapsada.fasta)
 5796 transcripts without ORFs

# Trinity Assembly - Germinated Spores
docker run --user $(id -u):$(id -g) -v /media/ext5tb/anajulia/montagem2/fungo_germinado:/media/ext5tb/anajulia/montagem2/fungo_germinado trinityrnaseq/trinityrnaseq Trinity \
--seqType fq \
--left /media/ext5tb/anajulia/montagem2/fungo_germinado/EG1_NGS396_S8_L001_R1_001.fastq.gz_cut_R1.fastq,/media/ext5tb/anajulia/montagem2/fungo_germinado/EG2_NGS396_S8_L001_R1_001.fastq.gz_cut_R1.fastq,/media/ext5tb/anajulia/montagem2/fungo_germinado/EG3_NGS396_S8_L001_R1_001.fastq.gz_cut_R1.fastq \
--right /media/ext5tb/anajulia/montagem2/fungo_germinado/EG1_NGS396_S8_L001_R1_001.fastq.gz_cut_R2.fastq,/media/ext5tb/anajulia/montagem2/fungo_germinado/EG2_NGS396_S8_L001_R1_001.fastq.gz_cut_R2.fastq,/media/ext5tb/anajulia/montagem2/fungo_germinado/EG3_NGS396_S8_L001_R1_001.fastq.gz_cut_R2.fastq \
--max_memory 150G --CPU 40 --SS_lib_type RF \
--output /media/ext5tb/anajulia/montagem2/fungo_germinado/trinity_montagem_germinado > trinity_germinado_run.log

# Trinity Stats
Counts of transcripts, etc.
Total trinity 'genes':	33658
Total trinity transcripts:	67520
Percent GC: 43.97

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 4278
	Contig N20: 3254
	Contig N30: 2698
	Contig N40: 2288
	Contig N50: 1946

	Median contig length: 784
	Average contig: 1183.57
	Total assembled bases: 79914310

 # BUSCO
BUSCO version is: 5.8.0 
The lineage dataset is: basidiomycota_odb10 (Creation date: 2020-09-10, number of genomes: 133, number of BUSCOs: 1764)
Summarized benchmarking in BUSCO notation for file /mnt/pulsar/files/staging/10379549/inputs/dataset_9f881ead-aa1b-4897-bce4-a97cf1361acc.dat
BUSCO was run in mode: euk_tran

	***** Results: *****

	C:90.7%[S:3.6%,D:87.1%],F:2.8%,M:6.5%,n:1764	   
	1600	Complete BUSCOs (C)			   
	63	Complete and single-copy BUSCOs (S)	   
	1537	Complete and duplicated BUSCOs (D)	   
	50	Fragmented BUSCOs (F)			   
	114	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched
 

# Mapping the transcriptomes using blast (WRONG!!!!)
The blastn was used to map the two assemblies, in order to identify sequences there are only expressed during interaction with the plant. So, first of all, the blast databank was created using the spore germinated transcriptome, because this assembly will be the reference for the blastn.

makeblastdb -in trinity_montagem_germinado.Trinity.fasta -dbtype nucl -out germinado_db

After this, we ran the blastn, using the interaction transcriptome to compare with the databank created early.

blastn -query montageminteracao_colapsed.fasta -db germinado_db -out germinado_vs_interacao.blastn -evalue 1e-05 -outfmt "6 std qcovs" -word_size 6 -num_threads 30

Best Hits: 30 (identidade ) e 80 (cobertura) 

The script below was used to select the best hits in the blast file: 
	#!/usr/bin/python3
	
	# Parâmetros de corte
	ident_cutoff = 30.0
	qcov_cutoff = 80.0
	
	# Dicionário para armazenar os melhores hits
	dic = {}
	
	# Leitura do arquivo de entrada
	with open("germinado_vs_interacao.blastn") as inp:
	    for line in inp:
	        if line[0] != '#':  # Ignorar linhas de cabeçalho
	            line = line.strip().split('\t')
	            try:
	                # line[2] = %ident, line[12] = %qcoverage
	                if float(line[2]) >= ident_cutoff and float(line[12]) >= qcov_cutoff:
	                    # Atualizar ou adicionar o melhor hit baseado no %ident
	                    if line[0] not in dic or float(dic[line[0]][2]) < float(line[2]):
	                        dic[line[0]] = line
	            except (KeyError, IndexError, ValueError):
	                continue
	
	# Escrita dos resultados no arquivo de saída
	with open('./germinado_vs_interacao_best_hits.tsv', 'w') as out:
	    # Cabeçalho do arquivo de saída
	    out.write('query\tsubject\t%ident\talignment_length\tmismatches\tgap_opens\tqstart\tqend\tsstart\tsend\tevalue\tbit_score\t%qcoverage\n')
	    for hit in dic.values():
	        out.write('\t'.join(hit) + '\n')

  After, a table file was done using the information in the best hits file, in order to register the corresponding IDs:
 	 #!/usr/bin/python3

	# Arquivos de entrada e saída
	input_file = "germinado_vs_interacao_best_hits.tsv"
	output_file = "correspondencias_ids.tsv"
	
	# Abrir o arquivo de entrada e criar o de saída
	with open(input_file) as inp, open(output_file, "w") as out:
	    # Escrever o cabeçalho no arquivo de saída
	    out.write("query_id subject_id\n")
	
	    # Ignorar o cabeçalho do arquivo de entrada
	    next(inp)
	
	    # Processar cada linha do arquivo
	    for line in inp:
	        # Dividir a linha por qualquer espaço ou tabulação
	        fields = line.strip().split()
	
	        if len(fields) >= 2:  # Garantir que há pelo menos duas colunas
	            query_id = fields[0]
	            subject_id = fields[1]
	
	            # Escrever os IDs no arquivo de saída separados por espaço
	            out.write(f"{query_id} {subject_id}\n")


# Quantification using SALMON

1- Changing the files names: /media/ext5tb/anajulia/montagem2/fungi_reads/fungi_cut_fastq/unmapped_fastq
The files were renamed according to their specific treatments and time points, as shown in the examples: inoc_0hpi_40_R1_PE1.fastq, inoc_0hpi_80_R1_PE2.fastq, inoc_7dpi_40_R1_PE1.fastq, inoc_15dpi_40_R2_PE2.fastq.

The same was done to the files located in: /media/ext5tb/anajulia/montagem2/fungo_germinado
Examples:  invitro_R1_PE1.fastq, invitro_R3_PE1.fastq, invitro_R2_PE2.fastq

2 - A table file was organized in order to give the exact file locations for Trinity. The table model is specified in topic Trinity Assembly

3 - The script to run Trinity was modified using the new files:
docker run --user $(id -u):$(id -g) -v $(pwd):$(pwd) trinityrnaseq/trinityrnaseq Trinity \
--seqType fq --samples_file /media/ext5tb/anajulia/montagem2/salmon/transcript_ref/invitro_vs_inoc_data.txt \
--max_memory 150G --CPU 40 --SS_lib_type RF --output /media/ext5tb/anajulia/montagem2/salmon/transcript_ref/trinity_ref_invitro_inoc > trinity_invitro_inoc_run.log

4 - Trinity Stats results
Total trinity 'genes':	90836
Total trinity transcripts:	138217
Percent GC: 44.29

	Contig N10: 3771
	Contig N20: 2772
	Contig N30: 2222
	Contig N40: 1823
	Contig N50: 1445

	Median contig length: 420
	Average contig: 805.38
	Total assembled bases: 111316687

5 - BUSCO Results 
C:91.5%[S:4.1%,D:87.4%],F:3.1%,M:5.4%,n:1764	   
	1614	Complete BUSCOs (C)			   
	73	Complete and single-copy BUSCOs (S)	   
	1541	Complete and duplicated BUSCOs (D)	   
	54	Fragmented BUSCOs (F)			   
	96	Missing BUSCOs (M)			   
	1764	Total BUSCO groups searched

 6 - Collapsing with CD-Hits
 Done in Galaxy. Results were not significant, so we chose to use the transcriptome done without CD-Hits.

7 - Run Salmon
	#!/bin/bash
	
	ref=/media/ext5tb/anajulia/montagem2/salmon/transcript_ref/trinity_ref_invitro_inoc.Trinity.fasta.gz 
	index=/media/ext5tb/anajulia/montagem2/salmon/transcript_ref/salmon_index_invitro_inoc 
	
	# Index reference
	/media/SSD1TB/pedro/salmon-latest_linux_x86_64/bin/salmon index -t "$ref" -p 20 -i "$index"
	
	# Quantify against reference
	for i in /media/ext5tb/anajulia/montagem2/salmon/invitro/*PE1.fastq 
	do
	    file2=$(echo "$i" | sed "s/PE1/PE2/g")
	    rep=$(basename "$i" | sed -E "s/invitro_R([0-9]+)_.*/\1/g") 
	
	    echo "Analyzing $i and $file2 for replicate $rep"
	
	/media/SSD1TB/pedro/salmon-latest_linux_x86_64/bin/salmon quant -i "$index" \
	-l A -1 "$i" -2 "$file2" \
	--validateMappings \
	 --threads 10 \
	 --seqBias \
	 --gcBias \
	 --minAssignedFrags 1 \
	-o /media/ext5tb/anajulia/montagem2/salmon/invitro/count_invitro_"$rep" 
	done

8 - Merging files
R scripts were used to merge the quant files given by Salmon, ir order to organize the multiple files in only one main file. In this way, it will be possible do identify the transcripts counts in each os the treatments. 

9 - Filter the transcripts that are expressed only in plant, not in in vitro condition. For this step, we used two R scripts with two differente filter conditions.
The first one, we separated the transcripts with 0 counts in in vitro condition and at least 1 positive one treatment must have a count ≥ 1 in every replicate (Script filter_inplanta_expressed_soft).

The second one, the transcript must have a count of 0 in all replicates and must have a positive count (≥ 1) in ALL replicates (Script filter_inplanta_expressed_hard)

# We discovered that there were still plant sequences in our files, which were interfering with the analyses of the fungal transcripts. Therefore, we had to adopt a different strategy. The sequencing reads, which had already been mapped against the grapevine genome (*Vitis labrusca* - var. Concord), were analyzed using Kraken2 with a plant database. Thus, the reads that showed no similarity to the database, meaning the unclassified reads, were selected again. The processes followed the same steps as described in the previous section.

Kraken2 - conda activate kraken2
kraken2_db - plant library and taxonomy
Krona - kraken2 reports combined 

# Finding Effector Candidates



 








