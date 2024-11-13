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

# Organizing files for Trinity
Trinity code has some specifics for the files, specially multiple files 
![image](https://github.com/user-attachments/assets/6b170023-2128-46bf-8ccf-0654f870b663) 

So I needed to put all the files in a table, and after, convert the table in a tab-delimited text file. I did this step for both of the datasets (interaction sequencing files and germinates spores sequencing files). 
For example: cond_A	cond_A_rep1	/media/ext5tb/anajulia/montagem2/fungi_reads/fungi_fasta/CV1-NGS488_S1_L001_cut_PE1.fasta	/media/ext5tb/anajulia/montagem2/fungi_reads/fungi_fasta/CV1-NGS488_S1_L001_cut_PE2.fasta
![image](https://github.com/user-attachments/assets/96fa3066-4506-41b8-bcd1-676479c23a6c) (line 1 in the table)

# Trinity Assembly
After organizing the files in the correct template for trinity, it is time for use it.

Trinity --seqType fa --max_memory 15G --samples_file trinitydata_interaction.txt --CPU 30

Trinity --seqType fq --max_memory 15G --samples_file trinitydata_germinated.txt --CPU 30


CONDA - eval "$(/media/ext5tb/anajulia/miniconda3/bin/conda shell.bash hook)"













