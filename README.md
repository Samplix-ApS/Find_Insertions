# Find Insertions
## Table of Contents
- [Description](#descript_)
 - [Requirements](#reqs_)
- [Find insertions](#find_insert_)
  - [Mask genome](#mask_genome_)
  - [Create custom genome](#custom_genome_)
  - [Identify reads of interest](#ident_reads_)
    - [Strict approach](#strict_)
    - [Loose approach](#loose_)
  - [Map extracted reads](#map_extract_)
  - [OPTIONAL: Karyoplot](#opt_plot_)
  - [Find mapped reads](#find_reads_)
  - [Create custom genome with insert](#insert_genome_)
  - [Map to custom genome with insert](#map_insert_genome_)
  - [Multiple insertions](#multi_insert_)
- [Authors](#authors_)

# <a name="descript_"></a> Description
This manual describes how to find insertions in the alignment.
## <a name="reqs_"></a> Requirements
* samtools
* seqkit
* python3
* addSAMtag.py - [found here](https://github.com/Samplix-ApS/add_SAM_tag)
* prep_reference.py - [found here](https://github.com/Samplix-ApS/Prepare_reference_for_pipeline)
* karyoploteR_read_density.R - [found here](https://github.com/Samplix-ApS/karyoploteR_read_density)

# <a name="find_insert_"></a> Find insertions
## <a name="mask_genome_"></a> Mask genome
Start out by masking regions in the genome where the construct maps.
Map the construct to the genome by using the following commands:
```
/usr/local/minimap/minimap2 -ax map-ont Genome.fasta.gz CONSTRUCT.fasta > CONSTRUCT_TO_GENOME.sam 
/home/ec2-user/samtools-1.9/samtools sort -@ 16 -O BAM -o  CONSTRUCT_TO_GENOME.bam  CONSTRUCT_TO_GENOME.sam 
/home/ec2-user/samtools-1.9/samtools index -@ 16 CONSTRUCT_TO_GENOME.bam
```

Then create a bed file for masking:
```
/usr/local/bedtools/bin/bedtools bamtobed -i CONSTRUCT_TO_GENOME.bam > CONSTRUCT_TO_GENOME.bed
```

Mask the genome based on construct mapping positions. Input must be in fasta format. ```-mc Specify a masking character. default N. ```

```
bedtools maskfasta [OPTIONS] -fi Genome.fasta -bed CONSTRUCT_TO_GENOME.bed -fo masked_Genome.fasta
```

## <a name="custom_genome_"></a> Create custom genome
Generate a costum genome containg the construct. Sequences need to be in same format e.g. either fa or fasta.gz (```gzip``` to compress or ```gzip -d``` to decompress).
Run the following command:
```
cat masked_Genome.fasta CONSTRUCT.fasta  > Masked_Genome_with_constructs.fasta
```

## <a name="ident_reads_"></a> Identify reads of interest
There are two diffirent options for identifying the reads of interest: strict or loose.
The strict approach is to map to the custom genome containing the construct and extract the reads mapping to the construct. The loose approach is to map directly to the construct and extract the aligned reads. The loose setting will give more false-positive integration sites, but also more reads on the potential good integrations sites. 
### <a name="strict_"></a> Strict approach
Map reads to masked custom genome with construct by using the following command:
```
/usr/local/minimap/minimap2 -ax map-ont Masked_Genome_with_constructs.fasta DATA.fastq > DATA_MAPPED_TO_MASKED_GENOME_CONSTRUCT.sam
/home/ec2-user/samtools-1.9/samtools sort -@ 16 -O BAM -o DATA_MAPPED_TO_MASKED_GENOME_CONSTRUCT.sort.bam DATA_MAPPED_TO_MASKED_GENOME_CONSTRUCT.sam
/home/ec2-user/samtools-1.9/samtools index DATA_MAPPED_TO_MASKED_GENOME_CONSTRUCT.sort.bam

```
Take a look in IGV and look at coverage profile in order to identify start and end positions.

Then extract reads aligned to construct:
´´´
samtools view -@ 16 -hb DATA_MAPPED_TO_MASKED_GENOME_CONSTRUCT.sort.bam CHR:start-end > DATA_MAPPED_TO_CONSTRUCT_xx-yy.bam
samtools view DATA_MAPPED_TO_CONSTRUCT_xx-yy.bam | awk '{print $1}' | sort | uniq -c | awk '{print $2}' > read.names
seqkit grep -f read.names DATA.fastq > READS_MAPPED_TO_CONSTRUCT_xx-yy.fastq
seqkit seq -n READS_FROM_INSERTION_REGION.fastq | wc -l
´´´
Note down the number of reads that mapped.

### <a name="loose_"></a> Loose approach
Map reads directly to the construct:
```
/usr/local/minimap/minimap2 -ax map-ont CONSTRUCT.fasta DATA.fastq > DATA_MAPPED_TO_CONSTRUCT.sam
/home/ec2-user/samtools-1.9/samtools sort -@ 16 -O BAM -o DATA_MAPPED_TO_CONSTRUCT.sort.bam DATA_MAPPED_TO_CONSTRUCT.sam
/home/ec2-user/samtools-1.9/samtools index DATA_MAPPED_TO_CONSTRUCT.sort.bam
```
Take a look in IGV and look at coverage profile in order to identify start and end positions.
Then extract reads aligned to construct:
```
/home/ec2-user/samtools-1.9/samtools view -@ 16 -hb [DATA_MAPPED_TO_CONSTRUCT.sort.bam OR DATA_MAPPED_TO_MASKED_GENOME_CONSTRUCT.sort.bam] CONSTRUCT_NAME:xx-yy > DATA_MAPPED_TO_CONSTRUCT_xx-yy.bam
/home/ec2-user/samtools-1.9/samtools bam2fq -@ 16 DATA_MAPPED_TO_CONSTRUCT_xx-yy.bam > READS_MAPPED_TO_CONSTRUCT_xx-yy.fastq
```
Note down the number of reads that mapped.

## <a name="map_extract_"></a> Map extracted reads
Map the extracted reads to the custom masked genome containing the construct.
Use the following command:
```
/usr/local/minimap/minimap2 -ax map-ont Masked_Genome_with_construct.fasta READS_MAPPED_TO_CONSTRUCT_xx-yy.fastq > DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.sam 
```
## <a name="add_tag_"></a> Add tags
Add tags to the alignment file, which will help identify in IGV all the chromosomes a read aligns to.
```
samtools sort -@ 16 -O BAM -o  DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.sort.bam DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.sam 
samtools index -@ 16 DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.sort.bam
python3 /home/ec2-user/data_CAJ/python_scripts/addSAMtag.py -i BAM/SAM_file <optional> 
```
## <a name="opt_plot_"></a> OPTIONAL: Karyoplot
This is optional. Plot ```DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_GENOME``` on chromosomes. It may be best to use a genome where the construct is not present, as the x-axis scale is set to highest density.
```
python3 /home/ec2-user/data_CAJ/python_scripts/prep_reference.py -i REFERENCE <optional>
Rscript /home/ec2-user/data_CAJ/python_scripts/karyoploteR_read_density.R -i DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.TAG.SORT.bam -p /home/ec2-user/Refseq/Standard/REFERENCE_GENOME/karyoploteR.genome.file.txt
```
## <a name="find_reads_"></a> Find mapped reads
Find where reads map in ```DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.TAG.SORT.bam```.
Find all positions with 1 coverage, merge continous regions within 100 bp. See regions using coverage, and the full command is:
```
bedtools genomecov -bg -ibam DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.TAG.SORT.bam | awk '$4>=1' | bedtools merge -d 100 -i stdin | bedtools coverage -a stdin -b DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.TAG.SORT.bam > Potentil_insert_sites.txt
```
Here for example there are 21 reads in region chr1:154566162-154566294. These 21 reads cover 132bp and the region itself is 132bp long, so 100% of the region is covered.
```
chr1	154566162	154566294	21	132	132	1.0000000
```
Use the different plots and positions to look at regions with coverage in ```DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.TAG.SORT.bam``` in IGV.

## <a name="insert_genome_"></a> Create custom genome with insert
Based on the results, create a new reference genome with the construct inserted. This can be done in CLC mainwork bench. 
Open the genome of interest, and rename. Find  the position for insert, right click and choose edit seletion. Add the insert and e.g. delete bases.
Export the new genome as a fasta sequence and transfer to the server.

## <a name="map_insert_genome_"></a> Map to custom genome with insert
Map all data to the new custom genome containing the insert either manually or use the pipeline.

For manual mapping, use the following:
```
/usr/local/minimap/minimap2 -ax map-ont NEW_GENOME.fasta DATA.fastq > DATA_MAPPED_TO_NEW_GENOME.sam
/home/ec2-user/samtools-1.9/samtools sort -@ 16 -O BAM -o DATA_MAPPED_TO_NEW_GENOME.sort.bam DATA_MAPPED_TO_NEW_GENOME.sam
/home/ec2-user/samtools-1.9/samtools index DATA_MAPPED_TO_NEW_GENOME.sort.bam
```
Take a look in IGV and look at coverage profile. Adjust the new custom genome with insert file if needed.  

## <a name="multi_insert_"></a> Multiple insertions
In case there are multiple insertion sites expected, it may be necessary to extract reads from each potential integration site and map again to the genome containing the construct.

Extract reads in region with coverage:
```
samtools view -@ 16 -hb DATA_MAPPED_TO_CONSTRUCT_xx-yy_TO_MASKED_GENOME_CONSTRUCT.TAG.SORT.bam CHR:INSERTION_REGION > READS_FROM_INSERTION_REGION.bam
samtools view READS_FROM_INSERTION_REGION.bam | awk '{print $1}' | sort | uniq -c | awk '{print $2}' > read.names
seqkit grep -f read.names DATA_MAPPED_TO_CONSTRUCT_xx-yy.fastq > READS_FROM_INSERTION_REGION.fastq
seqkit seq -n READS_FROM_INSERTION_REGION.fastq | wc -l
```
Note down the number of reads that mapped 

Then map to genome with construct:
```
/usr/local/minimap/minimap2 -ax map-ont Masked_Genome_with_constructs.fasta READS_FROM_INSERTION_REGION.fastq > READS_FROM_INSERTION_REGION_GENOME_WITH_CONSTRUCT.sam 
/home/ec2-user/samtools-1.9/samtools sort -@ 16 -O BAM -o  READS_FROM_INSERTION_REGION_GENOME_WITH_CONSTRUCT.sort.bam  READS_FROM_INSERTION_REGION_GENOME_WITH_CONSTRUCT.sam  
/home/ec2-user/samtools-1.9/samtools index -@ 16 READS_FROM_INSERTION_REGION_GENOME_WITH_CONSTRUCT.sort.bam
```

Then add tag as previously.

# <a name="authors_"></a> Authors
SOP from Services provided by CNA and adapted by CAJ
