STARsolo: mapping, demultiplexing and gene quantification for single cell RNA-seq
---------------------------------------------------------------------------------

First released in STAR 2.7.0a (Jan 23 2019)

STARsolo is a turnkey solution for analyzing droplet single cell RNA sequencing data (e.g. 10X Genomics Chromium System) built directly into STAR code.
STARsolo inputs the raw FASTQ reads files, and performs the following operations 
* error correction and demultiplexing of cell barcodes using user-input whitelist
* mapping the reads to the reference genome using the standard STAR spliced read alignment algorithm
* error correction and collapsing (deduplication) of Unique Molecular Identifiers (UMIa) 
* quantification of per-cell gene expression by counting the number of reads per gene

STARsolo output is designed to be a drop-in replacement for 10X CellRanger gene quantification output.
It follows CellRanger logic for cell barcode whitelisting and UMI deduplication, and produces nearly identical gene counts in the same format.
At the same time STARsolo is ~10 times faster than the CellRanger.

The STAR solo algorithm is turned on with:
```
--soloType Droplet 
```

Presently, the cell barcode whitelist has to be provided with:
``` 
--soloCBwhitelist /path/to/cell/barcode/whitelist
```

The 10X Chromium whitelist file can be found inside the CellRanger distribution, 
e.g. [10X-whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-).
Please make sure that the whitelist is compatible with the specific version of the 10X chemistry (V1,V2,V3 etc).

Importantly, in the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e.
```
--readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz
```

Important: the genome index has to be re-generated with the latest 2.7.0x release.

Other parameters that control STARsolo output are listed below. Note that default parameters are compatible with 10X Chromium V2 protocol.

```
soloCBstart                 1
    int>0: cell barcode start base

soloCBlen                   16
    int>0: cell barcode length

soloUMIstart                17
    int>0: UMI start base

soloUMIlen                  10
    int>0: UMI length

soloStrand                  Forward
    string: strandedness of the solo libraries:
                            Unstranded  ... no strand information
                            Forward     ... read strand same as the original RNA molecule
                            Reverse     ... read strand opposite to the original RNA molecule

soloFeatures                Gene
    string(s)               genomic features for which the UMI counts per Cell Barcode are collected
                            Gene            ... genes: reads match the gene transcript
                            SJ              ... splice junctions: reported in SJ.out.tab

soloUMIdedup                1MM_All
    string(s)               type of UMI deduplication (collapsing) algorithm
                            1MM_All             ... all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)
                            1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).
                            1MM_NotCollapsed      ... UMIs with 1 mismatch distance to others are not collapsed (i.e. all counted)

soloOutFileNames            Solo.out/ genes.tsv barcodes.tsv matrix.mtx matrixSJ.mtx
    string(s)               file names for STARsolo output
                            1st word    ... file name prefix
                            2nd word    ... barcode sequences
                            3rd word    ... gene IDs and names
                            4th word    ... cell/gene counts matrix
                            5th word    ... cell/splice junction counts matrix
```
