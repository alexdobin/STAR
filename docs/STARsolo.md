**STARsolo**: mapping, demultiplexing and quantification for single cell RNA-seq
=================================================================================

Major updates in STAR 2.7.5a (2020/06/16)
---------------------------------------
* [**Smart-seq scRNA-seq process:**](#plate-based-Smart-seq-scRNA-seq)
    * STARsolo now supports for the plate-based (a.k.a. Smart-seq) scRNAs-seq technologies.

Major updates in STAR 2.7.3a (Oct 8 2019)
-----------------------------------------
* **Output enhancements:**
    * Summary.csv statistics output for raw and filtered cells useful for quick run quality assessment.
    * --soloCellFilter option for basic filtering of the cells, similar to the methods used by CellRanger 2.2.x.
* [**Better compatibility with CellRanger 3.x.x:**](#matching-cellranger-3xx-results)
    * --soloUMIfiltering MultiGeneUMI option introduced in CellRanger 3.x.x for filtering UMI collisions between different genes.
    * --soloCBmatchWLtype 1MM_multi_pseudocounts option, introduced in CellRanger 3.x.x, which slightly changes the posterior probability calculation for CB with 1 mismatch.
* [**Velocyto spliced/unspliced/ambiguous quantification:**](#velocyto-splicedunsplicedambiguous-quantification)
    * --soloFeatures Velocyto option to produce Spliced, Unspliced, and Ambiguous counts similar to the [velocyto.py](http://velocyto.org/) tool developed by [LaManno et al](https://doi.org/10.1038/s41586-018-0414-6). This option is under active development and the results may change in the future versions.
* [**Support for complex barcodes, e.g. inDrop:**](#barcode-geometry)
    * Complex barcodes in STARsolo with --soloType CB_UMI_Complex, --soloCBmatchWLtype --soloAdapterSequence, --soloAdapterMismatchesNmax, --soloCBposition,--soloUMIposition
* [**BAM tags:**](#bam-tags)
    * CB/UB for corrected CellBarcode/UMI
    * GX/GN for gene ID/name

STARsolo
-------------
STARsolo is a turnkey solution for analyzing droplet single cell RNA sequencing data (e.g. 10X Genomics Chromium System) built directly into STAR code.
STARsolo inputs the raw FASTQ reads files, and performs the following operations
* error correction and demultiplexing of cell barcodes using user-input whitelist
* mapping the reads to the reference genome using the standard STAR spliced read alignment algorithm
* error correction and collapsing (deduplication) of Unique Molecular Identifiers (UMIa)
* quantification of per-cell gene expression by counting the number of reads per gene
* quantification of other transcriptomic features: splice junctions; pre-mRNA; spliced/unspliced reads similar to Velocyto

STARsolo output is designed to be a drop-in replacement for 10X CellRanger gene quantification output.
It follows CellRanger logic for cell barcode whitelisting and UMI deduplication, and produces nearly identical gene counts in the same format.
At the same time STARsolo is ~10 times faster than the CellRanger.

Running STARsolo for 10X Chromium scRNA-seq data
-------------------------------------
* STARsolo is run the same way as normal STAR run, with addition of several STARsolo parameters:
    ```
   /path/to/STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  [...other parameters...] --soloType ... --soloCBwhitelist ...
    ```
   The genome index is the same as for normal STAR runs. </br>
   The parameters required to run STARsolo on 10X Chromium data are described below:

* The STAR solo algorithm is turned on with:
    ```
    --soloType Droplet
    ```
    or, since 2.7.3a, with more descriptive:
    ```
    --soloType CB_UMI_Simple
    ```    

* The CellBarcode whitelist has to be provided with:

    ```
    --soloCBwhitelist /path/to/cell/barcode/whitelist
    ```
    The 10X Chromium [whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-) files can be found or inside the CellRanger distribution or on [GitHub/10XGenomics](https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes).
Please make sure that the whitelist is compatible with the specific version of the 10X chemistry: V2 or V3. For instance, in CellRanger 3.1.0, the *V2 whitelist* is:
    ```
    cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt

    https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
    ```
    and *V3 whitelist* (gunzip it for STAR):
    ```
    cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

    https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
    ```
* The default barcode lengths (CB=16b, UMI=10b) work for 10X Chromium V2. For V3, specify:
    ```
    --soloUMIlen 12
    ```

* Importantly, in the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e.
  ```
  --readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz
  ```
  For instance, standard 10X runs have cDNA as Read2 and barcode as Read1:
  ```
  --readFilesIn Read2.fastq.gz Read1.fastq.gz
  ```
  For multiple lanes, use commas separated lists for Read2 and Read1:
  ```
  --readFilesIn Read2_Lane1.fastq.gz,Read2_Lane2.fastq.gz,Read2_Lane3.fastq.gz  Read1_Lane1.fastq.gz,Read1_Lane2.fastq.gz,Read1_Lane3.fastq.gz
  ```
----------------------------------------
How to make STARsolo _raw_ gene counts (almost) identical to CellRanger's
----------------------------------------------------------------------
* CellRanger uses its own "filtered" version of annotations (GTF file) which is a subset of ENSEMBL annotations, with several gene biotypes removed (mostly small non-coding RNA). Annotations affect the counts, and to match CellRanger counts CellRanger annotations have to be used.

* 10X provides several versions of the CellRanger annotations:</br>
  [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest ](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)</br>
  For the best match, the annotations in CellRanger run and STARsolo run should be exactly the same.

* The FASTA and GTF files
    ```
    refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
    refdata-cellranger-GRCh38-3.0.0/genes/genome.fa
    ```
    have to be used in STAR genome index generation step before mapping:
    ```
    STAR  --runMode genomeGenerate --runThreadN ... --genomeDir ./ --genomeFastaFiles /path/to/genome.fa  --sjdbGTFfile /path/to/genes.gtf
    ```

* If you want to use your own GTF (e.g. newer version of ENSEMBL or GENCODE), you can generate the "filtered" GTF file using 10X's mkref tool:</br>
  [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references)

* To make the agreement between STARsolo and CellRanger even more perfect, you can add
    ```
    --genomeSAsparseD 3
    ```
    to the genome generation options, which is used by CellRanger to generate STAR genomes. It will generate sparse suffixs array whic has an additional benefit of fitting into 16GB of RAM. However, it also results in 30-50% reduction of speed.

* The considerations above are for *raw* counts, i.e. when cell filtering is not performed. To get *filtered* results, refer to [Basic cell filtering](#basic-cell-filtering) section.

#### Matching CellRanger 3.x.x results
* By default, cell barcode and UMI collapsing parameters are designed to give the best agreement with CellRanger 2.x.x. CellRanger 3.x.x introduced some minor changes to this algorithm. To get the best agreement between STARsolo and CellRanger 3.x.x, add these parameters:
    ```
    --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts
    ```
-----------------
Barcode geometry
-------------------
* Simple barcode lengths and start positions on barcode reads are described with
  ```
  --soloCBstart, --soloCBlen, --soloUMIstart, --soloUMIlen
  ```
  which works with
  ```
  --soloType CB_UMI_Simple (a.k.a Droplet)
  ```
* More complex barcodes are activated with ```--soloType CB_UMI_complex``` and are described with the following parameters
  ```
  soloCBposition              -
    strings(s)              position of Cell Barcode(s) on the barcode read.
                            Presently only works with --soloType CB_UMI_Complex, and barcodes are assumed to be on Read2.
                            Format for each barcode: startAnchor_startDistance_endAnchor_endDistance
                            start(end)Anchor defines the anchor base for the CB: 0: read start; 1: read end; 2: adapter start; 3: adapter end
                            start(end)Distance is the distance from the CB start(end) to the Anchor base
                            String for different barcodes are separated by space.
                            Example: inDrop (Zilionis et al, Nat. Protocols, 2017):
                            --soloCBposition  0_0_2_-1  3_1_3_8

   soloUMIposition             -
    string                  position of the UMI on the barcode read, same as soloCBposition
                            Example: inDrop (Zilionis et al, Nat. Protocols, 2017):
                            --soloUMIposition  3_9_3_14

   soloAdapterSequence         -
    string:                 adapter sequence to anchor barcodes.

   soloAdapterMismatchesNmax   1
    int>0:                  maximum number of mismatches allowed in adapter sequence
    ```

--------------------------------------
Basic cell filtering
--------------------
* Since 2.7.3a, in addition to raw, unfiltered output of gene/cell counts, STARsolo performs simple (knee-like) filtering of the cells, similar to the methods used by CellRanger 2.2.x. This is turned on by default and is controlled by:
   ```
    soloCellFilter              CellRanger2.2 3000 0.99 10
        string(s):              cell filtering type and parameters
                            CellRanger2.2   ... simple filtering of CellRanger 2.2, followed by three numbers: number of expected cells, robust  maximum percentile for UMI count, maximum to minimum ratio for UMI count
                            TopCells        ... only report top cells by UMI count, followed by the exact number of cells
                            None            ... do not output filtered cells
    ```
* This filtering is used to produce summary statistics for filtered cells in the Summary.csv file, which is similar to CellRanger's summary and is useful for Quality Control.
* Recent versions of CellRanger switched to more advanced filtering done with the EmptyDrop tool developed by [Lun et al](https://doi.org/10.1186/s13059-019-1662-y). To obtain filtered counts similar to recent CellRanger versions, we need to run this tools on **raw** STARsolo output


---------------------------------------------------
Quantification of different transcriptomic features
---------------------------------------------------
* In addition to the gene counts (deafult), STARsolo can calculate counts for other transcriptomic features:
    * pre-mRNA counts, useful for single-nucleus RNA-seq. This counts all read that overlap gene loci, i.e. included both exonic and intronic reads:
        ```
        --soloFeatures GeneFull
        ```
    * Counts for annotated and novel splice junctions:
        ```
        --soloFeatures SJ
        ```
    * #### Velocyto spliced/unspliced/ambiguous quantification
        This option will calculate Spliced, Unspliced, and Ambiguous counts per cell per gene similar to the [velocyto.py](http://velocyto.org/) tool developed by [LaManno et al](https://doi.org/10.1038/s41586-018-0414-6). This option is under active development and the results may change in the future versions.
        ```
        --soloFeatures Gene Velocyto
        ```
      Note that Velocyto quantification requires Gene features
    * All the features can be conveniently quantified in one run:
        ```
        --soloFeatures Gene GeneFull SJ Velocyto
        ```

--------------------------------------
BAM tags
-----------------
* To output BAM tags into SAM/BAM file, add them to the list of standard tags in
    ```
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
    ```
    Any combinations of tags can be used.
* CR/UR: **raw (uncorrected)** CellBarcode/UMI
* CY/UY: quality score for CellBarcode/UMI
* GX/GN: for gene ID/names
* sS/sQ: for sequence/quality combined CellBarcode and UMI;  sM for barcode match status.
* CB/UB: **corrected** CellBarcode/UMI. Note, that these tags require sorted BAM output, i.e. we need to add:
    ```
    --outSAMtype BAM SortedByCoordinate
    ```

--------------------------------
Different scRNA-seq technologies
--------------------------------
### Plate-based (Smart-seq) scRNA-seq
Plate-based (Smart-seq) scRNA-seq technologies produce separate FASTQ files for each cell. Cell barcodes are not incorporated in the read sequences, and there are no UMIs. Typical STAR command for mapping and quantification of these file will look like:
```
--soloType SmartSeq --readFilesManifest /path/to/manifest.tsv --soloUMIdedup Exact --soloStrand Unstranded
```

* STARsolo `--soloType SmartSeq` option produces cell/gene (and other [features](#quantification-of-different-transcriptomic-features))
count matrices, using rules similar to the droplet-based technologies. The differnces are (i) individual cells correspond to different FASTQ files,there are no Cell Barcode sequences, and "Cell IDs" have to be provided as input (ii) there are no UMI sequences, but reads can be deduplicated if they have identical start/end coordinates.

* The convenient way to list all the FASTQ files and Cell IDs is to create a file manifest and supply it in `--readFilesManifest /path/to/manifest.tsv`. The manifest file should contain 3 tab-separated columns. For paired-end reads:
```
Read1-file-name \t Read2-file-name \t Cell-id
```
For single-end reads, the 2nd column should contain the dash - :
```
Read1-file-name \t - \t Cell-id
```
Cell-id can be any string without spaces. Cell-id will be added as ReadGroup tag (*RG:Z:*) for each read in the SAM/BAM output. If Cell-id starts with *ID:*, it can contain several fields separated by tab, and all the fields will be copied verbatim into SAM *@RG* header line.
* Deduplication based on read start/end coordinates can be done with `--soloUMIdedup Exact` option. To avoid deduplication (e.g. for single-end reads) use `--soloUMIdedup NoDedup`. Both deduplication options can be used together `--soloUMIdedup Exact NoDedup` and will produce two columns in the *matrix.mtx* output.
* Common Smart-seq protocols are unstranded and thus will require `--soloStrand Unstranded` option. If your protocol is stranded, you can can choose the proper `--soloStrand Forward` (default) or `--soloStrand Reverse` options.

-------------------------------------------------------------
------------------------------------------------------
--------------------------------------------------
All parameters that control STARsolo output are listed again below with defaults and short descriptions:
---------------------------------------
```
soloType                    None
    string(s): type of single-cell RNA-seq
                            CB_UMI_Simple   ... (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
                            CB_UMI_Complex  ... one UMI of fixed length, but multiple Cell Barcodes of varying length, as well as adapters sequences are allowed in read2 only, e.g. inDrop.
                            CB_samTagOut    ... output Cell Barcode as CR and/or CB SAm tag. No UMI counting. --readFilesIn cDNA_read1 [cDNA_read2 if paired-end] CellBarcode_read . Requires --outSAMtype BAM Unsorted [and/or SortedByCoordinate]
                            SmartSeq        ... Smart-seq: each cell in a separate FASTQ (paired- or single-end), barcodes are corresponding read-groups, no UMI sequences, alignments deduplicated according to alignment start and end (after extending soft-clipped bases)

soloCBwhitelist             -
    string(s): file(s) with whitelist(s) of cell barcodes. Only --soloType CB_UMI_Complex allows more than one whitelist file.
                            None            ... no whitelist: all cell barcodes are allowed

soloCBstart                 1
    int>0: cell barcode start base

soloCBlen                   16
    int>0: cell barcode length

soloUMIstart                17
    int>0: UMI start base

soloUMIlen                  10
    int>0: UMI length

soloBarcodeReadLength       1
    int: length of the barcode read
                            1   ... equal to sum of soloCBlen+soloUMIlen
                            0   ... not defined, do not check

soloCBposition              -
    strings(s)              position of Cell Barcode(s) on the barcode read.
                            Presently only works with --soloType CB_UMI_Complex, and barcodes are assumed to be on Read2.
                            Format for each barcode: startAnchor_startPosition_endAnchor_endPosition
                            start(end)Anchor defines the Anchor Base for the CB: 0: read start; 1: read end; 2: adapter start; 3: adapter end
                            start(end)Position is the 0-based position with of the CB start(end) with respect to the Anchor Base
                            String for different barcodes are separated by space.
                            Example: inDrop (Zilionis et al, Nat. Protocols, 2017):
                            --soloCBposition  0_0_2_-1  3_1_3_8

soloUMIposition             -
    string                  position of the UMI on the barcode read, same as soloCBposition
                            Example: inDrop (Zilionis et al, Nat. Protocols, 2017):
                            --soloCBposition  3_9_3_14

soloAdapterSequence         -
    string:                 adapter sequence to anchor barcodes.

soloAdapterMismatchesNmax   1
    int>0:                  maximum number of mismatches allowed in adapter sequence.

soloCBmatchWLtype           1MM_multi
    string:                 matching the Cell Barcodes to the WhiteList
                            Exact                   ... only exact matches allowed
                            1MM                     ... only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match.
                            1MM_multi               ... multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches.
                                                        Allowed CBs have to have at least one read with exact match. Similar to CellRanger 2.2.0
                            1MM_multi_pseudocounts  ... same as 1MM_Multi, but pseudocounts of 1 are added to all whitelist barcodes.
                                                        Similar to CellRanger 3.x.x

soloStrand                  Forward
    string: strandedness of the solo libraries:
                            Unstranded  ... no strand information
                            Forward     ... read strand same as the original RNA molecule
                            Reverse     ... read strand opposite to the original RNA molecule

soloFeatures                Gene
    string(s):              genomic features for which the UMI counts per Cell Barcode are collected
                            Gene            ... genes: reads match the gene transcript
                            SJ              ... splice junctions: reported in SJ.out.tab
                            GeneFull        ... full genes: count all reads overlapping genes' exons and introns

soloUMIdedup                1MM_All
    string(s):              type of UMI deduplication (collapsing) algorithm
                            1MM_All             ... all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)
                            1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).
                            Exact               ... only exactly matching UMIs are collapsed
                            NoDedup             ... no deduplication of UMIs, count all reads. Allowed for --soloType SmartSeq

soloUMIfiltering            -
    string(s)               type of UMI filtering
                            -               ... basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)
                            MultiGeneUMI    ... remove lower-count UMIs that map to more than one gene (introduced in CellRanger 3.x.x)

soloOutFileNames            Solo.out/          features.tsv barcodes.tsv        matrix.mtx
    string(s)               file names for STARsolo output:
                            file_name_prefix   gene_names   barcode_sequences   cell_feature_count_matrix

soloCellFilter              CellRanger2.2 3000 0.99 10
    string(s):              cell filtering type and parameters
                            CellRanger2.2   ... simple filtering of CellRanger 2.2, followed by three numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count
                            TopCells        ... only report top cells by UMI count, followed by the exact number of cells
                            None            ... do not output filtered cells

soloOutFormatFeaturesGeneField3 "Gene Expression"
        string(s):                              field 3 in the Gene features.tsv file. If "-", then no 3rd field is output.
```
