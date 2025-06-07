# PDRF_PCR_Analysis
Script for analysis DNA seqeunces files types .fna and find to sites of the restrictions enzymes. By this script you can desing development the PDRF PCR for more quantity files with homologous genes.

Script main commands:

analyze Analyze files

  PDRF_PCR_Analysis analyze C:\Gene_files

enzymes List of default restrictase

  PDRF_PCR_Analysis enzymes

--help show this help message and exit

  PDRF_PCR_Analysis -- help

Script options commands:

-o OUTDIR, --outdir OUTDIR
                        Output path file for save data table. For save results default used directory that contains
                        script file
  -rs RESTRICTASE [RESTRICTASE ...], --restrictase RESTRICTASE [RESTRICTASE ...]
                        By default, it contains an internal dictionary of restriction sites. The list of restriction
                        sites can be viewed using the command list sites
  -ff FILTER_FILES [FILTER_FILES ...], --filter_files FILTER_FILES [FILTER_FILES ...]
                        Filter for name or type file, default .fna. Example: Seq .fna. You will analysis files that
                        contains in the name Seq and type .fna: Seq234.fna, Seq234.fna. Files name's 123.fna, Seq.gbk
                        not alnalysis, or others files not contains in the name Seq and .fna
  -fcs FILTER_COUNT_SITES [FILTER_COUNT_SITES ...], --filter_count_sites FILTER_COUNT_SITES [FILTER_COUNT_SITES ...]
                        Filter count sites of enzyme in each sequence. Exmaple 0 0. First number is count sites
                        enzymes in analyzed sequences, second number is max count sites positions restriction enzymes
                        in sequence.
  -fr FILTER_RESTRICTASE [FILTER_RESTRICTASE ...], --filter_restrictase FILTER_RESTRICTASE [FILTER_RESTRICTASE ...]


  Example:

  We have a folder with .fna files. sequences of genes of different families. And we need to find a suitable enzyme for rectification analysis with which we can analyze gene families.

  ![image](https://github.com/user-attachments/assets/bdf40760-ef2f-46c3-b1d6-f3ecb0c9f458)

  PDRF_PCR_Analysis analyze "D:\06_Python\py\Cry_genes_seq"

  ![image](https://github.com/user-attachments/assets/1378d547-035d-4da0-8b4d-80dddd1939b6)

  Result

  You get .csv file with data. In each column for each enzyme the number of sites of this enzyme in each sequence and its positions in square brackets are written.

  ![image](https://github.com/user-attachments/assets/360ebb1b-66f9-454d-baac-db24ab3d38dc)

  You can filter count sites by command -fcs.

  PDRF_PCR_Analysis analyze "D:\06_Python\py\Cry_genes_seq" -fcs 30 5

  Result

  ![image](https://github.com/user-attachments/assets/3a2add63-19d0-4e2f-b385-cd3927de631e)


  



