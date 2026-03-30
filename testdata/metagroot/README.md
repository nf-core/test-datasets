# MetagRoot DB
- HMM profile files from metagroot can be downloaded using [this link](https://pavlopoulos-lab.org/metagroot/DownloadHmm) from the [Downloads page](https://pavlopoulos-lab.org/metagroot/Downloads)
- downloads a .tar.gz (2GB) file which when unzipped is a folder "Skylign" of many .hmm files 

## Generating the test dataset: HMM.tar.gz 

1. download all fasta files from [samplesheet.csv](https://github.com/nf-core/test-datasets/blob/proteinannotator/samplesheet/samplesheet.csv) 

2. Unzipped HMM.tar.gz to get Skylign folder of hmm files



3. Concatenated all hmm files into one 
    
    e.g. 
    ```bash
    printf '%s\0' skylign/*.hmm | xargs -0 cat -- > metagrootdb.hmm
    ```

4. construct binary compressed datafiles 

    e.g. 
    ```bash
    hmmpress metagrootdb.hmm
    ```

5. scan 

    e.g.
    ```bash
    hmmscan --tblout result_l_arg.txt metagrootdb.hmm l_arginase.faa
    hmmscan --tblout result_T1024.txt metagrootdb.hmm T1024.fasta 
    hmmscan --tblout result_T1026.txt metagrootdb.hmm T1026.fasta 
    ```

6. kept subset of hit files for test dataset

    e.g. 
    ``` bash
    cp ~/original_folder/skylign/F422947.hmm ~/a_new_folder/skylign/
    cp ~/original_folder/skylign/F465379.hmm ~/a_new_folder/skylign/
    cp ~/original_folder/skylign/F451321.hmm ~/a_new_folder/skylign/
    ...
    ```

7. zipped to resemble original again

    e.g. 
    ```
    tar -czvf HMM.tar.gz Skylign
    ```
    and that is the HMM.tar.gz subset file uploaded here

## The other: subset_metagroot.hmm
Are the files from step 6 above concatenated into one hmm file.  But the HMM.tar.gz most resembles the metagroot database if retrieved from the actual [Downloads page](https://pavlopoulos-lab.org/metagroot/Downloads).