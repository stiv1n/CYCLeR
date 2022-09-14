# Working with CYCLeR

**CYCLeR** is a pipeline for reconstruction of circRNA transcripts from RNA-seq data and their subsequent quantification. The algorithm relies on comparison between control total RNA-seq samples and circRNA enriched samples to identify circRNA specific features. Then the selected circRNA features are used to infer the transcripts through a graph-based algorithm. Once the predicted transcript set is assembled, the transcript abundances are estimated through an EM algorithm with **kallisto**. **CYCLeR** takes as an input BAM files and back-splice junction (BSJ) files and outputs transcript infomation in different formats, including a FASTA output for abundance estimation. 

## Installation of CYCLeR

### Command line tools needed
The computation steps prior and post **CYCLeR** run are most efficiently run on HPC. It is very likely that any HPC in biological institute already has most of those tools installed. Just in case, a **Docker** image containing all the tools is provided.  
NOTE: prior to running **Docker** image, make sure that ***Docker** is indeed installed and working: https://docs.docker.com/get-started/

* **STAR** - https://github.com/alexdobin/STAR
* **samtools** - https://sourceforge.net/projects/samtools/files/samtools/
* **kallisto** - http://pachterlab.github.io/kallisto/download
* **bwa** (needed for CIRI2) - http://bio-bwa.sourceforge.net/bwa.shtml
* **CIRI2** - https://sourceforge.net/projects/ciri/files/CIRI2/
* **CIRCexplorer2** - https://circexplorer2.readthedocs.io/en/latest/
```
#Docker image with all command line tools  
sudo docker pull stiv1n/cycler.prerequisites
```
### R test run installation
For a test run, I suggest using a **Docker** container. There, all test input and all dependencies are provided. 
The **Docker** use requires you to mount a volume - a working directory (*<local_dir>*) where the output and input would be stored.
This container uses **RStudio server** and required login. In this case, the  username is *rstudio* the password is *guest*.
```
sudo docker pull stiv1n/cycler
sudo docker run --rm -ti -e PASSWORD=guest -v <local_dir>:/usr/workdir -p 8787:8787 stiv1n/cycler
```
## Full documentation
For more information on installation, pre-processing and core tool run please see the [vignette](https://raw.githubusercontent.com/stiv1n/CYCLeR/main/CYCLeR_workflow.pdf) and the [manual](https://raw.githubusercontent.com/stiv1n/CYCLeR/main/CYCLeR.pdf). 

## CYCLeR
To test the tool, use the R [script](https://raw.githubusercontent.com/stiv1n/CYCLeR/main/docker_test.R) provided and just run it in the RStudio server.

## Quantification 
The final step of the **CYCLeR** pipeline is running **kallisto** with the newly assembled transcriptome.    


```
sudo docker run  -v <local_dir>:/usr/local stiv1n/cycler.prerequisites \
  kallisto index -i kallisto_index -k 31 for_kallisto.fa
sudo docker run  -v <local_dir>:/usr/local stiv1n/cycler.prerequisites \
  kallisto quant -i kallisto_index -o ./ <sample_name>_1.fastq <sample_name>_2.fastq
```


