#!/bin/bash
SRA=$1
DIR=$2 
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l 200m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instan    t/reads/ByRun/sra/SRR/${SRA:4:6}/${SRA}/${SRA}.sra ${DIR}
