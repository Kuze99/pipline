#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=50G
#$ -l mem_req=50G

#------------------------------------------------------
# USAGE
#------------------------------------------------------
usage(){
	cat <<-EOF
    USAGE:

    Desc: 
          make STAR index
          ..../ver　の下で Fasta部分にFastaファイル(genoem.fa)を設置後
	  STAR dirを作成し、gtfファイル

    Options
	-l
　　　　
EOF
}

#------------------------------------------------------
# GET OPT
#------------------------------------------------------

LEN=100
FA=Fasta/genome.fa
GTF=Gene_Ref/genes.gtf

while getopts l::h OPT;do
	case $OPT in
		l )LEN=$OPTARG # read REN
			;;
		h )usage
	esac
done

#------------------------------------------------------
# CHECK PARM
#------------------------------------------------------

for i in $FA $GTF ;do
	if [ ! -e $i ];then
		echo -e "$i が存在しておりません";exit
	fi
done


exit
#------------------------------------------------------
#MAIN
#------------------------------------------------------
