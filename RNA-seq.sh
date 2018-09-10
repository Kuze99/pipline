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
	qsub -cwd $0 -f Read1.fq (-r Read2.fq)  -i STARindex -g GTF 
    Desc:
          STARのindexを作成するscript
          STAR dirを作成し、gtfファイル

    Options
        -f ... input Fastq  (default = Fasta/genome.fa)
        -g ... input Gff  (default = Gene_Ref/genes.gtf)
        -w ... outputディレクトリ (default = STAR)
        -l ... mappingするRead長(100でよい)　(default = 100)
EOF
exit;
}

#------------------------------------------------------
# PARM
#------------------------------------------------------

#------------------------------------------------------
# GET OPT
#------------------------------------------------------

. ./conf.txt

R1=
R2=
# PARE or Single
MODE=PE

INDEX=Fasta/genome.fa
GTF=Gene_Ref/genes.gtf
WORK=./STAR
NAME=

while getopts f:g:l:w::h OPT;do
        case $OPT in

                f )FA=$OPTARG # read REN
                        ;;
                g )GTF=$OPTARG # read REN
                        ;;
                w )WORK=$OPTARG # read REN
                        ;;
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


if [ ! -e $WORK ]; then
        mkdir -p $WORK
else
        printf "\tINDEX \"%s\" already exists.\n\tplease check\n" $WORK;
        exit 1;
fi

#------------------------------------------------------
#MAIN
#------------------------------------------------------
cd $WORK

STAR    --genomeDir $INDEX \
	--sjdbGTFfile $GTF \
	--outSAMtype BAM SortedByCoordinate\
	--outSAMstrandField intronMotif \
	--outFileNamePrefix  $NAME \
	--outFilterMultimapNmax 1 \
	--readFilesIn $R1 $R2

samtools index Aligned.sortedByCoord.out.bam


# Feature Count Genes
if [ $MODE = "SE" ] ;then
	

else if [ $MODE = "PE" ] ;then

else 
	echo "ERROR feature count に不明なoptionsがわたりました";
	exit
fi

# Calk RPKM


