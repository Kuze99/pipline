#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=10G
#$ -l mem_req=10G
#$ -pe def_slot 4

##
#このスクリプトは親スクリプトから呼び出すので
#conf の読み込みは親スクリプトだけにする(親スクリプトの変数を引き継ぐ)
#パラメータチェックの関数をパイプライン全体でくくりだし、親スクリプトでcheckさせる
#メモリとコア数も親スクリプトのqsub時に決定する。(qsubを入れ子にすることを想定)
#
#各step毎に異常終了を拾うようにする
#
#FASTQが圧縮されていた場合のSTAR option (→変数に入れ込んでif分で分岐するか、そもそものinputを統一させるか)

#------------------------------------------------------
# USAGE
#------------------------------------------------------

usage(){
        cat <<-EOF
    USAGE:
	qsub -cwd $0 -f Read1.fq (-r Read2.fq)  -n name -c config.txt
    Desc:
          STARのindexを作成するscript
          STAR dirを作成し、gtfファイル

    Options
        -f ... Read1
		-r ... Read2 (PE Only)
        -n ... Sample_name
		-c ... config.txt
        -h ... help
EOF
exit;
}

#------------------------------------------------------
# PARM
#------------------------------------------------------

#------------------------------------------------------
# GET OPT
#------------------------------------------------------

R1=Read1
NAME=SAMPLE
CONF=./conf.txt
# PARE or Single
MODE=PE
CORE=4

while getopts f:r:n::h OPT;do
        case $OPT in

                f )R1=$OPTARG # read REN
                        ;;
                r )R2=$OPTARG # read REN
                        ;;
                n )NAME=$OPTARG # read REN
                        ;;
				c )CONF=$OPTARG # read REN
						;;
                h )usage
        esac
done

#------------------------------------------------------
# CHECK PARM 後で関数かする
#------------------------------------------------------

if [ -z "$R2" ] ;then
	MODE=SE
elif [ ! -e $R2 ];then
	 echo -e "$R2 が存在しておりません";exit
fi

for i in $R1 $CONF ;do
        if [ ! -e $i ];then
                echo -e "$i が存在しておりません";exit
        fi
done

#------------------------------------------------------
#MAIN
#------------------------------------------------------
. ./conf.txt

STAR --genomeDir $STAR_INDEX \
	--runThreadN $CORE \
	--sjdbGTFfile $GTF \
	--outSAMtype BAM SortedByCoordinate\
	--outSAMstrandField intronMotif \
	--outFileNamePrefix  ${NAME}_ \
	--outFilterMultimapNmax 1 \
	--readFilesCommand zcat \
	--readFilesIn $R1 $R2

samtools index ${NAME}_Aligned.sortedByCoord.out.bam

# Feature Count Genes
if [ $MODE = "PE" ] ;then
	FCOUNT_OPT=" -p -B "
elif [ $MODE != "SE" ] ;then
	echo "ERROR feature count に不明なoptionsがわたりました";
	exit
fi

featureCounts ${FCOUNT_OPT} -a $GTF -o ${NAME}_unstarand ${NAME}_Aligned.sortedByCoord.out.bam -s 0 > COUNT_LOG.unstarand.txt 2>&1
featureCounts ${FCOUNT_OPT} -a $GTF -o ${NAME}_starand ${NAME}_Aligned.sortedByCoord.out.bam -s 1 > COUNT_LOG.starand.txt  2>&1
featureCounts ${FCOUNT_OPT} -a $GTF -o ${NAME}_reverse_starand ${NAME}_Aligned.sortedByCoord.out.bam -s 2 > COUNT_LOG.rev_unstarand.txt 2>&1


# MAKE STATS


# Calk RPKM


exit 0;
