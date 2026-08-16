#!/bin/bash
# This is the sub-pipe of 3'aQTL-pipe, here in this script we will perform fine-mapping of 3'aQTL detected by Matrix-eQTL with SuSieR
# @Xudong Zou, zouxd@szbl.ac.cn
# 2022-03-30

# -- Usage function
script_name=$0
function usage(){
	echo "#=============================="
	echo "Default usage:"
	echo "#=============================="
	echo "bash $script_name"
	echo "Options:"
	echo "        -w  integer,setting the window size around the genes for fine-mapping"
	echo "        -p  float, specify the minimum PIP for filtering fine mapped 3'aQTLs"
	echo "        -L  integer, specify the L value in susieR, default 10"
	echo "        -V  float, specify the variance used in susieR, default 0.2"
	echo "        -t  integer, setting threads to run susieR in parallel, default 1"
	echo "        -h  print the help information"
	exit 1
}

# define global variables from command parameters
currDir=`pwd`
sourceDir="/share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src"
PIP="0.1"
Variance="0.2"
L="10"
Threads="1"
window=`echo "1e6"|awk '{printf("%d",$0)}'`
Rscript="/share/pub/xingsl/shilai/project/APA/software/miniforge3/envs/R3.6.3/bin/Rscript"

while getopts :w:p:L:V:t:h opt
do
	case $opt in
		w)
			window=`echo "$OPTARG"|awk '{printf("%d",$0)}'`
		;;
		p)
			PIP="$OPTARG"
		;;
		L)
			L="$OPTARG"
		;;
		V)
			Variance="$OPTARG"
		;;
		t)
			Threads="$OPTARG"
		;;
		h)
			echo "Help message:"
			usage
		;;
		:)
			echo "The option -$OPTARG requires an argument."
			exit 1
		;;
		?)
			echo "Invalid option: $OPTARG"
			usage
			exit 2
		;;
	esac
done

# -- Basic settings
if [ ! -d "${currDir}/FineMapping/input" ]
then
	echo "${currDir}/FineMapping/input not found!"
	exit
fi

if [ ! -d "${currDir}/FineMapping/output" ]
then
	echo "${currDir}/FineMapping/output not found!"
	exit
fi

if [ ! -d "${currDir}/Matrix_eQTL" ]
then
        echo "No Matrix-eQTL output found!"
        exit
fi

# -- BLAS 线程限制（必须前置，否则 PBS 会 burst） --
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

# -- Main function --
function main(){
	echo "Running $script_name with the following parameters:"
	echo "*************************************************"
	echo "-w: $window"
	echo "-p: $PIP"
	echo "-L: $L"
	echo "-V: $Variance"
	echo "-t: $Threads"
	echo "*************************************************"
	echo "Start 3'aQTL fine mapping ..."
	date
	echo "Run fine-mapping analysis by susieR"
	run_fine_mapping $window $PIP $L $Variance $Threads
	echo "Done!"
	date
}

# -- Other functions --
function run_fine_mapping(){
	w=$1
	min_PIP=$2
	L=$3
	Var=$4
	threads=$5
	if [ ! -f "${currDir}/FineMapping/input/picked_asso_list.loc_${w}.txt" ]
	then
		echo "File ${currDir}/FineMapping/input/picked_asso_list.loc_${w}.txt not found!"
		exit
	fi

	if [ $threads -eq 1 ]
	then
		# ---- 单线程：顺序跑 ----
		for gene in `cat ${currDir}/FineMapping/input/picked_asso_list.loc_${w}.txt | cut -f1`
		do
			if [ -d "${currDir}/FineMapping/output/$gene" ]
			then
				GeneDir=${currDir}/FineMapping/output/$gene
				# 已完成：3aQTL.SuSiE.txt 存在 → 跳过
				if [ -f "${GeneDir}/3aQTL.SuSiE.txt" ]
				then
					echo "${gene}: already done, skip"
					continue
				elif [ -f "${GeneDir}/3aQTL.vcf" -a -f "${GeneDir}/expr.phen" ]
				then
					echo "Analyzing $gene"
					$Rscript ${sourceDir}/finemapping.R ${GeneDir} $L $Var $min_PIP
				else
					echo "${gene}: File 3aQTL.vcf and expr.phen not found!"
					continue
				fi
			else
				echo "${gene} not exits!"
				continue
			fi
		done
		cd $currDir
	else
		# ---- 多线程：xargs -P $threads 精确控制并发 ----
		# 不再用 split + 嵌套 for，也不再需要 tmp/finemap_task_* 中间文件
		cut -f1 ${currDir}/FineMapping/input/picked_asso_list.loc_${w}.txt | \
		xargs -P $threads -I {} bash -c '
			gene="{}"
			GeneDir="'"${currDir}"'/FineMapping/output/${gene}"
			if [ -d "${GeneDir}" ]; then
				# 已完成：3aQTL.SuSiE.txt 存在 → 跳过
				if [ -f "${GeneDir}/3aQTL.SuSiE.txt" ]; then
					echo "${gene}: already done, skip"
				elif [ -f "${GeneDir}/3aQTL.vcf" ] && [ -f "${GeneDir}/expr.phen" ]; then
					echo "Analyzing ${gene}"
					'"$Rscript"' '"${sourceDir}"'/finemapping.R ${GeneDir} '"$L"' '"$Var"' '"$min_PIP"'
				else
					echo "${gene}: File 3aQTL.vcf and expr.phen not found!"
				fi
			else
				echo "${gene} not exits!"
			fi
		'
		cd $currDir
	fi

}

# - main
main

