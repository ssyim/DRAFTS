#!/usr/bin/env bash

# Code to BBmerge of paired-end reads

search_dir="$1"
out_dir="$2"

if [ -z "${search_dir}" ]; then # if no args, output usgae and exit
	echo "Usage:"
	echo "./01_DRAFTS_merge_pairedreads.sh [dir_to_search] [out_dir (optional)]"
	exit 1
fi

if [ -z "${out_dir}" ]; then # if no out dir, then save in the dir being searched
	out_dir="${search_dir}"
fi

echo "Keep pre-merged data (y/n)? " # User input
read keep_data

echo "Merging paired-end reads"

if [ "${keep_data}" == "y" ]; then
	find "${search_dir}" -mindepth 1 -type f \( -not -iname ".*" \) | rev | cut -f1 -d"/" | rev | sed 's/_R/#/g' | cut -f1 -d"#" | sort -u |\
		while read -r fn; do
			sudo bash */bbmap/bbmerge.sh in1=$search_dir/"$fn"_R2.fastq in2=$search_dir/"$fn"_R1.fastq out=$out_dir/"$fn"-merged.fastq outu1=$out_dir/unmerged-"$fn"-R2.fastq outu2=$out_dir/unmerged-"$fn"-R1.fastq
		done
else
	find "${search_dir}" -mindepth 1 -type f \( -not -iname ".*" \) | rev | cut -f1 -d"/" | rev | sed 's/_R/#/g' | cut -f1 -d"#" | sort -u |\
		while read -r fn; do
			sudo bash */bbmap/bbmerge.sh in1=$search_dir/"$fn"_R2.fastq in2=$search_dir/"$fn"_R1.fastq out=$out_dir/"$fn"-merged.fastq outu1=$out_dir/unmerged-"$fn"-R2.fastq outu2=$out_dir/unmerged-"$fn"-R1.fastq
			rm -f $search_dir/"$fn"_R1.fastq $search_dir/"$fn"_R2.fastq
		done
fi

echo "Done"