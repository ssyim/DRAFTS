#!/usr/bin/env bash

# Code to 1) find and combine raw nextseq data in search_dir, which are generally separated in four folders of each lane,
# 2) gunzip them to the out_dir
# 3) then BBmerge of paired-end reads

search_dir="$1"
out_dir="$2"

if [ -z "${search_dir}" ]; then # if no args, output usgae and exit
	echo "Usage:"
	echo "./01_DRAFTS_find_unzip_merge_raw.sh [dir_to_search] [out_dir (optional)]"
	exit 1
fi

if [ -z "${out_dir}" ]; then # if no out dir, then save in the dir being searched
	out_dir="${search_dir}"
fi

echo "Combining files of same samples into ${out_dir}"
find "${search_dir}" -mindepth 1 -type d \( -not -iname ".*" \) | \
while read -r dn; do

	fn=$(echo $dn | rev | cut -f1 -d"/" | rev)
	fn=${fn/_L00/#}
	fn=$(echo "$fn" | cut -f1 -d"#" )
	# echo "$fn"

    cat "$dn"/*R1*.fastq.gz >> $out_dir/"$fn"_R1.fastq.gz
    cat "$dn"/*R2*.fastq.gz >> $out_dir/"$fn"_R2.fastq.gz
done

echo "Unzipping files and removing compressed data"
find "${out_dir}" -mindepth 1 -type f \( -not -iname ".*" \) | \
while read -r fn; do
	gzip -d "${fn}"
done

echo "Merging paired-end reads"
find "${out_dir}" -mindepth 1 -type f \( -not -iname ".*" \) | rev | cut -f1 -d"/" | rev | sed 's/_R/#/g' | cut -f1 -d"#" | sort -u |\
while read -r fn; do
	sudo bash */bbmap/bbmerge.sh in1=$out_dir/"$fn"_R2.fastq in2=$out_dir/"$fn"_R1.fastq out=$out_dir/"$fn"-merged.fastq outu1=$out_dir/unmerged-"$fn"-R2.fastq outu2=$out_dir/unmerged-"$fn"-R1.fastq
	rm -f $out_dir/"$fn"_R1.fastq $out_dir/"$fn"_R2.fastq
done

echo "Done"