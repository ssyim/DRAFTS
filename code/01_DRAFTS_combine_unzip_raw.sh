#!/usr/bin/env bash

# Code to 1) find and combine raw nextseq data in search_dir, which are generally separated in four folders of each lane,
# 2) unzip them to the out_dir

search_dir="$1"
out_dir="$2"

if [ -z "${search_dir}" ]; then # if no args, output usgae and exit
	echo "Usage:"
	echo "./01_DRAFTS_find_unzip_raw.sh [dir_to_search] [out_dir (optional)]"
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

echo "Done"