#!/bin/bash
echo "current wd"
pwd
echo "script start"
input=$1
echo "input" $1
output=$2
echo "output" $2
v=$3
dir_name=${input%/*}
file_name=${input##*/}
output_name=${output}/${file_name}copra.log
echo "files" $dir_name $file_name $output_name
script=chaldea/copra.jar
echo "script" $script
# preprocess
sed -i -e '1,3s/^/#/' $input
java -cp  $script COPRA $input  -w -mo -repeat 3 -v $v > $output_name
sed -i -e '/Invalid/d' $output_name
sed -i -e '/Ignoring/d' $output_name
cat $output_name
echo output $output/best-clusters-$file_name
mv best-clusters-$file_name $output/best-clusters-$file_name
# post process
sed -i -e '1,3s/^#//' $input


