for f in $(grep -rIPl "\t" ../../Source)
do
    echo "Converting tabs to spaces in $f"
    expand -t 8 $f > tmp_file
    mv tmp_file $f
done
