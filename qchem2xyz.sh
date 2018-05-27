#!/bin/bash

pattern='/molecule/,/end/{/molecule/!{/end/!p}}'
dirname='chem_database'

for filename in $(ls $dirname)
do
    mol="${filename%.*}"
    input="$mol".in
    output="$mol".xyz
    grep -c "^ " "$dirname/$input" > "$dirname/$output" && sed -n $pattern "$dirname/$input" >> "$dirname/$output"
done
