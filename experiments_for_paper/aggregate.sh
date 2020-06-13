#!/bin/bash
var=( $(ls *.csv) )


cat ${var[0]} | head -n1 > merged.csv
for f in *.csv; do cat "$f" | tail -n +2 >> merged.csv; done
