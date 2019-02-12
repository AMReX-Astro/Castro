#/bin/sh
# standard deviation is via sum of squares expression
grep -i "Coarse TimeStep" $1 | tail -5 | awk '{sum += $6; sumsq += $6^2; count +=1} END {print sum/count " " sqrt(sumsq/count - (sum/count)^2)}'
