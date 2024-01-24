awk -F'\t' '{print $1,$2}' OFS='\t' $1 | sed 's/\.\./\t/g' | awk -F'\t' '!seen[$2$3]++' OFS='\t' | sed '/primer/Id'
