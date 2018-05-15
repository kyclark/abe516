sample_table:
	awk -F"\t" 'BEGIN{OFS=","} NR > 1 {print $$12,$$14}' SraRunTable.txt | sed "s/_floral_development//" > sample_table.csv

genes:
	cut -f 1 -d "," gene_list.csv  | awk 'NR>1' | sed "s/\"//g"
