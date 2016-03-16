all: homework

clean: 
	rm -rf SamHinshawHomework.md SamHinshawHomework.html data_fixed.txt design_fixed.txt yeast_counts_fixed.tsv

data: 
	curl -o data.txt.gz 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_transcriptome_data.txt.gz?raw=true'
	gunzip -kq data.txt.gz

design: 
	curl -o design.txt 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_design.txt?raw=true'
	
yeast: 
	curl -o yeast.tsv 'http://stat540-ubc.github.io/homework/assignment/homework_data/GSE37599-data.tsv?raw=true'
	curl -o yeast_counts.tsv 'http://stat540-ubc.github.io/homework/assignment/homework_data/stampy.deep.counts.tsv?raw=true'
	curl -o yeast_counts_low.tsv 'http://stat540-ubc.github.io/homework/assignment/homework_data/stampy.low.counts.tsv?raw=true'

gitignore: data.txt.gz design.txt data.txt yeast.tsv
	grep -q -F 'data.txt' ../.gitignore || echo 'data.txt' >> ../.gitignore
	grep -q -F 'data.txt.gz' ../.gitignore || echo 'data.txt.gz' >> ../.gitignore
	grep -q -F 'design.txt' ../.gitignore || echo 'design.txt' >> ../.gitignore
	grep -q -F 'design_fixed.txt' ../.gitignore || echo 'design_fixed.txt' >> ../.gitignore
	grep -q -F 'data_fixed.txt' ../.gitignore || echo 'data_fixed.txt' >> ../.gitignore
	grep -q -F 'yeast.tsv' ../.gitignore || echo 'yeast.tsv' >> ../.gitignore

data_fixed: data.txt
	cp data.txt data_fixed.txt
	sed -i -e 's/^"GSE10718_Biomat_1"/"ProbeID"	"GSE10718_Biomat_1"/' data_fixed.txt

design_fixed: design.txt
	cp design.txt design_fixed.txt
	sed -i -e 's/^"ExternalID"/"InternalID"	"ExternalID"/' design_fixed.txt

counts_fixed: yeast_counts.tsv
	cp yeast_counts.tsv yeast_counts_fixed.tsv
	sed -i -e 's/^""/"gene"/' yeast_counts_fixed.tsv

low_counts_fixed: yeast_counts_low.tsv
	cp yeast_counts_low.tsv yeast_counts_low_fixed.tsv
	sed -i -e 's/^""/"gene"/' yeast_counts_low_fixed.tsv

homework:  SamHinshawHomework.Rmd gitignore data_fixed design_fixed counts_fixed
	Rscript -e "rmarkdown::render('SamHinshawHomework.Rmd')"
	