all: homework

clean: 
	rm -rf SamHinshawHomework.md SamHinshawHomework.html
	sed -i '/data.txt.gz/d' ../.gitignore
	sed -i '/design.txt/d' ../.gitignore

data: 
	curl -o data.txt.gz 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_transcriptome_data.txt.gz?raw=true'
	gunzip -k data.txt.gz

design: 
	curl -o design.txt 'http://stat540-ubc.github.io/homework/assignment/homework_data/NHBE_design.txt?raw=true'

gitignore:
	grep -q -F 'data.txt.gz' ../.gitignore || echo 'data.txt.gz' >> ../.gitignore
	grep -q -F 'design.txt' ../.gitignore || echo 'design.txt' >> ../.gitignore

homework:  SamHinshawHomework.Rmd data design gitignore
	Rscript -e "rmarkdown::render('SamHinshawHomework.Rmd')"
