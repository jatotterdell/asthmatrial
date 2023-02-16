# Definitions ----

OUT = ./out
SIM = ./sims
SIM_FILES := $(wildcard $(SIM)/sims*.R)

# Generic Targets ----

sims: $(OUT)/sims01_01.qs \ 
	$(OUT)/sims02_01.qs \ 
	$(OUT)/sims03_01.qs

# Run single simulations ----

$(OUT)/sims01_01.qs: sims/sims01.R
	mkdir -p $(OUT)
	Rscript $< -n 500 -c 22

$(OUT)/sims02_01.qs: sims/sims02.R
	mkdir -p $(OUT)
	Rscript $< -n 500 -c 22

$(OUT)/sims03_01.qs: sims/sims03.R
	mkdir -p $(OUT)
	Rscript $< -n 500 -c 22

# Notebooks ----

notebooks/asthmatrial-summarise-sims.pdf: notebooks/summarise-sims.qmd
	quarto render $< --to pdf

notebooks/asthmatrial-summarise-sims.docx: notebooks/summarise-sims.qmd
	quarto render $< --to docx