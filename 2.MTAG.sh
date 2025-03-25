# Shell script for meta-analysis with MTAG on linux
# MTAG was installed using conda, and requires at least 12G to run

conda activate mtag

python mtag.py --sumstats dyslexia_23andme_zscores_uncorrected.txt,word_reading_EUR_genlang_GCOFF_zscores.txt \
			--out ./mtag_reading_GCOFF_EUR \
			--n_min 0.0 \
			--stream_stdout \ 
			--fdr
			
# End