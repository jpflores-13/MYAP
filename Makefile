.PHONY: clean

objects :=\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.h5\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_sorb_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_cont_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_omega_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_loopCounts_10kb.h5\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_sorb_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_cont_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_omega_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_omega_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_cont_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_sorb_loops_10kb.h5\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.h5\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.rds\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_sorb_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_cont_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_omega_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_loopCounts_10kb.rds\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_sorb_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_cont_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/noDroso/noDroso_omega_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_omega_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_cont_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_sorb_loops_10kb.rds\
data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.rds\
data/processed/hic/hg38/diffLoops/bothDroso/diffLoops_bothDroso_10kb.rds\
data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds\
data/processed/hic/hg38/diffLoops/isDroso/diffLoops_isDroso_10kb.rds\

######################## figure out way to compact this more efficiently
data/processed/hic/dm6/loop_decay/dmLoops_1kb_max100_mh_index.rds\
data/processed/hic/hg38/loop_decay/bothDroso/gainedLoops_bothDroso_1kb_max50_mh_index.rds\
data/processed/hic/hg38/loop_decay/bothDroso/nullSet_bothDroso_1kb_max50_mh_index.rds\
data/processed/hic/hg38/loop_decay/bothDroso/ctcfLoops_bothDroso_1kb_max50_mh_index.rds\

data/processed/hic/hg38/loop_decay/bothDroso/gainedLoops_bothDroso_1kb_max100_mh_index.rds\
data/processed/hic/hg38/loop_decay/bothDroso/nullSet_bothDroso_1kb_max100_mh_index.rds\
data/processed/hic/hg38/loop_decay/bothDroso/ctcfLoops_bothDroso_1kb_max100_mh_index.rds\


data/processed/hic/hg38/loop_decay/noDroso/gainedLoops_noDroso_1kb_max50_mh_index.rds\
data/processed/hic/hg38/loop_decay/noDroso/nullSet_noDroso_1kb_max50_mh_index.rds\
data/processed/hic/hg38/loop_decay/noDroso/ctcfLoops_noDroso_1kb_max50_mh_index.rds\

data/processed/hic/hg38/loop_decay/noDroso/gainedLoops_noDroso_1kb_max100_mh_index.rds\
data/processed/hic/hg38/loop_decay/noDroso/nullSet_noDroso_1kb_max100_mh_index.rds\
data/processed/hic/hg38/loop_decay/noDroso/ctcfLoops_noDroso_1kb_max100_mh_index.rds\


plots/loop_decay_bothDroso_1kb_max50.pdf\
plots/loop_decay_bothDroso_1kb_max100.pdf\

plots/loop_decay_noDroso_1kb_max50.pdf\
plots/loop_decay_noDroso_1kb_max100.pdf\
###################

plots/tile_chr1.pdf\
plots/diffLoops_bothDroso_10kb_MA.pdf\
plots/diffLoops_noDroso_10kb_MA.pdf\
plots/diffLoops_isDroso_10kb_MA.pdf\
plots/drosoAnalysis_10kb.pdf\
plots/drosoAnalysis_5kb.pdf\
plots/EGFP-YAP_ChIP_surveyPlot_10kb.pdf\
plots/genotypeComparison_apa.pdf\
plots/motifEnrichment_ame.pdf

all: $(objects)
	echo done!
	
clean:
	rm -rf $(objects)

data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.h5:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_sorb_loops_10kb.h5:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_cont_loops_10kb.h5:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_omega_loops_10kb.h5:\
		scripts/processing/extract_bothDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
			Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/noDroso/noDroso_loopCounts_10kb.h5:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/noDroso
			Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/noDroso/noDroso_sorb_loops_10kb.h5:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/noDroso
			Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/noDroso/noDroso_cont_loops_10kb.h5:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/noDroso
			Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/noDroso/noDroso_omega_loops_10kb.h5:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/noDroso
			Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/isDroso/isDroso_omega_loops_10kb.h5:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/isDroso
			Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/isDroso/isDroso_cont_loops_10kb.h5:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/isDroso
			Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/isDroso/isDroso_sorb_loops_10kb.h5:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/isDroso
			Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.h5:	
	scripts/processing/extract_isDroso_loopCounts_10kb.R
			mkdir -p data/processed/hic/hg38/h5/loops/isDroso
			Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R
			
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.rds:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R
		
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_sorb_loops_10kb.rds:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R
		
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_cont_loops_10kb.rds:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R
		
data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_omega_loops_10kb.rds:\
	scripts/processing/extract_bothDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/bothDroso
		Rscript scripts/processing/extract_bothDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/noDroso/noDroso_loopCounts_10kb.rds:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/noDroso
		Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/noDroso/noDroso_sorb_loops_10kb.rds:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/noDroso
		Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/noDroso/noDroso_cont_loops_10kb.rds:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/noDroso
		Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R
		
data/processed/hic/hg38/h5/loops/noDroso/noDroso_omega_loops_10kb.rds:\
	scripts/processing/extract_noDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/noDroso
		Rscript scripts/processing/extract_noDroso_loopCounts_10kb.R
		
data/processed/hic/hg38/h5/loops/isDroso/isDroso_omega_loops_10kb.rds:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/isDroso
		Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/isDroso/isDroso_cont_loops_10kb.rds:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/isDroso
		Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/isDroso/isDroso_sorb_loops_10kb.rds:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/isDroso
		Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R

data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.rds:\
	scripts/processing/extract_isDroso_loopCounts_10kb.R\
	scripts/utils/saveLoopCounts.R
		mkdir -p data/processed/hic/hg38/h5/loops/isDroso
		Rscript scripts/processing/extract_isDroso_loopCounts_10kb.R

data/processed/hic/hg38/diffLoops/bothDroso/diffLoops_bothDroso_10kb.rds:\
	data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.rds\
	scripts/processing/call_diffLoops_bothDroso_10kb.R
		mkdir -p data/processed/hic/hg38/diffLoops/bothDroso
		Rscript scripts/processing/call_diffLoops_bothDroso_10kb.R
	
data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds:\
	data/processed/hic/hg38/h5/loops/noDroso/noDroso_loopCounts_10kb.rds\
	scripts/processing/call_diffLoops_noDroso_10kb.R
		mkdir -p data/processed/hic/hg38/diffLoops/noDroso
		Rscript scripts/processing/call_diffLoops_noDroso_10kb.R
		
data/processed/hic/hg38/diffLoops/isDroso/diffLoops_isDroso_10kb.rds:\
	data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.rds\
	scripts/processing/call_diffLoops_isDroso_10kb.R
		mkdir -p data/processed/hic/hg38/diffLoops/isDroso
		Rscript scripts/processing/call_diffLoops_isDroso_10kb.R

plots/tile_chr1.pdf:\
	scripts/analysis/tile_chr1.R
		mkdir -p plots
		Rscript scripts/analysis/tile_chr1.R
		
plots/diffLoops_bothDroso_10kb_MA.pdf:\
	data/processed/hic/hg38/h5/loops/bothDroso/bothDroso_loopCounts_10kb.rds\
	scripts/processing/call_diffLoops_bothDroso_10kb.R
		mkdir -p plots
		Rscript scripts/processing/call_diffLoops_bothDroso_10kb.R

plots/diffLoops_noDroso_10kb_MA.pdf:\
	data/processed/hic/hg38/h5/loops/noDroso/noDroso_loopCounts_10kb.rds\
	scripts/processing/call_diffLoops_noDroso_10kb.R
		mkdir -p plots
		Rscript scripts/processing/call_diffLoops_noDroso_10kb.R
		
plots/diffLoops_isDroso_10kb_MA.pdf:\
	data/processed/hic/hg38/h5/loops/isDroso/isDroso_loopCounts_10kb.rds\
	scripts/processing/call_diffLoops_isDroso_10kb.R
		mkdir -p plots
		Rscript scripts/processing/call_diffLoops_isDroso_10kb.R
		
plots/drosoAnalysis_10kb.pdf:\
	scripts/analysis/drosoAnalysis_10kb.R\
	data/processed/hic/hg38/diffLoops/*
		mkdir -p plots
		Rscript scripts/analysis/drosoAnalysis_10kb.R
		
plots/drosoAnalysis_5kb.pdf:\
	scripts/analysis/drosoAnalysis_5kb.R\
	data/processed/hic/hg38/diffLoops/*
		mkdir -p plots
		Rscript scripts/analysis/drosoAnalysis_5kb.R
		
plots/EGFP-YAP_ChIP_surveyPlot_10kb.pdf:\
	scripts/analysis/EGFP-YAP_ChIP_surveyPlot_10kb.R\
	data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds
		mkdir -p plots
		Rscript scripts/analysis/EGFP-YAP_ChIP_surveyPlot_10kb.R
		
plots/genotypeComparison_apa.pdf:\
	scripts/analysis/genotypeComparison_apa.R\
	data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds
		mkdir -p plots
		Rscript scripts/analysis/genotypeComparison_apa.R
		
plots/motifEnrichment_ame.pdf:\
	scripts/analysis/motifEnrichment_ame.R\
	data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds\
	data/processed/atac/hg38/diff_ATACcounts.rds
		mkdir -p plots
		Rscript scripts/analysis/motifEnrichment_ame.R
	