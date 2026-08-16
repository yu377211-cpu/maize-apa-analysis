
Rscript ./bin/facet_heatmap.r -i same_subpop/sameSubpop.ratio.txt -o sameSubpop
Rscript ./bin/facet_heatmap.r -i same_tissue/sameTissue.ratio.txt -o sameTissue 

python3 bin/stat_plot_NUM.py same_subpop/sameSubpop.NUMsta.txt sameSubpop
python3 bin/stat_plot_NUM.py same_tissue/sameTissue.NUMsta.txt sameTissue
