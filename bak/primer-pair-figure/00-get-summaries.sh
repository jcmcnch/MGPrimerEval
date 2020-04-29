mkdir -p tsv-summaries

for item in 515Y 926R; do

    cp ../../*/output-classify-workflow/*$item*0-mismatch.summary.tsv tsv-summaries

done
