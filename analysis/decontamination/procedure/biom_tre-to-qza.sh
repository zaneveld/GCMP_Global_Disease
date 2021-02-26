qiime tools import \
  --input-path physeq.noncont_feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path physeq.noncont-feature-table.qza
  
qiime tools import \
  --input-path physeq.noncont_tree-rooted.newick \
  --output-path physeq.noncton-rooted-tree.qza \
  --type 'Phylogeny[Rooted]'