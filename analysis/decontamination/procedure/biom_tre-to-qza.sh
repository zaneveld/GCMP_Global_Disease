qiime tools import \
  --input-path GCMP_decontaminated_feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path GCMP_decontaminated_feature-table.qza
  
qiime tools import \
  --input-path GCMP_decontaminated_1000_feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path GCMP_decontaminated_1000_feature-table.qza  

qiime tools import \
  --input-path GCMP_decontaminated_tree-rooted.newick \
  --output-path GCMP_decontaminated_tree-rooted.qza \
  --type 'Phylogeny[Rooted]'

qiime tools import \
  --input-path GCMP_decontaminated_1000_tree-rooted.newick \
  --output-path GCMP_decontaminated_1000_tree-rooted.qza \
  --type 'Phylogeny[Rooted]'

qiime tools import \
    --input-path GCMP_decontaminated_tax.txt \
    --output-path GCMP_decontaminated_tax.qza \
    --input-format HeaderlessTSVTaxonomyFormat \
    --type 'FeatureData[Taxonomy]'

qiime tools import \
    --input-path GCMP_decontaminated_1000_tax.txt \
    --output-path GCMP_decontaminated_1000_tax.qza \
    --input-format HeaderlessTSVTaxonomyFormat \
    --type 'FeatureData[Taxonomy]'

