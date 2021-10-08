#Copy input data into the input folder
cp ../../decontamination/output/GCMP_decontaminated_1000_tree-rooted.qza ../input/physeq_rooted_tree.qza
cp ../../decontamination/output/GCMP_decontaminated_1000_feature-table.qza ../input/physeq_feature_table.qza
cp ../../decontamination/output/GCMP_decontaminated_1000_sample-metadata.txt ../input/physeq_metadata.txt
cp ../../decontamination/output/GCMP_decontaminated_1000_tax.txt ../input/physeq_taxonomy.txt

cp ../../coral_disease_data/make_disease_graphs/output/merged_disease_table.tsv ../input/merged_disease_table.tsv
