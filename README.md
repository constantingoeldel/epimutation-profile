# Epimutation profile

## Extract CGs within the gbM genes

To extract the subset of CGs that are interesting for our analysis, run `bash subset.sh`

If there are any errors, verify that the symlinks point to the correct original files (the directory `MA3_new_total_original_methylome` and the annotation file `gbM_gene_anotation_extract_Arabidopsis.bed`), otherwise create them.

On the server, it takes around 5 minutes for the script to run.

To check the files created, run `head -x /mnt/extStorage/constantin/epimutation-profile/y/z` where `x` is the number of lines you want to view, `y` the state of analysis (`MA3_new_total_original_methylome` for analysis, `cg_bed` for all CGs of the methylome and `within_cs` the CGs within gbM genes) and `z` is the specific file your interested in. 

## View files in VS Code

A nice way to view the heads of files in VS Code: `head -400 gbM_gene_anotation_extract_Arabidopsis.bed | code -`