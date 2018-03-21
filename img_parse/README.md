# Scripts for working with IMG annotation output

The JGI IMG annotation pipeline results can be downloaded as tab-separated annotation tables and a zipped "download bundle".

These scripts are to help in format conversions for downstream analyses and sequence submissions. IMG does not submit sequences to the INSDC, so users will have to perform their own sequence submissions.

To see help message for any script, run `perl scriptname.pl` without arguments or `perl scriptname.pl --help`.

 * `img_tsv2gbk.pl` -- Convert tab-separated annotation table to Genbank format
 * `tidy_img_gff3.pl` -- Tidy up GFF and Fasta nucleotide (.fna) files from IMG to pass to the [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) utility.
 * `img_xls_extract_KO.pl` -- Extract KO numbers from IMG annotation table, for KEGG module reconstruction
 * `kegg_module_reconstructions_tabulate.pl` -- Tabulate KEGG module reconstructions from the online [Reconstruct Module tool](http://www.kegg.jp/kegg/tool/map_module.html).
