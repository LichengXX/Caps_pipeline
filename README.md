# Caps_pipeline
This is a simple pipeline to batch design Caps Marker via VCF.\
`git clone git@github.com:LichengXX/Caps_pipeline.git`
## Requirements:
  python3.X.X(conda install conda-forge::python)\
	cyvcf2(conda install bioconda::cyvcf2)\
	Biopython(conda install conda-forge::biopython)\
 or You can install it using the built-in script:\
 `cd caps-pipeline`\
 `bash ./install.sh`\
 It is recommended to use conda to install requirements.

## Quick start
Only vcf file, reference genome and restriction endonuclease file are needed:\
`python3 01_caps.py -v test.vcf -r reference_genome -e Enzyme_file -o output.tsv --filter`\
If do not want to filter vcf, you can add --no-filter parameter:\
`python3 01_caps.py -v test.vcf -r reference_genome -e Enzyme_file -o output.tsv --no-filter`\
The vcf files can be generated by variant detection software such as GATK and bcftools, and need to be in the standard vcf format.\
The reference genome should make index using `samtools faidx`.\
The restriction endonuclease file should have 3 colunms separated by tab, e.g.:\
AatII	GACGTC	GAC'GTC\
Acc65I	GGTACC	GGTAC'C\
The first column is the enzyme name, the second column is the enzyme sequence, and the third column is the enzyme site, with ' as the enzyme site mark.\
If you want to use other enzyme digestion files, you can make them yourself, as long as the format meets the above requirements.\
If the marker file is too large, you can choose to split the result file by chromosome:\
`bash ./02_split_tsv.sh output.tsv`
## Other
If you have any questions, please contact me at lichengeg999@gmail.com
