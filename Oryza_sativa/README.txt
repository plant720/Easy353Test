1. The simulated data is simulated by art_illumina, and the command to simulate data of Oryza_sativa is:
	art_illumina -ss HS25 -i Osativa_323_v7.0.cds.fa -o Osativa_1x -l 150 -f 10 -p -m 200 -s 10 -sam
	art_illumina -ss HS25 -i Osativa_323_v7.0.cds.fa -o Osativa_5x -l 150 -f 10 -p -m 200 -s 10 -sam
	art_illumina -ss HS25 -i Osativa_323_v7.0.cds.fa -o Osativa_10x -l 150 -f 10 -p -m 200 -s 10 -sam
	art_illumina -ss HS25 -i Osativa_323_v7.0.cds.fa -o Osativa_20x -l 150 -f 10 -p -m 200 -s 10 -sam
	art_illumina -ss HS25 -i Osativa_323_v7.0.cds.fa -o Osativa_50x -l 150 -f 10 -p -m 200 -s 10 -sam

1. HybPiper version is 
	hybpiper 2.0.1rc build 11

2. The command to run Hybpiper is 
	hybpiper assemble -r data_1.fq data_2.fq -t_dna ../353_Poaceae_hybpiper.fasta  --prefix test --bwa --cpu 4

3. The command to run Easy353 is 
	easy353.py -r ../353_ref_Poaceae -1 data_1.fq -2 data_2.fq -o test -k1 29 -k2 41 -kmer_limit 8

4. INSDC.PRJDB1747.Oryza_sativa.a353.fasta  
	The Angiosperms353 sequences of Oryza_sativa from Kew Tree of Life Explorer (https://treeoflife.kew.org)

5. Osativa_323_v7.0.cds.fa.gz
	The nucleotide FASTA format data of all CDS of Oryza sativa from Phytozome (https://phytozome-next.jgi.doe.gov)

6. Oryza_sativa_gold.fasta
	The gold standard for comparison, which is the Angiosperms353 sequences corrected by CDS of Oryza sativa

7. 353_Poaceae.zip
	The Angiosperms353 sequences of the Poaceae family, which are used as reference sequences.

