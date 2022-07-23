# 1. simulate NGS data of Oryza_sativa CDS with different sequencing depths
cat depth.list | while read a; do art_illumina -ss HS25 -i Osativa_323_v7.0.cds.fa -o Osativa_${a}_ -l 150 -f 10 -p -m 200 -s 10 -sam

# 2. run the HybPiper
cat depth.list | while read a; do hybpiper assemble -r Osativa_${a}_1.fq Osativa_${a}_2.fq -t_dna 353_Poaceae_hybpiper.fasta --prefix hyb_${a} --bwa --cpu 4;done

# 3. run the Easy353
cat depth.list | while read a; do easy353.py -r 353_ref_Poaceae -1 Osativa_${a}_1.fq -2 Osativa_${a}_2.fq -o easy_${a} -k1 29 -k2 41 -kmer_limit 8;done
