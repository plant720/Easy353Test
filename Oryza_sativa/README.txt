1. The simulated data is simulated by art_illumina (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

2. INSDC.PRJDB1747.Oryza_sativa.a353.fasta  
	The Angiosperms353 sequences of Oryza_sativa from Kew Tree of Life Explorer (https://treeoflife.kew.org)

3. Osativa_323_v7.0.cds.fa.gz
	The nucleotide FASTA format data of all CDS of Oryza sativa from Phytozome (https://phytozome-next.jgi.doe.gov)

4. Oryza_sativa_gold.fasta
	The gold standard for comparison, which is the Angiosperms353 sequences corrected by CDS of Oryza sativa

5. 353_Poaceae
	The directory to store Angiosperms353 sequences of the Poaceae family, are used as reference sequences for Easy353

6. 353_Poaceae_hybpiper.fasta
	The reference sequences used for HybPiper, are modified and merged from 353_Poaceae

7. depth.list
	The file records the sequencing depth
	
8. work.sh
	The file records the command to run the test (Note: you should install ART, HybPiper 2.0, Easy353)

9. compare.py
	A script for calculating the identity and coverage between the recovered sequences and the gold standard sequences.
	
	
Reference
	1. Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a next-generation sequencing read simulator, Bioinformatics (2012) 28 (4): 593-594
	2. Johnson MG, Gardner EM, Liu Y, Medina R, Goffinet B, Shaw AJ, Zerega NJ, Wickett NJ. HybPiper: Extracting coding sequence and introns for phylogenetics from high-throughput sequencing reads using target enrichment. Appl Plant Sci. 2016 Jul 12;4(7):apps.1600016. doi: 10.3732/apps.1600016. PMID: 27437175; PMCID: PMC4948903.
	3. Baker W.J., Bailey P., Barber V., Barker A., Bellot S., Bishop D., Botigue L.R., Brewer G., Carruthers T., Clarkson J.J., Cook J., Cowan R.S., Dodsworth S., Epitawalage N., Francoso E., Gallego B., Johnson M., Kim J.T., Leempoel K., Maurin O., McGinnie C., Pokorny L., Roy S., Stone M., Toledo E., Wickett N.J., Zuntini A.R., Eiserhardt W.L., Kersey P.J., Leitch I.J. & Forest F. 2022. A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life. Systematic Biology 71: 301–319. https://doi.org/10.1093/sysbio/syab035
	4. David M. Goodstein, Shengqiang Shu, Russell Howson, Rochak Neupane, Richard D. Hayes, Joni Fazo, Therese Mitros, William Dirks, Uffe Hellsten, Nicholas Putnam, and Daniel S. Rokhsar, Phytozome: a comparative platform for green plant genomics, Nucleic Acids Res. 2012 40 (D1): D1178-D1186



