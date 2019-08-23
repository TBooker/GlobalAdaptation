## A li'l pipeline to call SNPs from trees and to use VCFtools to do a widowed analysis on the data

##	for i in $(ls */)
##		do

w_size=10000
w_step=10000
s_size=50

python3 /home/booker/work/FishersWaveFst/slim/bin/recapSprinkler.py \
	--input $1  \
	--output $1.vcf \
	--sample $s_size \
	--just_vcf
#exit 0
#		rep_10.generation_120000.shat_0.001.M_5.pA_0.txt

vcftools --vcf $1.vcf \
		--weir-fst-pop ../vcftoolsConfigs/p1.$s_size.individuals.txt \
		--weir-fst-pop ../vcftoolsConfigs/p2.$s_size.individuals.txt  \
		--out $1.w10000.s500.n50 \
		--fst-window-size 10000 \
		--fst-window-step 500

vcftools --vcf $1.vcf \
		--weir-fst-pop ../vcftoolsConfigs/p1.$s_size.individuals.txt \
		--weir-fst-pop ../vcftoolsConfigs/p2.$s_size.individuals.txt  \
		--out $1.w5000.s500.n50 \
		--fst-window-size 5000 \
		--fst-window-step 500

vcftools --vcf $1.vcf \
		--weir-fst-pop ../vcftoolsConfigs/p1.$s_size.individuals.txt \
		--weir-fst-pop ../vcftoolsConfigs/p2.$s_size.individuals.txt  \
		--out $1.w50000.s500.n50 \
		--fst-window-size 50000 \
		--fst-window-step 500

vcftools --vcf $1.vcf \
		--weir-fst-pop ../vcftoolsConfigs/p1.$s_size.individuals.txt \
		--weir-fst-pop ../vcftoolsConfigs/p2.$s_size.individuals.txt  \
		--out $1.w100000.s500.n50 \
		--fst-window-size 100000 \
		--fst-window-step 500

vcftools --vcf $1.vcf \
		--weir-fst-pop ../vcftoolsConfigs/p1.$s_size.individuals.txt \
		--weir-fst-pop ../vcftoolsConfigs/p2.$s_size.individuals.txt  \
		--out $1.n$s_size

python3 /home/booker/work/FishersWaveFst/slim/bin/mutsFromTrees.py --input $1 --output $1.segSites.txt

gzip $1*fst

rm $1.vcf

gzip $1
