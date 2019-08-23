## A li'l pipeline to call SNPs from trees and to use VCFtools to do a widowed analysis on the data

##	for i in $(ls */)
##		do

w_size=10000
w_step=10000
s_size=50


## Keep an eye on the mutation rate in the recapSprinkler.py script, it varies by a factor of 2.5x in the demoPlot runs compared to the big batch
python3 /home/booker/work/FishersWaveFst/slim/bin/recapSprinkler.py \
	--input $1  \
	--output $1.vcf \
	--sample $s_size \
	--just_vcf

vcftools --vcf $1.vcf \
		--weir-fst-pop ../vcftoolsConfigs/p1.$s_size.individuals.txt \
		--weir-fst-pop ../vcftoolsConfigs/p2.$s_size.individuals.txt  \
		--out $1.w10000.s500.n50 \
		--fst-window-size 10000 \
		--fst-window-step 500

python3 /home/booker/work/FishersWaveFst/slim/bin/mutsFromTrees.py --input $1 --output $1.segSites.txt

gzip $1*fst

rm $1.vcf