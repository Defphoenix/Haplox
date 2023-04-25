#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin '$Bin';

my ($outdir,$samplename,$cftype,$gtype,$cfDNAfq1,$cfDNAfq2,$gDNAfq1,$gDNAfq2,$bed,$bed2,$t,$s,$germ,$panel) = @ARGV;
=head1 Uasge

        perl  HapOnco_pipeline_pairV2.0.pl <outdir> <samplename> <cfDNAtype> <gDNAtype> <cfDNAfq1> <cfDNAfq2> <gDNAfq1> <gDNAfq2> <bed> <primary bed> <thread> <1 or 2 > <Y or N>  <451 or 605 or 605v2>
	
	<outdir>		output dir
	<samplename>		sample  prefix
	<cfDNAtype>		Duplex/nonUMI
	<gDNAtype>		Duplex/nonUMI
	<cfDNAfq1>		cfDNA fq1
	<cfDNAfq2>		cfDNA fq2
	<gDNAfq1>		gDNA fq1
	<gDNAfq2>		gDNA fq2
	<bed>			bed file  </thinker/net/ctDNA/Bed/605panel_v2/605panel_v2_capture_targets_sorted_merge_gene.bed>
	<primary bed>		bed file to stat QC</thinker/net/ctDNA/Bed/605panel_v2/605panel_v2_primary_targets_genes.bed>
	<thread>		thread use
	<genecore -s>		if  set 1 for all moles ;if  set 2 only for pcr moles 
	<germlline>		if run germline  set Y ; if you do not need run germline set N 
	<panel>			if 451panel set 451;if 605panel set 605; if 605panel set 605v2


=cut

die `pod2text $0` if (@ARGV != 14);

########################cfDNA#######################
system("mkdir -m 755 -p $outdir") if (!-d "$outdir");
system("mkdir -m 755 -p $outdir/mutscan") if (!-d "$outdir/mutscan");
system("mkdir -m 755 -p $outdir/fusion") if (!-d "$outdir/fusion");
system("mkdir -m 755 -p $outdir/cnv") if (!-d "$outdir/cnv");
system("mkdir -m 755 -p $outdir/bam") if (!-d "$outdir/bam");
system("mkdir -m 755 -p $outdir/shell") if (!-d "$outdir/shell");
system("mkdir -m 755 -p $outdir/somatic") if (!-d "$outdir/somatic");


#if ($cfDNAfq1=~/gdna/ && $cfDNAfq2=~/gdna/){print "cfdna  or ffpe fastq is not right!\n";die;}# add by chenyr2018/9/14


open OUT1,">$outdir/shell/$samplename\_cfDNA.sh";
print OUT1 "export SENTIEON_LICENSE=192.168.1.29:8990\n";
if($cftype eq "nonUMI"){
	print OUT1 "$Bin/fastp -f 1 -F 1 -t 1 -T 1 -3 -W 4 -M 25 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $cfDNAfq1 -o $outdir/$samplename\_cfDNA_good_R1.fastq.gz -I $cfDNAfq2 -O $outdir/$samplename\_cfDNA_good_R2.fastq.gz  -h $outdir/$samplename\_cfDNA_fastp.html -j $outdir/$samplename\_cfDNA_fastp.json > $outdir/$samplename\_cfDNA_fastp.stat 2>&1\n";
}elsif($cftype eq "Duplex"){
	print OUT1 "$Bin/fastp --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $cfDNAfq1 -o $outdir/$samplename\_cfDNA_good_R1.fastq.gz -I $cfDNAfq2 -O $outdir/$samplename\_cfDNA_good_R2.fastq.gz  -h $outdir/$samplename\_cfDNA_fastp.html -j $outdir/$samplename\_cfDNA_fastp.json  -U --umi_loc=per_read --umi_skip=4 --umi_prefix=UMI  --umi_len=3  > $outdir/$samplename\_cfDNA_fastp.stat 2>&1\n";
}
# print OUT1 "nohup $Bin/mutscan -1 $outdir/$samplename\_cfDNA_good_R1.fastq.gz -2 $outdir/$samplename\_cfDNA_good_R2.fastq.gz -h $outdir/mutscan/$samplename\_mutscan.html >$outdir/mutscan/$samplename\_mutscan.log 2>&1 &\n";

print OUT1 "nohup $Bin/genefuse -t $t -r $Bin/WES_ref/hg19.fa -1 $outdir/$samplename\_cfDNA_good_R1.fastq.gz -2 $outdir/$samplename\_cfDNA_good_R2.fastq.gz -f $Bin/GeneFuse/genes/druggable.hg19.csv -h $outdir/fusion/$samplename\_fusion.html  -j $outdir/fusion/$samplename\_fusion.json >$outdir/fusion/$samplename.fusion.log 2>&1  && perl $Bin/genefuse.filter_v2.pl  -c  $Bin/Comicfusiongenepair.txt -b  $Bin/genefuse.gene.list  -t $outdir/fusion/$samplename\_fusion.html -o $outdir/fusion/$samplename\_fusion.result &\n"; ## change by chenyr 2018/11/1

print OUT1 "$Bin/sentieon-genomics-201808.05/bin/bwa mem -k 32 -R \"\@RG\\tID:$samplename\_tumor\\tLB:$samplename\_tumor\\tPL:ILLUMINA\\tSM:$samplename\_tumor\" -t $t -M $Bin/WES_ref/hg19.fa $outdir/$samplename\_cfDNA_good_R1.fastq.gz $outdir/$samplename\_cfDNA_good_R2.fastq.gz| $Bin/sentieon-genomics-201808.05/bin/sentieon  util sort -o $outdir/bam/$samplename.cfDNA.sort.bam -t $t --sam2bam -i - \n";

if($cftype eq "Duplex"){
	print OUT1 "$Bin/gencore -i  $outdir/bam/$samplename.cfDNA.sort.bam  -o  $outdir/bam/$samplename.cfDNA.gencore.sort.bam -r $Bin/WES_ref/hg19.fa -u  UMI -s $s -j $outdir/bam/$samplename.cfDNA.gencore.json  -h $outdir/bam/$samplename.cfDNA.gencore.html \n";
}elsif($cftype eq "nonUMI"){
	print OUT1 "$Bin/gencore -i  $outdir/bam/$samplename.cfDNA.sort.bam  -o  $outdir/bam/$samplename.cfDNA.gencore.sort.bam -r $Bin/WES_ref/hg19.fa   -s $s -j $outdir/bam/$samplename.cfDNA.gencore.json -h $outdir/bam/$samplename.cfDNA.gencore.html \n";
}
print OUT1 "$Bin/samtools-1.3/bin/samtools index $outdir/bam/$samplename.cfDNA.gencore.sort.bam\n";
print OUT1 "wait\n";
print OUT1 "$Bin/samtools-1.3/bin/samtools depth -aa -q 25 -d 100000 -b $bed2 $outdir/bam/$samplename.cfDNA.gencore.sort.bam > $outdir/bam/$samplename.cfDNA.gencore.sort.depth & \n";
print OUT1 "$Bin/samtools-0.1.19/samtools mpileup -AB -Q 25 -q 30  -d 100000 -f $Bin/WES_ref/hg19.fa -l $bed  $outdir/bam/$samplename.cfDNA.gencore.sort.bam >$outdir/bam/$samplename.cfDNA.gencore.mpileup & \n";
print OUT1 "rm $outdir/bam/$samplename.cfDNA.sort.bam  & \n";
##########################
print OUT1 "wait\n";
print OUT1 "perl $Bin/stat_umi4.pl $outdir/$samplename\_cfDNA_fastp.stat  $outdir/bam/$samplename.cfDNA.gencore.json $outdir/bam/$samplename.cfDNA.gencore.mpileup $outdir/bam/$samplename.cfDNA.gencore.sort.depth >$outdir/$samplename.cfDNA.stat\n";
print OUT1 "perl $Bin/script/QC_csv.pl $cfDNAfq1 $outdir/$samplename.cfDNA.stat >  $outdir/$samplename.cfDNA.QC.csv \n"; #add QC file by xuyun@20190123
print OUT1 "curl haplab.haplox.net/api/report/depth-new -F \"import_file=\@$outdir/$samplename.cfDNA.QC.csv\"\n"; #add upload by chenyr@2019/9/18
##########################

print OUT1 "echo;date\n";
close OUT1;
######################################################################

open OUT2,">$outdir/shell/$samplename\_gDNA.sh";
print OUT2 "export SENTIEON_LICENSE=192.168.1.29:8990\n";
if($gtype eq "nonUMI"){
	print OUT2 "$Bin/fastp -f 1 -F 1 -t 1 -T 1 -3 -W 4 -M 25 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $gDNAfq1 -o $outdir/$samplename\_gDNA_good_R1.fastq.gz -I $gDNAfq2 -O $outdir/$samplename\_gDNA_good_R2.fastq.gz -h $outdir/$samplename\_gDNA_fastp.html -j $outdir/$samplename\_gDNA_fastp.json > $outdir/$samplename\_gDNA_fastp.stat 2>&1\n";
}elsif($gtype eq "Duplex"){
	print OUT2 "$Bin/fastp --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $gDNAfq1 -o $outdir/$samplename\_gDNA_good_R1.fastq.gz -I $gDNAfq2 -O $outdir/$samplename\_gDNA_good_R2.fastq.gz -h $outdir/$samplename\_gDNA_fastp.html -j $outdir/$samplename\_gDNA_fastp.json  -U --umi_loc=per_read --umi_skip=4 --umi_prefix=UMI  --umi_len=3 > $outdir/$samplename\_gDNA_fastp.stat 2>&1\n";
}

print OUT2 "$Bin/sentieon-genomics-201808.05/bin/bwa mem -k 32 -R \"\@RG\\tID:$samplename\_normal\\tLB:$samplename\_normal\\tPL:ILLUMINA\\tSM:$samplename\_normal\" -t $t -M $Bin/WES_ref/hg19.fa $outdir/$samplename\_gDNA_good_R1.fastq.gz $outdir/$samplename\_gDNA_good_R2.fastq.gz  | $Bin/sentieon-genomics-201808.05/bin/sentieon  util sort -o $outdir/bam/$samplename.gDNA.sort.bam -t $t --sam2bam -i - \n";

if($gtype eq "Duplex"){
       print OUT2 "$Bin/gencore -i  $outdir/bam/$samplename.gDNA.sort.bam  -o  $outdir/bam/$samplename.gDNA.gencore.sort.bam -r $Bin/WES_ref/hg19.fa -u UMI -s 1 -j $outdir/bam/$samplename.gDNA.gencore.json  -h $outdir/bam/$samplename.gDNA.gencore.html \n";
}elsif($gtype eq "nonUMI"){
        print OUT2 "$Bin/gencore -i  $outdir/bam/$samplename.gDNA.sort.bam  -o  $outdir/bam/$samplename.gDNA.gencore.sort.bam -r $Bin/WES_ref/hg19.fa  -s 1 -j $outdir/bam/$samplename.gDNA.gencore.json  -h $outdir/bam/$samplename.gDNA.gencore.html \n";
}
print OUT2 "$Bin/samtools-1.3/bin/samtools index $outdir/bam/$samplename.gDNA.gencore.sort.bam\n";
print OUT2 "wait\n";
print OUT2 "$Bin/samtools-1.3/bin/samtools depth -aa -q 25 -d 100000 -b   $bed2  $outdir/bam/$samplename.gDNA.gencore.sort.bam > $outdir/bam/$samplename.gDNA.gencore.sort.depth & \n";
print OUT2 "$Bin/samtools-0.1.19/samtools mpileup -AB -Q 25 -q 30 -d 100000 -f $Bin/WES_ref/hg19.fa -l $bed  $outdir/bam/$samplename.gDNA.gencore.sort.bam >$outdir/bam/$samplename.gDNA.gencore.mpileup &\n";
print OUT2 "rm  $outdir/bam/$samplename.gDNA.sort.bam  & \n";
###################
print OUT2 "wait\n";
print OUT2 "perl $Bin/stat_umi4.pl $outdir/$samplename\_gDNA_fastp.stat  $outdir/bam/$samplename.gDNA.gencore.json  $outdir/bam/$samplename.gDNA.gencore.mpileup $outdir/bam/$samplename.gDNA.gencore.sort.depth  >$outdir/$samplename.gDNA.stat\n";
print OUT2 "perl $Bin/script/QC_csv.pl $gDNAfq1 $outdir/$samplename.gDNA.stat >  $outdir/$samplename.gDNA.QC.csv \n"; #add QC file by xuyun@20190123
print OUT2 "curl haplab.haplox.net/api/report/depth-new -F \"import_file=\@$outdir/$samplename.gDNA.QC.csv\"\n"; #add upload by chenyr@2019/9/18
####################
print OUT2 "echo;date\n";
close OUT2;

########################Germline####################################
my $dbsnp="$Bin/germline_c/dbsnp_138.hg19.vcf";
my $millsIndels="$Bin/germline_c/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
my $oneKgIndels="$Bin/germline_c/1000G_phase1.indels.hg19.sites.vcf";



if ($germ eq "Y") {
	open OUT_G,">$outdir/shell/$samplename\_germline.sh";
	print  OUT_G "mkdir -m 755 -p $outdir/germline/tmp\n";
	print  OUT_G "mkdir -m 755 -p $outdir/germline/bamdst_stat\n";
	print OUT_G "export SENTIEON_LICENSE=192.168.1.29:8990\n";
	print OUT_G "$Bin/bamdst   -p $bed2 -o  $outdir/germline/bamdst_stat   $outdir/bam/$samplename.gDNA.gencore.sort.bam  &\n" ;

	print OUT_G  "$Bin/sentieon-genomics-201808.05/bin/sentieon driver -r  $Bin/WES_ref/hg19.fa  -t $t -i  $outdir/bam/$samplename.gDNA.gencore.sort.bam --algo Realigner -k $oneKgIndels -k $millsIndels   --interval_list   $bed $outdir/bam/$samplename.normal_realigned.bam 2>$outdir/germline/tmp/${samplename}_gatk.log\n" ;

	print OUT_G  "$Bin/sentieon-genomics-201808.05/bin/sentieon driver -t $t -r  $Bin/WES_ref/hg19.fa  -i  $outdir/bam/$samplename.normal_realigned.bam     --interval  $bed --algo Haplotyper --emit_conf  30.0  --call_conf  10      $outdir/germline/${samplename}.raw.vcf 2>>$outdir/germline/tmp/${samplename}_gatk.log\n"    ;

	print OUT_G  "java -Djava.io.tmpdir=$outdir/germline/tmp/${samplename}_gatk    -jar $Bin/GATK/GenomeAnalysisTK.jar    -R  $Bin/WES_ref/hg19.fa    -T SelectVariants   -V $outdir/germline/${samplename}.raw.vcf   -selectType SNP    -o $outdir/germline/${samplename}.raw.snp.vcf   2>>$outdir/germline/tmp/${samplename}_gatk.log\n" ;
	print OUT_G  "java -Djava.io.tmpdir=$outdir/germline/tmp/${samplename}_gatk    -jar $Bin/GATK/GenomeAnalysisTK.jar     -R  $Bin/WES_ref/hg19.fa    -T SelectVariants   -V $outdir/germline/${samplename}.raw.vcf   -selectType INDEL     -o $outdir/germline/${samplename}.raw.indel.vcf    2>>$outdir/germline/tmp/${samplename}_gatk.log\n"  ;
	print OUT_G   "java -Djava.io.tmpdir=$outdir/germline/tmp/${samplename}_gatk     -jar $Bin/GATK/GenomeAnalysisTK.jar     -R  $Bin/WES_ref/hg19.fa    -T VariantFiltration   -o $outdir/germline/${samplename}.raw.snp-tmp.vcf  -V $outdir/germline/${samplename}.raw.snp.vcf    -filter \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"     -filterName \"SNP_filter\"       2>>$outdir/germline/tmp/${samplename}_gatk.log\n"    ;
	print OUT_G   "java -Djava.io.tmpdir=$outdir/germline/tmp/${samplename}_gatk   -jar $Bin/GATK/GenomeAnalysisTK.jar  -R  $Bin/WES_ref/hg19.fa   -T VariantFiltration  -o $outdir/germline/${samplename}.raw.indel-tmp.vcf  -V $outdir/germline/${samplename}.raw.indel.vcf  -filter \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\"  -filterName \"INDEL_filter\"  2>>$outdir/germline/tmp/${samplename}_gatk.log\n";
	print OUT_G   "grep '^#' $outdir/germline/${samplename}.raw.snp-tmp.vcf > $outdir/germline/${samplename}.filter.vcf \n"  ;
	print OUT_G   "grep PASS $outdir/germline/${samplename}.raw.snp-tmp.vcf >> $outdir/germline/${samplename}.filter.vcf \n"  ;
	print OUT_G   "grep PASS $outdir/germline/${samplename}.raw.indel-tmp.vcf>> $outdir/germline/${samplename}.filter.vcf\n"  ;
	print OUT_G  "$Bin/annovar/table_annovar.pl $outdir/germline/${samplename}.filter.vcf  $Bin/annovar/humandb/ -buildver hg19  -out   $outdir/germline/${samplename}.filter     -remove    -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,snp138,ljb26_all    -operation g,r,r,f,f,f,f,f,f,f     -arg '-splicing_threshold 2',,,,,,,,,     -nastring .    -vcfinput    2>>$outdir/germline/tmp/${samplename}.annovar_snp.log   \n"    ;
	print OUT_G   "perl $Bin/germline_c/chem_sen.pl  $outdir/germline/${samplename}.filter.vcf  $Bin/germline_c/Chem_drug_sen.list  $outdir/germline/${samplename}.chem_genotype.txt\n"    ;
	print OUT_G   "perl $Bin/germline_c/NewPan_chem.pl $outdir/germline/${samplename}.filter.vcf $Bin/germline_c/Chem451_gene.txt  $outdir/germline/${samplename}.chem_451.txt\n"    ;
	print OUT_G   "perl $Bin/germline_c/NewPan_chem.pl $outdir/germline/${samplename}.filter.vcf $Bin/germline_c/Target451_gene.txt  $outdir/germline/${samplename}.Target_451.txt\n"    ;
	print OUT_G   "perl $Bin/germline_c/PathCall.pl -i $outdir/germline/${samplename}.filter.hg19_multianno.txt -g $Bin/germline_c/Heart_Cancer_Incidental.txt  -o  $outdir/germline/${samplename}.germline.txt\n"    ;
	print OUT_G   "perl $Bin/germline_c/germline_trans.pl   $outdir/germline/${samplename}.germline.txt   $outdir/germline/${samplename}.trans_germ.txt\n"    ;
	print OUT_G  "perl $Bin/germline_c/Gene2Disease.pl    $Bin/germline_c/female_cancer.list  $outdir/germline/${samplename}.germline.txt  $outdir/germline/${samplename}.cancer.female.txt   $outdir/germline/${samplename}.nocancer.female.txt  \n" ;
	print OUT_G  "perl $Bin/germline_c/Gene2Disease.pl    $Bin/germline_c/male_cancer.list  $outdir/germline/${samplename}.germline.txt  $outdir/germline/${samplename}.cancer.male.txt   $outdir/germline/${samplename}.nocancer.male.txt \n" ;
	print OUT_G  "perl $Bin/germline_c/germline_trans_v2.pl  $outdir/germline/${samplename}.cancer.male.txt  $outdir/germline/${samplename}.trans.cancer.male.txt\n" ;
	print OUT_G  "perl $Bin/germline_c/germline_trans_v2.pl  $outdir/germline/${samplename}.cancer.female.txt  $outdir/germline/${samplename}.trans.cancer.female.txt\n" ;
	print OUT_G  "perl $Bin/germline_c/QC_exome.v2.pl   $outdir/germline/bamdst_stat  $outdir/bam/$samplename.gDNA.gencore.json   $samplename \n";

}


#########################SOMATIC##################################

open OUT3,">$outdir/shell/$samplename\_somatic.sh";
print OUT3 "echo;date\n";
print OUT3 "perl$Bin/Bed/add_geneinfo_to_bed.pl -i$Bin/GRCh37.p13_2017_01_13.gff  -i1 $Bin/Gene_transcript_new.list -i2  $bed -o  $outdir/cnv/analysis.bed   \n";
print OUT3 "python $Bin/cnvkit-0.9.6/cnvkit-0.9.6/cnvkit.py  batch  $outdir/bam/$samplename.cfDNA.gencore.sort.bam   -n  $outdir/bam/$samplename.gDNA.gencore.sort.bam -t $outdir/cnv/analysis.bed  -f   $Bin/WES_ref/hg19.fa   --access   $Bin/access-5k-mappable.hg19.bed  -d  $outdir/cnv  > $outdir/$samplename.cnv.log 2>&1 && perl $Bin/cnv_filter_wesplus.pl   -l $Bin/final.cnv.genelist  -c $outdir/cnv/${samplename}.cfDNA.gencore.sort.cns  -o $outdir/cnv/${samplename}.cnv.result && perl /haplox/users/xuyun/script/result_extract.pl -i $outdir -s $samplename && rm    $outdir/cnv/analysis.bed & \n";#add extract CNV and fusion result by xuyun@2018/12/04
print OUT3 "java -jar $Bin/VarScan.v2.4.3.jar somatic $outdir/bam/$samplename.gDNA.gencore.mpileup $outdir/bam/$samplename.cfDNA.gencore.mpileup --output-snp $outdir/somatic/$samplename.snp.vcf --output-indel $outdir/somatic/$samplename.indel.vcf --min-coverage-tumor 80 --min-coverage-normal 60 --min-var-freq 0.0001 --min-freq-for-hom 0.90 --normal-purity 1 --tumor-purity 0.01 --p-value 0.99 --somatic-p-value 0.05 --strand-filter 1 --output-vcf 1\n";
print OUT3 "$Bin/annovar/table_annovar.pl $outdir/somatic/$samplename.snp.vcf $Bin/annovar/humandb/ -buildver hg19 -out $outdir/somatic/$samplename.snp -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,snp138,ljb26_all,cosmic87,clinvar_20180603 -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n";
print OUT3 "$Bin/annovar/table_annovar.pl $outdir/somatic/$samplename.indel.vcf $Bin/annovar/humandb/ -buildver hg19 -out $outdir/somatic/$samplename.indel -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,snp138,ljb26_all,cosmic87,clinvar_20180603 -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n";
#sample check with germline Homozygous SNV# By xuyun@2018.11.29
print OUT3 "perl $Bin/check_snp.v1.pl $outdir/somatic/$samplename.snp.hg19_multianno.txt > $outdir/somatic/$samplename.check.result\n";
print OUT3 "perl $Bin/check.pl $outdir/somatic/$samplename.snp.hg19_multianno.txt > $outdir/somatic/$samplename.check.contaminate\n";
print OUT3 "perl $Bin/filter_pair_step1.pl -i $outdir/somatic/$samplename.snp.hg19_multianno.txt -o $outdir/somatic/$samplename.snp.hg19_multianno.txt.filter\n";
print OUT3 "perl $Bin/filter_pair_step1.pl -i $outdir/somatic/$samplename.indel.hg19_multianno.txt -o $outdir/somatic/$samplename.indel.hg19_multianno.txt.filter\n";
print OUT3 "$Bin/MrBam_new/paraMrBam $t $outdir/somatic/$samplename.snp.hg19_multianno.txt.filter  --cfdna $outdir/bam/$samplename.cfDNA.gencore.sort.bam --gdna $outdir/bam/$samplename.gDNA.gencore.sort.bam -m 2 -q 25 --fast --snp >$outdir/somatic/$samplename.snp_Test_MrBam.txt\n";
print OUT3 "$Bin/MrBam_new/paraMrBam $t $outdir/somatic/$samplename.indel.hg19_multianno.txt.filter  --cfdna $outdir/bam/$samplename.cfDNA.gencore.sort.bam --gdna $outdir/bam/$samplename.gDNA.gencore.sort.bam -m 2 -q 25 --fast --indel >$outdir/somatic/$samplename.indel_Test_MrBam.txt\n";
#print OUT3 "$Bin/baselineanno -i $outdir/somatic/$samplename.snp_Test_MrBam.txt -o $outdir/somatic/$samplename.snp_Test_MrBam.baseline -t 1 -redis 192.168.1.10:6379 -anno 1\n"; #this is old health_cfdna baseline
#print OUT3 "$Bin/baselineanno -i $outdir/somatic/$samplename.indel_Test_MrBam.txt -o $outdir/somatic/$samplename.indel_Test_MrBam.baseline -t 1 -redis 192.168.1.10:6379 -anno 1\n"; #this is old health_cfdna baseline
print OUT3 "$Bin/baselineanno4 -i $outdir/somatic/$samplename.snp_Test_MrBam.txt -o $outdir/somatic/$samplename.snp_Test_MrBam.baseline -t 1 -redis  Haplox2022*#\@10.6.12.2:6379 -anno 6\n"; # change db ip 2023/2/24
print OUT3 "$Bin/baselineanno4 -i $outdir/somatic/$samplename.indel_Test_MrBam.txt -o $outdir/somatic/$samplename.indel_Test_MrBam.baseline -t 1 -redis  Haplox2022*#\@10.6.12.2:6379 -anno 6\n"; # change db ip 2023/2/24
print OUT3 "$Bin/extract.v2.pl $outdir/somatic/$samplename.snp_Test_MrBam.baseline $outdir/somatic/$samplename.snp_Test.txt\n";
print OUT3 "$Bin/extract.v2.pl $outdir/somatic/$samplename.indel_Test_MrBam.baseline $outdir/somatic/$samplename.indel_Test.txt\n";
#add 451 or 605panel genelist filter# by xuyun@2019.01.09
print OUT3 "$Bin/filter2_v2.py true $panel $outdir/somatic/$samplename.snp_Test.txt $outdir/somatic/$samplename.indel_Test.txt $outdir/somatic/$samplename.filter_Test.log >$outdir/somatic/$samplename.filter_Test.result\n";
#605panle baseline annovation# By xuyun@2018.11.29
print OUT3 "perl $Bin/baseline_anno_v2.pl $outdir/somatic/$samplename.filter_Test.result $s >  $outdir/somatic/$samplename.filter_Test.result.baseline\n";
print OUT3 "rm $outdir/bam/$samplename.cfDNA.gencore.mpileup  $outdir/bam/*DNA.gencore.bam $outdir/*gz &&  rm $outdir/bam/*depth "; # add by chenyr2018/8/9
close OUT3;

open OUT4,">$outdir/shell/$samplename\_run.sh";
print OUT4 "###################################################################\n";
print OUT4 "#Pipeline\t/thinker/net/ctDNA/HapOnco_pipeline_pairV2.0.pl\n";
print OUT4 "#Version\tV2.0\n";
print OUT4 "#Update\t\t2020/2/28\n";
print OUT4 "#Author\t\tchenyr\n";
print OUT4 "#Email\t\tchenyr\@haplox.com\n";
print OUT4 "###################################################################\n";
print OUT4 "#Run  Command\n";
print OUT4 "#perl  /thinker/net/ctDNA/HapOnco_pipeline_pairV2.0.pl $outdir $samplename $cftype $gtype $cfDNAfq1 $cfDNAfq2 $gDNAfq1 $gDNAfq2 $bed $bed2 $t $s $germ $panel\n";
print OUT4 "###################################################################\n";
print OUT4 "nohup sh $outdir/shell/$samplename\_gDNA.sh >$outdir/shell/$samplename\_gDNA.sh.o 2>&1 &\n";
print OUT4 "nohup sh $outdir/shell/$samplename\_cfDNA.sh >$outdir/shell/$samplename\_cfDNA.sh.o 2>&1 &\n";
print OUT4 "wait\n";
print OUT4 "sh $outdir/shell/$samplename\_somatic.sh >$outdir/shell/$samplename\_somatic.sh.o 2>&1 & \n";
if ($germ eq "Y") {
	print OUT4 "sh $outdir/shell/$samplename\_germline.sh >$outdir/shell/$samplename\_germline.sh.o 2>&1 &\n";
}
# Test for New MrBam # By wangzy@haplox.com
close OUT4;

print "cd $outdir\n";
print "nohup sh $outdir/shell/$samplename\_run.sh >$outdir/shell/$samplename\_run.sh.o 2>&1 &\n";
