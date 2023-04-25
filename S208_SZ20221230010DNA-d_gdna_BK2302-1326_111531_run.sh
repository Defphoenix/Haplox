###################################################################
#Pipeline	/thinker/net/ctDNA/HapOnco_pipeline_pairV2.0.pl
#Version	V2.0
#Update		2020/2/28
#Author		chenyr
#Email		chenyr@haplox.com
###################################################################
#Run  Command
#perl  /thinker/net/ctDNA/HapOnco_pipeline_pairV2.0.pl /haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531 S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531 nonUMI nonUMI /haplox/users/chenyr/Project/1326Panel/rawfq/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_S10_L003_R1_001.fastq.gz /haplox/users/chenyr/Project/1326Panel/rawfq/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_S10_L003_R2_001.fastq.gz /haplox/users/chenyr/Project/1326Panel/rawfq/S209_SZ20221230009DNA-a_gdna_BK2302-1326_111532_S11_L003_R1_001.fastq.gz /haplox/users/chenyr/Project/1326Panel/rawfq/S209_SZ20221230009DNA-a_gdna_BK2302-1326_111532_S11_L003_R2_001.fastq.gz /haplox/users/chenyr/Panel_Design/HapOnco1326/HapOnco1326_probe.exp100.bed /haplox/users/chenyr/Panel_Design/HapOnco1326/HapOnco1326_probe.bed 10 1 N 605v2
###################################################################
nohup sh /haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531/shell/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_gDNA.sh >/haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531/shell/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_gDNA.sh.o 2>&1 &
nohup sh /haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531/shell/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_cfDNA.sh >/haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531/shell/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_cfDNA.sh.o 2>&1 &
wait
sh /haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531/shell/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_somatic.sh >/haplox/users/chenyr/Project/1326Panel/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531/shell/S208_SZ20221230010DNA-d_gdna_BK2302-1326_111531_somatic.sh.o 2>&1 & 
