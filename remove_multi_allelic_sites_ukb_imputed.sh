#Goal: go from bgen files to a vcf file per chromosome, excluding variants with more than two alleles and keeping reference allele as described in the .bim file
#Imput: .bgen and plink (.bim) files for UK Biobank imputed data, list of samples one wishes to include.
#Output: compressed vcf file without multi-allelic SNPs and preserving REF allele identity.
#

for chr in {1..22];
do
grep -F -f tmp${chr} UKB/imputed/ukb_imp_chr${chr}_afr.bim > out${chr}.txt #tmp is a file with SNP IDs I want to look at. This makes this process faster since I don't need all variants.
awk '{print $2}' out${chr}.txt |sort|uniq -u > tmp_${chr}.txt #keep only bi-allelic SNPs
grep -F -f tmp_22.txt out22.txt |awk 'OFS="\t"{print $2,$5}' > ref_chr${chr}.txt #retain SNP ID and REFERENCE allele
#convert bgen to vcf for a list of samples, for the SNPs I am interested in, and keeping reference allele as described.
plink2 --bgen ukb_imp_chr${chr}_afr.bgen --sample my_samples.sample --extract tmp_${chr}.txt --ref-allele force ref_chr${chr}.txt --recode vcf --out chr${chr}
bgzip ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/chr${chr}.vcf
echo ${chr}
echo 'done'
done
