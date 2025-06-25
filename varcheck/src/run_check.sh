
# Extract matching known variants
bcftools view -R checklist.tsv toy_patient.vcf -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > matched.tsv

