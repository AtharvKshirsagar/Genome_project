import hail as hl
import matplotlib.pyplot as plt

# Initialize Hail
hl.init()

# Load the VCF file
vcf_file = "your_vcf_file.vcf"
mt = hl.import_vcf(vcf_file)

# Calculate allele frequencies
mt = mt.annotate_rows(allele_freq=hl.agg.call_stats(mt.GT).AF[1])

# Collect the top 10 variants with the highest allele frequencies
top_variants = mt.rows().key_by().top(10, mt.allele_freq, reverse=True)

# Print and plot the top variants
for variant in top_variants:
    print(f"Variant: {variant.locus}:{variant.alleles} - Allele Frequency: {variant.allele_freq}")

# Create a bar plot of allele frequencies
allele_freqs = top_variants.variant.collect()
variant_labels = [f"{v.locus}:{v.alleles}" for v in allele_freqs]
allele_frequencies = top_variants.allele_freq.collect()

plt.barh(variant_labels, allele_frequencies)
plt.xlabel("Allele Frequency")
plt.ylabel("Variant")
plt.title("Top 10 Variants with Highest Allele Frequencies")
plt.show()

# Stop the Hail session
hl.stop()
