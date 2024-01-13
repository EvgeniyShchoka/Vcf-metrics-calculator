# VCF Metrics Calculator

This Python script calculates various metrics from a Variant Call Format (VCF) file and generates an updated VCF file with additional information. Additionally, it creates a JSON file summarizing the calculated metrics.

## Files

- `input.vcf`: Input VCF file containing genetic variant information.
- `output.vcf`: Output VCF file with added Allele Frequency (AF) fields.
- `output.json`: Output JSON file containing calculated metrics.

## Metrics Calculated

1. **Allele Frequency (AF):**
   - AF is calculated for each variant and added as a new field in the output VCF file.

2. **Chromosome-wise Variant Count:**
   - Counts the number of variants per chromosome.

3. **Heterozygous to Homozygous Ratio:**
   - Calculates the ratio of heterozygous to homozygous variants.

4. **Transition to Transversion Ratio (Ti/Tv):**
   - Calculates the ratio of transitions (A<>G or C<>T) to transversions (A<>C, A<>T, G<>C, G<>T).

## Usage

```bash
python3 generate_vcf_metrics.py input.vcf output.vcf output.json

