import sys
import argparse

def process_vcf(vcf_file):
    """
    Processes the input VCF file to count SNPs and INDELs, 
    and calculate INDEL lengths.

    Args:
        vcf_file (str): Path to the input VCF file.

    Returns:
        dict: A dictionary containing overall statistics and sample-specific statistics.
    """
    with open(vcf_file, 'r') as file:
        lines = file.readlines()

    # Initialize counters for overall statistics
    snp_count = 0
    indel_count = 0
    indel_lengths = []

    # Sample-specific counters
    samples = ["bcftools", "Deepvariant", "GATK"]
    sample_snp_counts = {sample: 0 for sample in samples}
    sample_indel_counts = {sample: 0 for sample in samples}
    sample_indel_lengths = {sample: [] for sample in samples}

    # Process each line in the VCF file
    for line in lines:
        if line.startswith('#'):
            continue
        fields = line.split()
        ref = fields[2]
        alt = fields[3]
        genotypes = fields[4:]

        # Determine if the variant is SNP or INDEL
        if len(ref) == len(alt) == 1:
            snp_count += 1
            for i, genotype in enumerate(genotypes):
                if genotype != "./.":
                    sample_snp_counts[samples[i]] += 1
        else:
            indel_count += 1
            indel_length = len(alt) - len(ref)
            indel_lengths.append(indel_length)
            for i, genotype in enumerate(genotypes):
                if genotype != "./.":
                    sample_indel_counts[samples[i]] += 1
                    sample_indel_lengths[samples[i]].append(indel_length)

    return {
        "total_snps": snp_count,
        "total_indels": indel_count,
        "indel_lengths": indel_lengths,
        "sample_snp_counts": sample_snp_counts,
        "sample_indel_counts": sample_indel_counts,
        "sample_indel_lengths": sample_indel_lengths
    }

def write_statistics(output_file, data, sample_name=None):
    """
    Writes the statistics to an output file.

    Args:
        output_file (str): Path to the output file.
        data (dict): Dictionary containing the statistics data.
        sample_name (str, optional): Name of the sample. If None, writes overall statistics.
    """
    with open(output_file, 'w') as file:
        if sample_name:
            file.write(f"\n{sample_name} statistics:\n")
            file.write(f"SNPs: {data['sample_snp_counts'][sample_name]}\n")
            file.write(f"INDELs: {data['sample_indel_counts'][sample_name]}\n")
            file.write(f"INDEL lengths: {data['sample_indel_lengths'][sample_name]}\n")
        else:
            file.write(f"Total SNPs: {data['total_snps']}\n")
            file.write(f"Total INDELs: {data['total_indels']}\n")
            file.write(f"INDEL lengths distribution: {data['indel_lengths']}\n")

def main():
    """
    Main function to parse command-line arguments and process the VCF file.
    Generates output files with overall and sample-specific statistics.
    """
    parser = argparse.ArgumentParser(description='Process VCF file and generate statistics.')
    parser.add_argument('vcf_file', help='Input VCF file')
    parser.add_argument('output_prefix', help='Output file prefix')

    args = parser.parse_args()

    data = process_vcf(args.vcf_file)
    
    # Write overall statistics
    write_statistics(f"{args.output_prefix}_overall.txt", data)

    # Write statistics for each sample
    samples = ["bcftools", "Deepvariant", "GATK"]
    for sample in samples:
        write_statistics(f"{args.output_prefix}_{sample}.txt", data, sample)

if __name__ == "__main__":
    main()
