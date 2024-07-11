import vcf
import sys

def calculate_precision_recall_f1(vcf_file, output_file):
    # Initialize counts for SNP and INDEL
    snp_tp = 0
    snp_fp = 0
    snp_fn = 0
    indel_tp = 0
    indel_fp = 0
    indel_fn = 0

    # Open the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        truth_sample = record.samples[0]
        test_sample = record.samples[1]

        truth_gt = truth_sample['GT']
        test_gt = test_sample['GT']

        # Determine variant type (SNP or INDEL)
        if len(record.REF) == len(record.ALT[0]):
            variant_type = 'SNP'
        else:
            variant_type = 'INDEL'

        # Calculate TP, FP, FN for SNP and INDEL separately
        if truth_gt == '1/1' and test_gt == '1/1':
            if variant_type == 'SNP':
                snp_tp += 1
            elif variant_type == 'INDEL':
                indel_tp += 1
        elif truth_gt == '1/1' and test_gt == './.':
            if variant_type == 'SNP':
                snp_fn += 1
            elif variant_type == 'INDEL':
                indel_fn += 1
        elif truth_gt == './.' and test_gt == '1/1':
            if variant_type == 'SNP':
                snp_fp += 1
            elif variant_type == 'INDEL':
                indel_fp += 1

    # Calculate precision, recall, and F1-score for SNP
    snp_precision = snp_tp / (snp_tp + snp_fp) if (snp_tp + snp_fp) > 0 else 0
    snp_recall = snp_tp / (snp_tp + snp_fn) if (snp_tp + snp_fn) > 0 else 0
    snp_f1_score = 2 * (snp_precision * snp_recall) / (snp_precision + snp_recall) if (snp_precision + snp_recall) > 0 else 0

    # Calculate precision, recall, and F1-score for INDEL
    indel_precision = indel_tp / (indel_tp + indel_fp) if (indel_tp + indel_fp) > 0 else 0
    indel_recall = indel_tp / (indel_tp + indel_fn) if (indel_tp + indel_fn) > 0 else 0
    indel_f1_score = 2 * (indel_precision * indel_recall) / (indel_precision + indel_recall) if (indel_precision + indel_recall) > 0 else 0

    # Write the results to the output file
    with open(output_file, 'w') as out_f:
        out_f.write("SNP Metrics:\n")
        out_f.write(f"Precision: {snp_precision:.4f}\n")
        out_f.write(f"Recall: {snp_recall:.4f}\n")
        out_f.write(f"F1-score: {snp_f1_score:.4f}\n\n")
        
        out_f.write("INDEL Metrics:\n")
        out_f.write(f"Precision: {indel_precision:.4f}\n")
        out_f.write(f"Recall: {indel_recall:.4f}\n")
        out_f.write(f"F1-score: {indel_f1_score:.4f}\n")

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <output_file>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_file = sys.argv[2]

    calculate_precision_recall_f1(input_vcf, output_file)
