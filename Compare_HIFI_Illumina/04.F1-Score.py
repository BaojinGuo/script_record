import vcf
import sys

def calculate_precision_recall_f1(vcf_file, output_file):
    # Initialize counts for TP, FP, and FN
    tp_count = 0
    fp_count = 0
    fn_count = 0

    # Open the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        truth_sample = record.samples[0]
        test_sample = record.samples[1]

        truth_gt = truth_sample['GT']
        test_gt = test_sample['GT']

        # Check for non-missing genotype (1/1)
        if truth_gt == '1/1' and test_gt == '1/1':
            tp_count += 1
        elif truth_gt == '1/1' and test_gt == './.':
            fn_count += 1
        elif truth_gt == './.' and test_gt == '1/1':
            fp_count += 1

    # Calculate precision, recall, and F1-score
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    # Write the results to the output file
    with open(output_file, 'w') as out_f:
        out_f.write(f"Precision: {precision:.4f}\n")
        out_f.write(f"Recall: {recall:.4f}\n")
        out_f.write(f"F1-score: {f1_score:.4f}\n")

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <output_file>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_file = sys.argv[2]

    calculate_precision_recall_f1(input_vcf, output_file)
