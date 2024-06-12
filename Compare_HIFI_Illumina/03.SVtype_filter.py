import argparse
import sys

def filter_sv(input_vcf, output_vcf):
    try:
        # Open the input VCF file
        with open(input_vcf, 'r') as f_input:
            # Open the output VCF file
            with open(output_vcf, 'w') as f_output:
                # Initialize line number
                line_num = 0

                # Iterate over each line in the input VCF file
                for line in f_input:
                    # Increment line number
                    line_num += 1

                    # Check if the line is a header line
                    if line.startswith('#'):
                        # Write the header line to the output VCF file
                        f_output.write(line)
                        continue

                    # Split the line into fields
                    fields = line.strip().split('\t')

                    # Check if there are enough fields
                    if len(fields) < 8:
                        print(f"Warning: Line {line_num} does not contain enough fields.")
                        continue

                    # Extract SV type and SV length from the INFO field
                    info_fields = fields[7].split(';')
                    svtype = None
                    svlen = None
                    for field in info_fields:
                        parts = field.split('=')
                        if len(parts) == 2:
                            key, value = parts
                            if key == 'SVTYPE':
                                svtype = value
                            elif key == 'SVLEN':
                                svlen = int(value)

                    # Check if SV type and SV length were found
                    if svtype is None or svlen is None:
                        print(f"Warning: Line {line_num} does not contain SVTYPE or SVLEN in the INFO field.")
                        continue

                    # Check if SV type is INS, DEL, INV, or DUP
                    if svtype in ['INS', 'DEL', 'INV', 'DUP']:
                        # Check if SV length is within the specified range
                        if -100000 < svlen < 100000:
                            f_output.write(line)
    except Exception as e:
        print("Error:", e)
        sys.exit(1)

if __name__ == "__main__":
    # Create a command-line parser
    parser = argparse.ArgumentParser(description='Filter VCF file for specific SV types and lengths.')
    parser.add_argument('input_vcf', help='Input VCF file path')
    parser.add_argument('output_vcf', help='Output VCF file path')
    args = parser.parse_args()

    # Run the filtering function
    filter_sv(args.input_vcf, args.output_vcf)

                                                        
