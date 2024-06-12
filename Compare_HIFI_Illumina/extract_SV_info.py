###Usage: python script.py input.vcf output.summary
import sys

def extract_info(input_file, output_file):
    """
    Extracts information from a VCF file and writes it to an output file.

Args:
    input_file (str): The name of the input VCF file.
    output_file (str): The name of the output file.

    The extracted information includes chromosome, position, SVLEN, SVTYPE,
    and the first part of the information in columns 10 to the last column
    of the VCF file. Lines with SVLEN greater than 100000 or less than -100000
    are excluded.
    """
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as out_file:
            for line in f:
                if line.startswith('#'):  # Skip header lines
                    continue
                columns = line.strip().split('\t')
                sv_len = ""
                sv_type = ""
                info_parts = []

                # Extract data from the 8th column
                eighth_column = columns[7]
                entries = eighth_column.split(';')
                for entry in entries:
                    if '=' in entry:
                        key, value = entry.split('=')
                        if key == "SVLEN":
                            sv_len = value
                        elif key == "SVTYPE":
                            sv_type = value

                # Check if SVLEN is within the specified range
                if sv_len:
                    sv_len = int(sv_len)
                    if sv_len > 100000 or sv_len < -100000:
                        continue

                # Extract data from columns 10 to the last column
                info_columns = columns[9:]
                for col in info_columns:
                    info_parts.append(col.split(':')[0])

                # Write extracted information to the output file
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(
                    columns[0], columns[1], sv_len, sv_type, '\t'.join(info_parts)))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.vcf output.summary")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_info(input_file, output_file)
    print("Information extracted and written to the file.")

