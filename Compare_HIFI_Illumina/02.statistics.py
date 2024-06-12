###raw name: statistics.py
###Usage: python script.py input.txt SVtype_total_number.txt SV_per_chr.txt SV_len.txt

import sys

def analyze_output(input_file, svtype_total_file, sv_per_chr_file, sv_len_file):
    """
    Analyze structural variation (SV) data and generate statistical results.
    
    Parameters:
    input_file (str): Path to the input file containing SV data.
    svtype_total_file (str): Path to the output file to save the total count of each SV type.
    sv_per_chr_file (str): Path to the output file to save the count of each SV type per chromosome.
    sv_len_file (str): Path to the output file to save the length distribution of each SV type.

    Output:
    Generates three files containing the total count of each SV type, the count of each SV type per chromosome, 
    and the length distribution of each SV type.
    """

    svtype_count = {'DEL': 0, 'INS': 0, 'INV': 0, 'DUP': 0}
    chromosome_svtype_count = {}
    svtype_length = {'DEL': [], 'INS': [], 'INV': [], 'DUP': []}

    with open(input_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            chromosome = columns[0]
            svtype = columns[3]
            length_str = columns[2]

            if svtype not in svtype_count:
                continue  # Skip SV types that are not DEL, INS, INV, or DUP

            # Ensure the length is a valid integer
            if not length_str.lstrip('-').isdigit():
                continue

            length = int(length_str)

            # Count the number of SVs per SVTYPE
            svtype_count[svtype] += 1

            # Count the number of SVs per chromosome per SVTYPE
            if chromosome not in chromosome_svtype_count:
                chromosome_svtype_count[chromosome] = {'DEL': 0, 'INS': 0, 'INV': 0, 'DUP': 0}
            chromosome_svtype_count[chromosome][svtype] += 1

            # Collect the lengths of SVs per SVTYPE
            svtype_length[svtype].append(length)

    # Write the results to the respective files
    with open(svtype_total_file, 'w') as out_file:
        out_file.write("Total SVTYPE counts:\n")
        for svtype, count in svtype_count.items():
            out_file.write(f"{svtype}: {count}\n")

    with open(sv_per_chr_file, 'w') as out_file:
        out_file.write("Chromosome SVTYPE counts:\n")
        for chromosome, svtypes in chromosome_svtype_count.items():
            out_file.write(f"{chromosome}:\n")
            for svtype, count in svtypes.items():
                out_file.write(f"  {svtype}: {count}\n")

    with open(sv_len_file, 'w') as out_file:
        out_file.write("SVTYPE length distribution:\n")
        for svtype, lengths in svtype_length.items():
            if lengths:
                avg_length = sum(lengths) / len(lengths)
                out_file.write(f"{svtype}: average length = {avg_length:.2f}, max length = {max(lengths)}, min length = {min(lengths)}\n")
            else:
                out_file.write(f"{svtype}: No data\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py input.txt SVtype_total_number.txt SV_per_chr.txt SV_len.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    svtype_total_file = sys.argv[2]
    sv_per_chr_file = sys.argv[3]
    sv_len_file = sys.argv[4]

    analyze_output(input_file, svtype_total_file, sv_per_chr_file, sv_len_file)
    print("Information analyzed and written to the files.")



      
