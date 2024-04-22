#this python script aims to extract one assembly unique insertion region compared to other assemblies which input files were from syri outputs
#due to the large assembly, i split the chromosomes to one chromosome per file, and aligned all the assemblies to Chinese Spring v2.1 reference using minimap2 with the parameters "-ax asm5 -I100g -f100 -t 128 --eqx", output must be sam becaose of the length (see command_minimap2.sh).
#Then using syri to get strcture variations with -k -F S --nc 1 (see command_syri.sh)
#here i wanna compare wheat variety 3504 with 13 wheat assemblies
#list.txt 
#Arinalrfor.1A/Arinalrfor.1A_CS21syri.out
#CS21.1A/CS21.1A_CS21syri.out
#Fielder.1A/Fielder.1A_CS21syri.out
#Jagger.1A/Jagger.1A_CS21syri.out
#Julius.1A/Julius.1A_CS21syri.out
#Kariega.1A/Kariega.1A_CS21syri.out
#Lancer.1A/Lancer.1A_CS21syri.out
#Landmark.1A/Landmark.1A_CS21syri.out
#Mace.1A/Mace.1A_CS21syri.out
#Mattis.1A/Mattis.1A_CS21syri.out
#Norin61.1A/Norin61.1A_CS21syri.out
#Renan.1A/Renan.1A_CS21syri.out
#Stanley.1A/Stanley.1A_CS21syri.out
#usage: python unique_insertion.py -F 3504_v1.1A/3504_v1.1A_CS21syri.out -L list.txt -O output.txt
import argparse

def extract_ins_records(file_list):
    ins_records = set()
    for file_path in file_list:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) >= 11 and columns[10] == 'INS':
                    ins_records.add((columns[0], columns[1]))
    return sorted(list(ins_records))

def extract_unique_records(file_path, unique_records):
    unique_lines = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 11 and columns[10] == 'INS' and (columns[0], columns[1]) in unique_records:
                unique_lines.append(line)
    return unique_lines

def main():
    parser = argparse.ArgumentParser(description='Extract unique INS records from file1.txt compared to other files.')
    parser.add_argument('-F', '--firstfile', required=True, help='Path to the first file (file1.txt)')
    parser.add_argument('-O', '--output', required=True, help='Output file path')
    parser.add_argument('-L', '--listfile', required=True, help='Path to the file containing list of other files')
    args = parser.parse_args()

    # 从文件列表中读取其他文件的路径
    with open(args.listfile, 'r') as listfile:
        other_files = [line.strip() for line in listfile]

    # 提取所有文件中的 INS 记录的第一列和第二列
    all_ins_records = extract_ins_records(other_files)

    # 提取第一个文件的 INS 记录的第一列和第二列
    firstfile_ins_records = set(extract_ins_records([args.firstfile]))

    # 找出只有第一个文件有而其他文件没有的 INS 记录
    unique_records = firstfile_ins_records.difference(all_ins_records)

    # 提取第一个文件中独特的记录
    unique_lines = extract_unique_records(args.firstfile, unique_records)

    # 将结果写入输出文件
    with open(args.output, 'w') as output_file:
        for line in unique_lines:
            output_file.write(line)

if __name__ == "__main__":
    main()
