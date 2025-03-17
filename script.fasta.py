import argparse #导入模块，解析命令行参数
import gzip #压缩输出
import re #如理正则表达式

from Bio import SeqIO #读写fasta文件
from Bio.Data import CodonTable
from Bio.Seq import Seq #表示生物序列
from Bio.SeqRecord import SeqRecord #seqIO读写以后的东西，包括序列和源数据

def read_id_list(file_path):#获得一个只有id的集合
    """读取 ID 列表文件"""
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f)

def filter_by_id(records, target_id):
    """根据目标 ID 过滤序列"""
    for record in records:
        if record.id == target_id:
            yield record

def exclude_by_name(records, excluded_name):
    """根据排除名称过滤序列"""
    for record in records:
        if excluded_name not in record.id:
            yield record

def filter_by_regex(records, regex):
    """根据正则表达式过滤序列"""
    pattern = re.compile(regex)
    for record in records:
        if pattern.search(record.id):
            yield record

def filter_by_length(records, min_length=None, max_length=None):
    """根据长度过滤序列"""
    for record in records:
        seq_length = len(record.seq)
        if (min_length is None or seq_length >= min_length) and \
           (max_length is None or seq_length <= max_length):
            yield record

def add_prefix_to_id(records, prefix, separator, position):
    """为 ID 添加前缀"""
    for record in records:
        parts = record.id.split(separator)
        parts.insert(position, prefix)
        record.id = separator.join(parts)
        yield record

def aligned_to_unaligned(records):
    """将对齐序列转换为未对齐格式"""
    for record in records:
        record.seq = record.seq.ungap('-')
        yield record

def multi_to_single_line(records):
    """将多行序列转换为单行格式"""
    for record in records:
        record.seq = Seq(str(record.seq).replace('\n', ''))
        yield record

def translate_to_protein(records, codon_table):
    """将 CDS 序列翻译为蛋白质"""
    table = CodonTable.ambiguous_dna_by_name[codon_table]
    for record in records:
        protein_seq = record.seq.translate(table=table)
        yield SeqRecord(protein_seq, id=record.id, description="")

def write_output(records, output_file, compress=False):
    """将序列写入输出文件"""
    if compress:
        with gzip.open(output_file, 'wt') as f:
            SeqIO.write(records, f, 'fasta')
    else:
        with open(output_file, 'w') as f:
            SeqIO.write(records, f, 'fasta')

def basic_statistics(records):
    """计算基本统计信息"""
    num_sequences = 0
    total_length = 0
    for record in records:
        num_sequences += 1
        total_length += len(record.seq)
    print(f"Number of sequences: {num_sequences}")
    print(f"Total length: {total_length}")
    print(f"Average length: {total_length / num_sequences if num_sequences > 0 else 0}")
def cut_by_file_number(records, num_files, output_prefix):
    """将序列按文件数切割"""
    records = list(records)  # 将生成器转换为列表
    chunk_size = len(records) // num_files
    for i in range(num_files):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < num_files - 1 else len(records)
        output_file = f"{output_prefix}_{i + 1}.fasta"
        write_output(records[start:end], output_file)

def cut_by_sequence_number(records, num_sequences, output_prefix):
    """将序列按序列数切割"""
    records = list(records)  # 将生成器转换为列表
    num_files = len(records) // num_sequences
    for i in range(num_files):
        start = i * num_sequences
        end = (i + 1) * num_sequences if i < num_files - 1 else len(records)
        output_file = f"{output_prefix}_{i + 1}.fasta"
        write_output(records[start:end], output_file)

#设置参数，即快捷键
def main():
    parser = argparse.ArgumentParser(description="Process FASTA/FASTQ files with various options.")
    parser.add_argument('-n', '--target_id', type=str, help="Target ID exactly")
    parser.add_argument('-ex', '--excluded_name', type=str, help="Excluded name")
    parser.add_argument('-rg', '--regex', type=str, help="Target regex word to find sequences")
    parser.add_argument('-nl', '--id_list_file', type=str, help="Select sequences by an ID list")
    parser.add_argument('-exl', '--exclude_list_file', type=str, help="Select sequences by an excluding list")
    parser.add_argument('-rgl', '--regex_list_file', type=str, help="Select sequences by a regular expression list")
    parser.add_argument('-ap', '--add_prefix', type=str, help="Add a specific prefix to ID")
    parser.add_argument('-sep', '--separator', type=str, default='|', help="Separator for adding a specific prefix to ID")
    parser.add_argument('-pOS', '--position', type=int, default=0, help="Position of tag to add")
    parser.add_argument('-lgt', '--min_length', type=int, help="Select sequences larger than [INT] bp")
    parser.add_argument('-lle', '--max_length', type=int, help="Select sequences smaller than [INT] bp")
    parser.add_argument('-a2u', '--aligned_to_unaligned', action='store_true', help="Change format, aligned sequence to unaligned")
    parser.add_argument('-rW', '--multi_to_single_line', action='store_true', help="Change format, multi-line to single line")
    parser.add_argument('-tr', '--translate', action='store_true', help="Translate FASTA into protein (CDS only)")
    parser.add_argument('-codon', '--codon_table', type=str, default='Standard', help="Codon table usage (CDS only, with -tr)")
    parser.add_argument('-cut', '--cut_by_file_number', type=int, help="Cut FASTA into files by file number")
    parser.add_argument('-cuts', '--cut_by_sequence_number', type=int, help="Cut FASTA into files by specific sequence number")
    parser.add_argument('-o', '--output', type=str, help="Output file or outdoor for -cut/-cuts")
    parser.add_argument('-z', '--compress', action='store_true', help="Output a gzip type file")
    parser.add_argument('-sta', '--statistics', action='store_true', help="Basic statistics of FASTA")
    parser.add_argument('input_file', type=str, help="Input FASTA/FASTQ file")

    args = parser.parse_args()#把设置的参数整合到args这个盒子中

    # 读取输入文件
    records = SeqIO.parse(args.input_file, 'fasta')#引入待处理文件
#使用参数调用函数
    # 应用过滤器
    if args.target_id:
        records = filter_by_id(records, args.target_id)
    if args.excluded_name:
        records = exclude_by_name(records, args.excluded_name)
    if args.regex:
        records = filter_by_regex(records, args.regex)
    if args.id_list_file:
        id_list = read_id_list(args.id_list_file)
        records = (record for record in records if record.id in id_list)
    if args.exclude_list_file:
        exclude_list = read_id_list(args.exclude_list_file)
        records = (record for record in records if record.id not in exclude_list)
    if args.regex_list_file:
        with open(args.regex_list_file, 'r') as f:
            regex_list = [line.strip() for line in f]
        records = (record for record in records if any(re.search(regex, record.id) for regex in regex_list))
    if args.min_length or args.max_length:
        records = filter_by_length(records, args.min_length, args.max_length)
    if args.add_prefix:
        records = add_prefix_to_id(records, args.add_prefix, args.separator, args.position)
    if args.aligned_to_unaligned:
        records = aligned_to_unaligned(records)
    if args.multi_to_single_line:
        records = multi_to_single_line(records)
    if args.translate:
        records = translate_to_protein(records, args.codon_table)
    if args.cut_by_file_number:
        cut_by_file_number(records, args.cut_by_file_number, args.output)
    if args.cut_by_sequence_number:
        cut_by_sequence_number(records, args.cut_by_sequence_number, args.output)
    # 输出结果
    if args.output:
        write_output(records, args.output, args.compress)
    if args.statistics:
        basic_statistics(records)

#加一个开关，控制main函数的使用，如果是导入这个脚本就不会使用这个函数，如果直接运行这个函数会正常使用
if __name__ == "__main__":
    main()
#1.先导入需要使用的模块和包
#2.根据所需功能设计函数
#3.argparse模块可以用来设置快捷键即参数
#4.用parser.add_argument的方法设置参数，短名称即为快捷键，长名称联系上功能函数的形参
#5.最后用if语句当用户使用arg.指定参数时，调用相关的函数