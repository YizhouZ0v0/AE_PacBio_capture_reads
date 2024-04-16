import gzip
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
from collections import Counter
import subprocess

def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
    reverse_seq = dna_sequence[::-1]
    reverse_complement_seq = ''.join(complement[base] for base in reverse_seq)
    return reverse_complement_seq


source_fastq_gz_path = sys.argv[1]
target_fastq_gz_path = sys.argv[2]


# Define the regular expression to match
first_part = sys.argv[3]
second_part = sys.argv[4]
min_len = sys.argv[5]
max_len = sys.argv[6]
pattern1 = re.compile(r'%s([ATCG]{%s,%s})%s'%(first_part,min_len,max_len,second_part))
pattern2 = re.compile(r'%s([ATCG]{%s,%s})%s'%(reverse_complement(second_part),min_len,max_len,reverse_complement(first_part)))


# Initialize list to store all matching records
matched_records = []
each_read_split_num = []
# Open source file for reading
with gzip.open(source_fastq_gz_path, "rt") as source_handle:
    for record in SeqIO.parse(source_handle, "fastq"):
        # Search for all matching patterns in each sequence
        matches1 = list(pattern1.finditer(str(record.seq)))
        matches2 = list(pattern2.finditer(str(record.seq)))
        
        # Identify positive and negative chains
        if len(matches1) > len(matches2) and len(matches2)==0:
            matches = matches1
        elif len(matches1) < len(matches2) and len(matches1)==0:
            matches = matches2
        elif len(matches1) == 0 and len(matches2) ==0:
            each_read_split_num.append(0)
            continue
        else:
            print("ERROR match {} {}".format(len(matches1),len(matches2)))
            sys.exit(0)
        i = 0
        for match in matches:
            i+=1
            # Extract matching base sequences and corresponding quality scores
            matched_sequence = match.group(1)
            start, end = match.span(1)
            matched_quality = record.letter_annotations["phred_quality"][start:end]
            
            # Create a new SeqRecord object and add a suffix to the ID to distinguish different matching results
            new_record_id = f"{record.id}_match_{i}"
            new_record = SeqRecord(Seq(matched_sequence),
                                   id=new_record_id,
                                   description="{}_{}".format(start, end))
            new_record.letter_annotations["phred_quality"] = matched_quality
            
            # Add a new SeqRecord object to the list
            matched_records.append(new_record)
        each_read_split_num.append(i)
# Open the target file and write all matching records
with gzip.open(target_fastq_gz_path, "wt") as target_handle:
    SeqIO.write(matched_records, target_handle, "fastq")

# Count the number of times each read is divided
print(min(each_read_split_num),max(each_read_split_num))
counts = Counter(each_read_split_num)
print(counts)
