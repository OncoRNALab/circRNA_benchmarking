#!/usr/bin/python3

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main script')
parser.add_argument('-i', nargs=1, required=True, help='input fastq dir')
parser.add_argument('-o', nargs=1, required=True, help='output dir')
args = parser.parse_args()
input_dir = args.i[0]
output_dir = args.o[0]

FWD_primer_file = open('/scratch/gent/vo/000/gvo00027/projects/circRNA/MVR1902/exp12_amplicon_seq/FWD_primers.txt')
REV_primer_file = open('/scratch/gent/vo/000/gvo00027/projects/circRNA/MVR1902/exp12_amplicon_seq/REV_primers.txt')


# make FWD primer dict containing the first 16 nts of the FWD primer

FWD_dict = {}

for line in FWD_primer_file:
	FWD_dict[line.rstrip().split('\t')[0][0:16]] =  line.rstrip().split('\t')[0]

FWD_primer_file.close()

# make FWD primmer annotation dict

FWD_an = {}
FWD_primer_file = open('/scratch/gent/vo/000/gvo00027/projects/circRNA/MVR1902/exp12_amplicon_seq/FWD_primers.txt')


for line in FWD_primer_file:
	primer = line.rstrip().split('\t')[0]
	nr = line.rstrip().split('\t')[1]
	FWD_an[primer] = nr

FWD_primer_file.close()

# make REV primer dict containing the first 16 nts of the REV primer

REV_dict = {}

for line in REV_primer_file:
	REV_dict[line.rstrip().split('\t')[0][0:16]] =  line.rstrip().split('\t')[0]

REV_primer_file.close()

# make FWD primmer annotation dict

REV_an = {}
REV_primer_file = open('/scratch/gent/vo/000/gvo00027/projects/circRNA/MVR1902/exp12_amplicon_seq/REV_primers.txt')


for line in REV_primer_file:
	primer = line.rstrip().split('\t')[0]
	nr = line.rstrip().split('\t')[1]
	REV_an[primer] = nr

REV_primer_file.close()

# add some statistics
nr_total = 0
nr_both = 0
nr_R1 = 0
nr_R2 = 0
nr_none = 0
nr_expected = 0

# make list of expected primer combos
primers = open('/scratch/gent/vo/000/gvo00027/projects/circRNA/MVR1902/exp12_amplicon_seq/primers.txt')

primer_combos = []
for combo in primers:
	primer_combos.append(combo.rstrip().replace('\t', '_'))

primers.close()

# iterate through reads and compare with primer seqences

fastq_R1 = open(input_dir + '_R1.fastq')
fastq_R2 = open(input_dir + '_R2.fastq')



for line_R1, line_R2 in zip(fastq_R1, fastq_R2):
	# first fastq line: identifier
	if line_R1[0] == '@':
		new_seq_R1 = line_R1
		new_seq_R2 = line_R2
		line_index = 0

	# second fastq line: sequence => check or primers
	elif line_index == 1:

		## check first 16 nts of R1
		seq_R1 = line_R1
		seq_sub_R1 = seq_R1[0:16]


		FWD_match_R1 = ""
		REV_match_R1 = ""

		if seq_sub_R1 in FWD_dict:
			FWD_match_R1 = FWD_dict[seq_sub_R1]
			new_seq_R1 = new_seq_R1 + seq_R1[len(FWD_match_R1):]
			detected_combo = FWD_match_R1

		elif seq_sub_R1 in REV_dict:
			REV_match_R1 = REV_dict[seq_sub_R1]
			new_seq_R1 = new_seq_R1 + seq_R1[len(REV_match_R1):]
			detected_combo = REV_match_R1


		else:
			new_seq_R1 = new_seq_R1 + line_R1
			detected_combo = 'unknown'

		## check R2
		seq_R2 = line_R2
		seq_sub_R2 = seq_R2[0:16]


		FWD_match_R2 = ""
		REV_match_R2 = ""

		if seq_sub_R2 in FWD_dict:
			FWD_match_R2 = FWD_dict[seq_sub_R2]
			new_seq_R2 = new_seq_R2 + seq_R2[len(FWD_match_R2):]
			detected_combo = FWD_match_R2 + '_' + detected_combo

		elif seq_sub_R2 in REV_dict:
			REV_match_R2 = REV_dict[seq_sub_R2]
			new_seq_R2 = new_seq_R2 + seq_R2[len(REV_match_R2):]
			detected_combo = detected_combo + '_' + REV_match_R2

		else:
			new_seq_R2 = new_seq_R2 + line_R2
			detected_combo = detected_combo + 'unknown'


		# for the stats:
		nr_total += 1
		if (FWD_match_R1 != "" or REV_match_R1 != "") & (FWD_match_R2 != "" or REV_match_R2 != ""):
			nr_both += 1
		elif (FWD_match_R1 != "" or REV_match_R1 != ""):
			nr_R1 += 1
		elif(FWD_match_R2 != "" or REV_match_R2 != ""):
			nr_R2 += 1
		else:
			nr_none += 1



	# third fastq line: +
	elif line_index == 2:
		new_seq_R1 = new_seq_R1 + line_R1
		new_seq_R2 = new_seq_R2 + line_R2

	# fourth fastq line: qual control (trim needed?) + print everything to new fastq file
	elif line_index == 3:

		## R1
		if FWD_match_R1 != "":
			new_seq_R1 = new_seq_R1 + line_R1[len(FWD_dict[seq_sub_R1]):]
			R1_out = 'FWD_' + FWD_an[FWD_match_R1]
			inverted = 'no'

		elif REV_match_R1 != "":
			new_seq_R1 = new_seq_R1 + line_R1[len(REV_dict[seq_sub_R1]):]
			R1_out = "REV_" + REV_an[REV_match_R1]
			inverted = 'yes'

		else:
			new_seq_R1 = new_seq_R1 + line_R1
			R1_out = "unknown"
			inverted = 'no'


		## R2
		if FWD_match_R2 != "":
			new_seq_R2 = new_seq_R2 + line_R2[len(FWD_dict[seq_sub_R2]):]
			R2_out = 'FWD_' + FWD_an[FWD_match_R2]
			inverted = 'yes'

		elif REV_match_R2 != "":
			new_seq_R2 = new_seq_R2 + line_R2[len(REV_dict[seq_sub_R2]):]
			R2_out = "REV_" + REV_an[REV_match_R2]
			inverted = 'no'

		else:
			new_seq_R2 = new_seq_R2 + line_R2
			R2_out = "unknown"
			inverted = 'no'


		if inverted == 'yes':
			tmp = R1_out
			R1_out = R2_out
			R2_out = tmp

		# check if combo was expected
		expected = "not_expected/"
		if detected_combo in primer_combos:
			expected = 'expected/'
			nr_expected += 1


		# paste to files
		fastq_out_R1 = open(output_dir + expected + R1_out + "_" + R2_out + '_R1.fastq', 'a+')
		fastq_out_R2 = open(output_dir + expected + R1_out + "_" + R2_out + '_R2.fastq', 'a+')

		fastq_out_R1.write(new_seq_R1)
		fastq_out_R2.write(new_seq_R2)


	line_index += 1


# print summary

print("total nr of read pairs: " + str(nr_total) + '\nnr of read pairs where both reads start with a primer: ' + str(nr_both) + '\nnr of read pairs where only R1 starts with a primer: ' + str(nr_R1) + '\nnr of read pairs where only R2 starts with a primer: ' + str(nr_R2) + '\nnr of read pairs where both reads do not start with a primer: ' + str(nr_none) + '\nnr of read pairs that match an existing primer combo: ' + str(nr_expected) + '\nnr of read pairs that do not match an existing primer combo: ' + str(nr_total - nr_expected) )
