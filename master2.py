import argparse
import sys
import os
from degenerate import degen
from translation import translate
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbitblastnCommandline

#Command Line Arguments
parser = argparse.ArgumentParser(description = "Get translations from a genbank file")
parser.add_argument('-g', '--genome', help = "Path to genbank file", required = False, default = "sbw25.gbff")
parser.add_argument('-t', '--tag', help = "Gene tag", required = False, default = ["PFLU_RS02030"], nargs = '*')
parser.add_argument('--use_old', action = 'store_true', help = 'Use old locus tag entry')
args = parser.parse_args()

#New/Old Tag
tag_type = 'locus_tag'
if args.use_old:
	tag_type = 'old_locus_tag'

#Translation/DNA Dictionaries
locusToTranslation = dict()
locusToDNA = dict()
def get_cds_feature_with_qualifier_value(genbank_path, tag_type, tag):
	for feature in genbank_path.features:
		if feature.type == "CDS":
			if feature.qualifiers.get(tag_type) is not None and feature.qualifiers.get('translation') is not None:
				locusToTranslation.update({feature.qualifiers.get(tag_type)[0] : feature.qualifiers.get('translation')[0]})
				gene_sequence = feature.extract(genome_record.seq)
				locusToDNA[feature.qualifiers.get(tag_type)[0]] = gene_sequence

if __name__ == "__main__":
	for tag in args.tag:
		genome_record = SeqIO.read(args.genome, "genbank")
		cds_feature = get_cds_feature_with_qualifier_value(genome_record, tag_type, tag)

#Mapping Keys/Values & Creating FASTA Files
def write_translations(tagInterest):
	interest = open('translationInterest.fasta', 'w')
	notInterest = open('translationOthers.fasta', 'w')
	degenInterest = open('degenInterest.fasta', 'w')
	translationTag = locusToTranslation.keys()
	for tag in translationTag:
		degenTag = degen(locusToTranslation[tag])
		if tag == tagInterest:
			interest.write('>' + tag + '\n' + locusToTranslation[tag] + '\n')
			degenInterest.write('>' + tag + '\n' + degenTag + '\n')
		else:
			notInterest.write('>' + tag + '\n' + locusToTranslation[tag] + '\n')
write_translations(tag)

def write_DNA(tagInterest):
	interest = open('dnaInterest.fasta', 'w')
	notInterest = open('dnaOthers.fasta', 'w')
	dnaTag = locusToDNA.keys()
	for key in dnaTag:
		if key == tagInterest:
			interest.write('>' + str(key) + '\n' + str(locusToDNA[key]) + '\n')
		else:
			notInterest.write('>' + str(key) + '\n' + str(locusToDNA[key]) + '\n')
write_DNA(tag)

RF = open('RF.fasta', 'w')
def read_frame(seq, locus_tag, frame = 1):
	if frame >= 0:
		seq = seq[frame:]
	else:
		seq = seq.reverse_complement()[abs(frame):]

	trim_char = len(seq)%3
	if trim_char != 0:
		seq = seq[:-trim_char]

	RF.write('>' + locus_tag + str(frame) + '\n' + str(seq) + '\n')

def write_RF(tagInterest):
	# RF = open('RF.fasta', 'w')
	rfTag = locusToDNA.keys()
	for key in rfTag:
		if key == tagInterest:
			read_frame(locusToDNA[key], key, frame = 1)
			read_frame(locusToDNA[key], key, frame = 2)
			read_frame(locusToDNA[key], key, frame = -0)
			read_frame(locusToDNA[key], key, frame = -1)
			read_frame(locusToDNA[key], key, frame = -2)

			# RF1 = locusToDNA[key][1:]
			# trim_char1 = len(RF1)%3
			# if trim_char1 != 0:
			# 	trim_char1 = RF1[:-trim_char1]
			# 	RF.write('>' + str(key) + '+1' + '\n' + translate(trim_char1) + '\n')
			# else:
			# 	RF.write('>' + str(key) + '+1' + '\n' + translate(RF1) + '\n')
			# RF2 = locusToDNA[key][2:]
			# trim_char2 = len(RF2)%3
			# if trim_char2 != 0:
			# 	trim_char2 = RF2[:-trim_char2]
			# 	RF.write('>' + str(key) + '+2' + '\n' + translate(trim_char2) + '\n')
			# else:
			# 	RF.write('>' + str(key) + '+2' + '\n' + translate(RF2) + '\n')
			# #Opposite Strand
			# rev = locusToDNA[key].reverse_complement()
			
			# RF0 = rev[0:]
			# trim_char0 = len(RF0)%3
			# if trim_char0 != 0:
			# 	trim_char0 = RF0[:-trim_char0]
			# 	RF.write('>' + str(key) + '-0' + '\n' + translate(trim_char0) + '\n')
			# else:
			# 	RF.write('>' + str(key) + '-0' + '\n' + translate(RF0) + '\n')
			# RF4 = rev[1:]
			# trim_char4 = len(RF4)%3
			# if trim_char4 != 0:
			# 	trim_char4 = RF4[:-trim_char4]
			# 	RF.write('>' + str(key) + '-1' + '\n' + translate(trim_char4) + '\n')
			# else:
			# 	RF.write('>' + str(key) + '-1' + '\n' + translate(RF4) + '\n')
			# RF5 = rev[2:]
			# trim_char5 = len(RF5)%3
			# if trim_char5 != 0:
			# 	trim_char5 = RF5[:-trim_char5]
			# 	RF.write('>' + str(key) + '-2' + '\n' + translate(trim_char5) + '\n')
			# else:
			# 	RF.write('>' + str(key) + '-2' + '\n' + translate(RF5) + '\n')

write_RF(tag)
print('Output in translationInterest.fasta, translationOthers.fasta, degenInterest.fasta, dnaInterest.fasta, dnaOthers.fasta, RF.fasta, blastp.out, tblastn.out')

#Local Blast
query = os.path.join('RF.fasta')
db = os.path.join('blastdb2')
blastout = os.path.join('blastp.out')

cmd_blastp = NcbiblastpCommandline(
    query=query, out=blastout, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', db=db)
cmd_blastp()

query2 = os.path.join('translationInterest.fasta')
db2 = os.path.join('blastdb')
blastout2 = os.path.join('tblastn.out')

cmd_tblastn = NcbitblastnCommandline(
    query=query2, out=blastout2, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore frames', db=db2)
cmd_tblastn()
