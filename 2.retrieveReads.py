#!/usr/bin/env python
import sys
import pysam
import time
import argparse


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument('-b', '--bam', help='input bam file', required=True)
	group.add_argument('-v', '--vcf', help='input vcf file', required=True)
	group.add_argument('-q', '--quality', type=int, help='threshold of mapping quality, default=0', default=0)
	group.add_argument('-o', '--out', help='output file', required=True)

	return group.parse_args()


def read_vcf(in_vcf):
	snp_db = {}
	snp_pos_db = {}
	with open(in_vcf) as f_in:
		for line in f_in:
			if line[0] == '#':
				continue
			data = line.strip().split()
			if data[-1].split(':')[0] != '0/1':
				continue
			chrn = data[0]
			pos = int(data[1])
			ref = data[3]
			alt = data[4]
			if len(ref) != 1 or len(alt) != 1:
				continue
			if chrn not in snp_db:
				snp_db[chrn] = {}
				snp_pos_db[chrn] = []
			snp_pos_db[chrn].append(pos)
			snp_db[chrn][pos] = [ref, alt]
	for chrn in snp_pos_db:
		snp_pos_db[chrn] = sorted(snp_pos_db[chrn])
	return snp_db, snp_pos_db


def bin_search(pos_list, v, pos):
	s = 0
	e = len(pos_list) - 1
	while s<=e:
		mid = (s+e)/2
		if pos_list[mid] > v:
			e = mid - 1
		elif pos_list[mid] < v:
			s = mid + 1
		else:
			if pos == 'r':
				return mid+1
			else:
				return mid
	if pos == 'l':
		return s
	else:
		return e+1


def get_pos_in_range(pos_list, sp, ep):
	sr = bin_search(pos_list, sp, 'l')
	er = bin_search(pos_list, ep, 'r')
	return pos_list[sr: er]


def reverse_base(base):
	if base.lower() == 'a':
		return 't'
	elif base.lower() == 't':
		return 'a'
	elif base.lower() == 'g':
		return 'c'
	elif base.lower() == 'c':
		return 'g'
	else:
		return base


def get_base_pos_with_offset(cigar, query_sequence, offset, is_rev):
	off_ref = 0
	off_query = 0
	for t, n in cigar:
		if t in [0, 1, 4, 7, 8]:	# means M/I/S/=/X, these markers consume query
			off_query += n
		if t in [0, 2, 3, 7, 8]:	# means M/D/N/=/X, these markers consume reference
			off_ref += n
		if off_ref > offset:
			off_ref -= n
			off_query -= n
			break
	off_query += offset-off_ref
	if off_query >= len(query_sequence):
		return ''
	if is_rev:
		return reverse_base(query_sequence[off_query])
	else:
		return query_sequence[off_query]


def retrieve_reads(in_bam, in_vcf, quality, out_list):
	print("\033[32m%s\033[0m Reading VCF"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	snp_db, snp_pos_db = read_vcf(in_vcf)
	reads_db = {}
	for chrn in snp_pos_db:
		reads_db[chrn] = {}
		for pos in snp_pos_db[chrn]:
			reads_db[chrn][pos] = {'ref': {}, 'alt': {}}
	bamfile = pysam.AlignmentFile(in_bam, 'rb')
	
	print("\033[32m%s\033[0m Reading bam"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	for line in bamfile:
		read_name = line.query_name
		flag = bin(line.flag)
		#if line.mapping_quality < quality or line.mapping_quality == 255: # filter data with mapping quality
		#	continue
		if flag[-3] == '1':	# this flag means segment unmapped
			continue
		if len(flag) > 7 and flag[-5] == '1':	# this flag means query seq is reverse complemented
			is_rev = True
		else:
			is_rev = False
		chrn = line.reference_name
		sp = line.reference_start
		ep = sp+line.reference_length
		search_pos = get_pos_in_range(snp_pos_db[chrn], sp, ep)
		cigar = line.cigartuples
		query_sequence = line.query_sequence
		alignment_length = line.query_alignment_length
		mapping_quality = line.mapping_quality
		for i in search_pos:
			offset = i - sp -1
			if offset > len(query_sequence):
				continue
			ref = snp_db[chrn][i][0]
			alt = snp_db[chrn][i][1]
			query_base = get_base_pos_with_offset(cigar, query_sequence, offset, is_rev)
			if query_base.lower() == ref.lower():
				if read_name not in reads_db[chrn][i]['alt']:
					reads_db[chrn][i]['ref'][read_name] = [mapping_quality, alignment_length]
				else:
					amapq, amapl = reads_db[chrn][i]['alt'][read_name]
					if mapping_quality > amapq:
						reads_db[chrn][i]['alt'].pop(read_name)
						reads_db[chrn][i]['ref'][read_name] = [mapping_quality, alignment_length]
					elif alignment_length > amapl:
						reads_db[chrn][i]['alt'].pop(read_name)
						reads_db[chrn][i]['ref'][read_name] = [mapping_quality, alignment_length]
			elif query_base.lower() == alt.lower():
				if read_name not in reads_db[chrn][i]['ref']:
					reads_db[chrn][i]['alt'][read_name] = [mapping_quality, alignment_length]
				else:
					rmapq, rmapl = reads_db[chrn][i]['ref'][read_name]
					if mapping_quality > rmapq:
						reads_db[chrn][i]['ref'].pop(read_name)
						reads_db[chrn][i]['alt'][read_name] = [mapping_quality, alignment_length]
					elif alignment_length > rmapl:
						reads_db[chrn][i]['ref'].pop(read_name)
						reads_db[chrn][i]['alt'][read_name] = [mapping_quality, alignment_length]

	bamfile.close()
	print("\033[32m%s\033[0m Writing result"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	with open(out_list, 'w') as f_out:
		for chrn in sorted(reads_db):
			for pos in sorted(reads_db[chrn]):
				f_out.write("%s,%d,%s,"%(chrn, pos, '|'.join(snp_db[chrn][pos])))
				ref_list = []
				for read_name in sorted(reads_db[chrn][pos]['ref']):
					ref_list.append("%s(%d;%d)"%(read_name, reads_db[chrn][pos]['ref'][read_name][0], reads_db[chrn][pos]['ref'][read_name][1]))
				alt_list = []
				for read_name in sorted(reads_db[chrn][pos]['alt']):
					alt_list.append("%s(%d;%d)"%(read_name, reads_db[chrn][pos]['alt'][read_name][0], reads_db[chrn][pos]['alt'][read_name][1]))
				f_out.write("%s,%s\n"%('|'.join(ref_list), '|'.join(alt_list)))
	print("\033[32m%s\033[0m Finished"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))


if __name__ == "__main__":
	opts = get_opts()
	in_bam = opts.bam
	in_vcf = opts.vcf
	quality = opts.quality
	out_list = opts.out
	retrieve_reads(in_bam, in_vcf, quality, out_list)
