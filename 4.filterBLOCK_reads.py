#!/usr/bin/env python
import sys
import argparse


def GetOpts():
	group = argparse.ArgumentParser()
	group.add_argument('-i', '--input', help='input block file', required=True)
	group.add_argument('-o', '--out', help='output file prefix', required=True)

	return group.parse_args()


def FilterBlock(inFile, outPrefix):
	reads_db = {}
	with open(inFile, 'r') as fin:
		for line in fin:
			if line[:2] == "ID":
				continue
			data = line.strip().split(',')
			x = int(data[5].split('|')[0]) + 7
			y = int(data[6].split('|')[0]) + 7
			for read in data[x].split('|'):
				if read not in reads_db:
					reads_db[read] = {'x': 0, 'y': 0}
				reads_db[read]['x'] += 1
			for read in data[y].split('|'):
				if read not in reads_db:
					reads_db[read] = {'x': 0, 'y': 0}
				reads_db[read]['y'] += 1
	
	fx = open(outPrefix+"-x.list", 'w')
	fy = open(outPrefix+"-y.list", 'w')
	for read in reads_db:
		if read == '':
			continue
		if reads_db[read]['x'] > reads_db[read]['y']:
			fx.write(read+"\n")
		elif reads_db[read]['x'] < reads_db[read]['y']:
			fy.write(read+"\n")
	fx.close()
	fy.close()


if __name__ == "__main__":
	opts = GetOpts()
	inFile = opts.input
	outPrefix = opts.out
	FilterBlock(inFile, outPrefix)
