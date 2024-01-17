#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='Input file')
parser.add_argument('--output', '-o', help='Output file')
parser.add_argument('what_shall_i_write_here', nargs=argparse.REMAINDER)
args = parser.parse_args()

def main():
	distribution = parse_file( args.input )
	print(distribution)
	total_amount = 0
	for i in distribution:
		total_amount += int(i["count"])

	print(total_amount)

	outFile = open(os.path.realpath( args.output ), 'w')
	outFile.write( "%s\t%s\t%s\n" % (
		"length", "count", "percent"
	))
	for i in distribution:
		percent = ( int(i["count"]) * 100 )/total_amount
		outFile.write( "%s\t%s\t%s\n" % (
			i["length"],
			i["count"],
			str(percent)
		) )

	outFile.write( "Total amount of siRNAs: " + str(total_amount) )




def parse_file( file ):
	dis = []
	with open ( file, 'r' ) as bed:
		for line in bed:
			line = line.strip()
			(count, length) = line.split()
			dis.append({"count": count, "length": length})
	return dis

def write_file(wig, fileName):
	print("Writing file:\t" + fileName)

	outFile = open(os.path.realpath( fileName ), 'w')


main()
