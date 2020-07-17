#!/usr/bin/env python3
#Emma Kirkegaard (ekirkega)

'''
Run by specifying an input file, output file, input type (-P33in, -P64in, -P64Bin, or -P64SOLin)
and ouput type  (P33out or P64out)
Example:
		$python3 fastqConverter.py <infile.fastq >outfile.fastq -P33in -P64out
'''

import sys
import math

class FastQreader :
	''' 
	Define objects to read FastQ files. 
	'''
	def __init__ (self, fname=''):
		'''contructor: saves attribute fname '''
		self.fname = fname
			
	def doOpen (self):
		''' Handle file opens, allowing STDIN.'''
		if self.fname is '':
			return sys.stdin
		else:
			return open(self.fname)
		
	def readFastQ (self):
		''' 
		Read an entire FastQ record and return the header sequence 
		quality header and quality score sequence
		'''
		header = ''
		sequence = ''
		qHeader = ''
		qScore = ''
		
		with self.doOpen() as fileH:
			
			header = ''
			sequence = ''
			qHeader = ''
			qScore = ''
			
			# skip to first fasta header
			line = fileH.readline()

			while not line.startswith('@') :
				line = fileH.readline()

			header = line[1:].rstrip()

			inSeq = True #are we in the ID/sequence portion? 
			inQscore = False #are we in the Qscore portion?

			for line in fileH:
				if line.startswith ('@'): #new sequnce is starting 
					inSeq = True 
					inQscore = False
					yield header,sequence, qHeader, qScore #return prev sequence
					header = line.rstrip()
					sequence = ''
				elif line.startswith ('+'): #qscore portion starting 
					inQscore = True
					inSeq = False
					qHeader = line.rstrip()
					qScore = ''
				elif inSeq:
					sequence += ''.join(line.rstrip().split()).upper().replace('*','N').replace('.','N')
				elif inQscore:
					qScore += ''.join(line.rstrip().split())

		yield header,sequence, qHeader, qScore


class FastQ :
	'''FASTQ document object'''

	def __init__(self, qScore):
		'''initialize to hold a string of qscore and the file type'''
		self.qScore = qScore

		#dictionaries for conversion 
		self.p33top64 = {chr(q):chr(q+31) for q in range(33, 73 + 1)}
		self.p64top33 = {chr(q):chr(q-31) for q in range(64, 104 + 1)}
		self.Soltop64 = {';':'A', '<':'A', '=':'B', '>':'B', '?':'C', '@':'C', 'A':'D', 'B':'D',
						 'C':'E', 'D':'E', 'E':'F', 'F':'G', 'G':'H', 'H':'I', 'I':'J', 'J':'J'
		 				}


########################################################################
# CommandLine
########################################################################
class CommandLine() :
	'''
	Handle the command line, usage and help requests.

	CommandLine uses argparse, now standard in 2.7 and beyond. 
	it implements a standard command line argument parser with various argument options,
	a standard usage and help.

	attributes:
	all arguments received from the commandline using .add_argument will be
	avalable within the .args attribute of object instantiated from CommandLine.
	For example, if myCommandLine is an object of the class, and requiredbool was
	set as an option using add_argument, then myCommandLine.args.requiredbool will
	name that option.
 
	'''


	
	def __init__(self, inOpts=None) :
		'''
		Implement a parser to interpret the command line argv string using argparse.
		'''
		
		import argparse
		self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
											 epilog = 'Program epilog - some other stuff you feel compelled to say', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s [options] -option1[default] <input >output'
											 )
		self.parser.add_argument('-P33in', action = 'store_true', help='specify input file type')
		self.parser.add_argument('-P64in', action = 'store_true', help='specify input file type') 
		self.parser.add_argument('-P64Bin', action = 'store_true', help='specify input file type')
		self.parser.add_argument('-P64SOLin', action = 'store_true', help='specify input file type')
		self.parser.add_argument('-P33out', action = 'store_true', help='specify output file type') 
		self.parser.add_argument('-P64out', action = 'store_true', help='specify output file type') 
		if inOpts is None :
			self.args = self.parser.parse_args()
		else :
			self.args = self.parser.parse_args(inOpts)

########################################################################
#helper function
########################################################################

def convert(info,func):
	for h, s, qh, qs in info:
		print(h)
		print(s)
		print(qh)
		score = FastQ(qs)
		func(score)

def convert33to64(self):
	scores = ''
	for q in self.qScore:
		scores+=self.p33top64.get(q)
	print(scores)

def convert64to33(self):
	scores = ''
	for q in self.qScore:
		scores+=self.p64top33.get(q)
	print(scores)

def convertSoltop64(self):
	scores = ''
	for q in self.qScore:
			if ord(q) < 74:
				q = self.Soltop64.get(q)
			scores += q
	print(scores)

def convertSoltop33(self):
	scores = ''
	for q in self.qScore:
		if ord(q) < 74:
			q = self.Soltop64.get(q)
		q = self.p64top33.get(q)
		scores += q
	print(scores)

def convert64Bto64(self):
	scores = ''
	for q in self.qScore:
		if q == 'B':
			q = '@'
		scores+=q
	print(scores)

def convert64Bto33(self):
	scores = ''
	for q in self.qScore:
		if q == 'B':
			q = '@'
		q = self.p64top33(q)
		scores+=q
	print(scores)
########################################################################
#main
########################################################################
def main(inCL=None):
	'''print out the resulting converted file'''
	if inCL is None:
		myCommandLine = CommandLine()
	else :
		myCommandLine = CommandLine(inCL)

	options = myCommandLine.args
	myReader = FastQreader()
	info = myReader.readFastQ()
	if options.P64out:
		if options.P33in:
			convert(info, convert33to64)
		elif options.P64Bin:
			convert(info, convert64Bto64)
		elif options.P64SOLin:
			convert(info, convertSoltop64)
		else: #return P64
			for h,s,qh,qs in info:
				print(h)
				print(s)
				print(qh)
				print(qs)
	#P33out
	elif options.P33out:
		if options.P64in:
			convert(info,convert64to33)
		elif options.P64SOLin:
			convert(info,convertSoltop33)
		elif options.P64Bin:
			convert(info,convert64Bto33)
		else: #return P33
			for h,s,qh,qs in info:
				print(h)
				print(s)
				print(qh)
				print(qs)
		

if __name__ == "__main__":
	main()




