from datetime import datetime

###############################
# Class BioPolymer
###############################
class BioPolymer (object):
	seq = ""
	name = ""

	# placeholder values...
	filter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*"
	masses = { 'A' : 0.0, 'B' : 0.0, 'C' : 0.0, 'D': 0.0, 'E' : 0.0, 'F' : 0.0,
	'G' : 0.0, 'H' : 0.0, 'I' : 0.0, 'J' : 0.0, 'K' : 0.0, 'L' : 0.0, 'M' : 0.0,
	'N' : 0.0, 'O' : 0.0, 'P' : 0.0, 'Q' : 0.0, 'R' : 0.0, 'S' : 0.0, 'T' : 0.0, 
	'U' : 0.0, 'V' : 0.0, 'W' : 0.0, 'X' : 0.0, 'Y' : 0.0, 'Z' : 0.0, '*' : 0.0 }
	
	def __init__ (self, seq, name = ("unnamed-" + str(datetime.now()))):
		self.name = name
	
		s = seq.upper ()
		newseq = ""
		for c in s:
			if c in self.filter:
				newseq += c
			else:
				print "rejecting char:", c
		self.seq = newseq
		print "New", type(self), "created:", self.name, "; length =", len(self.seq)
		
	def seqname (self):
		return self.name

	def length (self):
		return len(self.seq)
		
	def molwt (self):
		mw = 0.0
		for c in self.seq:
			mw += self.masses [c]
		return mw


	def __repr__ (self):
		return self.seq


###############################
# Class DNA
###############################
class DNA (BioPolymer):
	filter = "ACGT"
	masses = { 'A' : 313.2, 'C' : 289.2, 'G' : 329.2, 'T' : 304.2 }

	def transcribe (self):
		rna_seq = ""
		for c in self.seq:
			if c == 'T' : rna_seq += 'U'
			else: rna_seq += c 
		name = self.name + ".rna"
		rna = RNA (rna_seq, name)
		return rna 

	# reverse_complement() returns the reverse complement of seq in 5'-3' direction 
	def reverse_complement (self):
		s = self.seq[::-1]
		rseq = ""

		for c in s:
			if c == 'A':
				rseq += 'T'
			elif c == 'C':
				rseq += 'G'
			elif c == 'G':
				rseq += 'C'
			elif c == 'T':
				rseq += 'A'
			else:
				print "Error: reverse_complement: invalid character:", c
		name = self.name + ".reverse_complement"
		return DNA (rseq, name)


###############################
# Class RNA
###############################
class RNA (BioPolymer):
	filter = "ACGU"
	masses = { 'A' : 329.2, 'C' : 305.2, 'G' : 345.2, 'U' : 306.2 }
	
###############################
# Class Protein
###############################
class Protein (BioPolymer):
	filter = "ACDEFGHIKLMNPQRSTUVWY*"

	### masses from http://web.expasy.org/findmod/findmod_masses.html#AA
	### use 18 for stop codon to acct for addition of one water mol
	masses = { 'A' : 71.0788, 'C' : 103.1388, 'D' : 115.0886, 'E' : 129.1155, 
	'F' : 147.1766, 'G' : 57.0519, 'H' : 137.1411, 'I' : 113.1594, 'K' : 128.1741, 
	'L' : 113.1594, 'M' : 131.1926, 'N' : 114.1038, 'O' : 237.3018, 'P' : 97.1167, 
	'Q' : 128.1307, 'R' : 156.1875, 'S' : 87.0782, 'T' : 101.1051, 'U' : 150.0388, 
	'V' : 99.1326, 'W' : 186.2132, 'Y' : 163.1760, '*' : 18.0 }

	## subtract one to account for stop codon
	def length (self):
		return len(self.seq) - 1
		

	