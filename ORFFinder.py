from GeneticCode import *

########################
# Class ORFFinder
########################
class ORFFinder (object):
	def __init__ (self, gencode):
		self.gencode = gencode


	def find (self, dna):
		seqlen = dna.length ()
		starts = []

		seq = str(dna)

		# find start and stop codons for this genetic code...
		start_codons = [ ]
		stop_codons = [ ]
		for i,c in enumerate (str(self.gencode)):
			if c == 'M': start_codons.append (codon_table [i])
			elif c == '*': stop_codons.append (codon_table [i])

		# print "find_ORFs: start_codons =", start_codons
		# print "find_ORFs: stop_codons =", stop_codons

		for i in range (seqlen-2):
			if seq[i:i+3] in start_codons:
				# Found a start codon at position 'i'
				starts.append (i)

		ORFs = [ ]

		for start in starts:
			# l is length from start codon to end of the sequence
			l = seqlen - start

			# i is position pointer in sequence
			i = 0
			# l-2 so no reading past end of the sequence
			dontadd = False
			while i < l-2:
				if seq[start+i:start+i+3] in stop_codons:
					# print "\nORF starts at pos", start, "and ends at pos", start+i+3
					for ORF in ORFs:
						# if we've used this stop codon before, don't include ORF
						if ORF[1] == start+i+3: dontadd = True	
					if (not dontadd): ORFs.append ([start, start+i+3])
					# print_sequence (seq[start:start+i+3], "ORF-"+str(nORFs), gencode)
					starts.remove (start)
					break;
				i += 3

		# for ORF in ORFs:
		# 	print "\nORF starts at pos", ORF[0], "and ends at pos", ORF[1]	
		return ORFs

