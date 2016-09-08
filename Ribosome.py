from GeneticCode import three_letter_code

###############################
# helper functions
###############################
def _lookup (nt):
	if nt == 'A':
		return 0
	elif nt == 'C':
		return 1
	elif nt == 'G':
		return 2
	return 3	

###############################
# Class Ribosome
###############################
class Ribosome (object):

	def __init__ (self, genetic_code):
		self.genetic_code = genetic_code

	def translate_codon (self, codon, mode = 0):
		### assumes a previously cleaned up sequence (upper case, no invalid chars)
		### if (mode) returns 3 letter translation; otherwise, single letter
	
		translation = ""
		# gencode = self.genetic_code.code()
		gencode = str(self.genetic_code)

		if len(codon) != 3:
			print "Error: translate_codon: not a valid codon:", codon
		else:
			first = _lookup (codon[0])
			second = _lookup (codon[1])
			third = _lookup (codon[2])
			translation = gencode [(first * 16) + (second * 4) + third]
			if (mode):
				translation = three_letter_code[translation]
		return translation
	
	def translate_seq (self, seq, mode = 0):
		### assumes a previously cleaned up sequence
		### if (mode) returns 3 letter translation; otherwise, single letter
		
		translation = [ ]
		codons = [ ]

		ncodons = len(seq) / 3		
	
# 		seq = seq[offset:]
# 			if (debug): print "translate_seq: working seq = ", seq
		for i in range(ncodons):			
			codon = seq[(i*3)]+seq[(i*3)+1]+seq[(i*3)+2]
#				if (debug) : print "translate_seq: found codon",codon
			codons.append (codon) 
		for i in codons:
			codon = self.translate_codon (i, mode)
			translation.append (codon)
		translation = ''.join(translation)
#			if (debug): print "translate_seq: translation =", translation
		return translation
