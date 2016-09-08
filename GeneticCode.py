###############################
# Class GeneticCode
#
# AAA, AAC, AAG, AAU, ACA, ACC, ACG, ACU, AGA ...
# A = 0, C = 1, G = 2, T/U = 3
# Index = (16 * first) + (4 * second) + third 
##############################
std_genetic_code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CYCLFLF"

three_letter_code = { 'A' : 'Ala', 'B' : '???', 'C' : 'Cys', 'D' : 'Asp', 'E' : 'Glu', 'F' : 'Phe', \
'G' : 'Gly', 'H' : 'His', 'I' : 'Ile', 'J' : '???', 'K' : 'Lys', 'L' : 'Leu', 'M' : 'Met', 'N' : 'Asn', \
'O' : 'Pyr', 'P' : 'Pro', 'Q' : 'Gln', 'R' : 'Arg', 'S' : 'Ser', 'T' : 'Thr', 'U' : 'Sel', 'V' : 'Val', 'W' : 'Trp', \
'X' : '???', 'Y' : 'Tyr', 'Z' : '???', '*' : '***' }

codon_table = [ 'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG',
'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 
'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 
'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 
'TGT', 'TTA', 'TTC', 'TTG', 'TTT' ]

class GeneticCode (object):
	def __init__ (self, gencode = std_genetic_code, name = "Standard Genetic Code"):
		self.name = name
		self.gencode = gencode
		print "Using genetic code %s" % (self.name)
		
	def __repr__ (self):
		return self.gencode
		
	def code(self):
		for i,c in enumerate (self.gencode):
			print codon_table[i], c, three_letter_code [c]
		# return self.gencode

	# given a codon (eg, 'AAA', returns translation, 'Lys')
	def get_codon (self, codon):
		if len(codon) != 3:
			print "Error: get_codon: not a valid codon:", codon
		else:
			first = _lookup (codon[0])
			second = _lookup (codon [1])
			third = _lookup (codon [2])
			aa = self.gencode [(first * 16) + (second *4) + third]
			return three_letter_code [aa]

	def set_codon (self, codon, assignment, three_letter):
		if len(codon) != 3:
			print "Error: translate_codon: not a valid codon:", codon
		if not assignment.isalpha() or len(assignment) != 1:
			print "Error: invalid input:", assignment
		else:
			first = _lookup (codon[0])
			second = _lookup (codon[1])
			third = _lookup (codon[2])
			g = list (self.gencode)
			g[(first * 16) + (second * 4) + third] = assignment
			self.gencode = ''.join(g)
			self.three_letter_code[assignment] = three_letter
			print "Codon %s now defined as %s, %s" % (codon, assignment, three_letter)
			print self.gencode


