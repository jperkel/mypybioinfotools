## bioinformatics tools v2

std_genetic_code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CYCLFLF"

three_letter_code = { 'A' : 'Ala', 'B' : '???', 'C' : 'Cys', 'D' : 'Asp', 'E' : 'Glu', 'F' : 'Phe', \
'G' : 'Gly', 'H' : 'His', 'I' : 'Ile', 'J' : '???', 'K' : 'Lys', 'L' : 'Leu', 'M' : 'Met', 'N' : 'Asn', \
'O' : '???', 'P' : 'Pro', 'Q' : 'Gln', 'R' : 'Arg', 'S' : 'Ser', 'T' : 'Thr', 'U' : '???', 'V' : 'Val', 'W' : 'Trp', \
'X' : '???', 'Y' : 'Tyr', 'Z' : '???', '*' : '***' }


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
# Class BioPolymer
###############################
class BioPolymer (object):
	seq = ""
	name = ""
	filter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*"
#	molwt = [ ]
	
	def __init__ (self, name, seq):
		self.name = name
	
		s = seq.upper ()
		newseq = ""
		for c in s:
			if c in self.filter:
				newseq += c
			else:
				print "rejecting char:", c
		self.seq = newseq
		print "New biopolymer created: ", self.name, self.seq
		
	def seqname (self):
		return self.name

# 	def sequence (self):
# 		return self.seq
		
	def length (self):
		return len(self.seq)
		
	def __repr__ (self):
		return self.seq


###############################
# Class DNA
###############################
class DNA (BioPolymer):
	filter = "ACGT"

	
###############################
# Class RNA
###############################
class RNA (BioPolymer):
	filter = "ACGU"

	
###############################
# Class Protein
###############################
class Protein (BioPolymer):
	filter = "ACDEFGHIKLMNPQRSTVWY*"

	
###############################
# Class GeneticCode
#
# AAA, AAC, AAG, AAU, ACA, ACC, ACG, ACU, AGA ...
# A = 0, C = 1, G = 2, T/U = 3
# Index = (16 * first) + (4 * second) + third 
##############################
class GeneticCode (object):
	def __init__ (self, gencode = std_genetic_code, name = "Standard Genetic Code"):
		self.name = name
		self.gencode = gencode
		print "Using genetic code %s" % (self.name)
		
	def __repr__ (self):
		return "(%s: %s)" % (self.name, self.gencode)
		
	def code(self):
		return self.gencode
		
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


###############################
# Class Ribosome
###############################
class Ribosome (object):

	def __init__ (self, gencode):
		self.gencode = gencode

	def translate_codon (self, codon, mode = 0):
		### assumes a previously cleaned up sequence (upper case, no invalid chars)
		### if (mode) returns 3 letter translation; otherwise, single letter
	
		translation = ""
		if len(codon) != 3:
			print "Error: translate_codon: not a valid codon:", codon
		else:
			first = _lookup (codon[0])
			second = _lookup (codon[1])
			third = _lookup (codon[2])
			translation = self.gencode [(first * 16) + (second * 4) + third]
			if (mode):
				translation = self.three_letter_code[translation]
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

###############################
# functions
###############################
#
# print_ruler, print_pretty_seq and find_ORFs from bioinfo.py
#
# given start = '0', prints ruler labeled 1 to 80 ...
def print_ruler (start, linelen=80):
	base = "---------+"
	ndecades = linelen / 10
	ruler = base * ndecades
	# add 3 spaces to acct for 5', 3'...
	print start+1, '\t   ' + ruler, (start + linelen)

def print_pretty_seq (seq, name, linelen=80):
	seqlen = len (seq)
	nlines = seqlen / linelen
	if ((seqlen % linelen) != 0): nlines += 1

	start = 0
	for i in range (nlines):
		c = complement (seq[start:start + linelen])
		print '\n\t5\'-' + seq [start:start + linelen]
		print_ruler (start, linelen)
		print '\t3\'-' + c
		start += linelen

	print "\nSequence:", name
	print "SeqLen:", seqlen, "bases"

def find_ORFs (seq):
	seqlen = len (seq)
	starts = []
	stop_codons = [ "TAG", "TAA", "TGA"]
	start_codon = "ATG"

	for i in range (seqlen-2):
		if seq[i:i+3] == start_codon:
			# print "Found a start codon at pos", i
			starts.append (i)

	nORFs = 0
	for start in starts:
		l = seqlen - start
		# ntriplets = l / 3
		# if (l % 3): ntriplets += 1

		i = 0
		while i < l-2:
			if seq[start+i:start+i+3] in stop_codons:
				print "\nORF found!"
				print "ORF starts at pos", start, "and ends at pos", start+i+3
				print_pretty_seq (seq[start:start+i+3], "ORF-"+str(nORFs))
				starts.remove (start)
				nORFs += 1
				break;
			i += 3

	# anything left in 'starts' is potential ORF; print to end
	for start in starts:
		print "\nPotential incomplete ORF..."
		print_pretty_seq (seq[start:], "ORF-"+str(nORFs))
		nORFs += 1

###########################
###########################

def print_sequence (seq, gencode, translate = 0, frames = [0,0,0,0,0,0], mode = 1):
### not yet printing frames 3-5 (reverse frames)
###
	labels = [ "F1", "F2", "F3", "R1", "R2", "R3" ]
	trans = [ "","","","","","" ]
	offsets = [ 0, 1, 2, 0, 1, 2 ]
#	seqlen = dna.length() ## changed input from 'dna' to 'seq'
	seqlen = len(seq)
#	seq = dna.sequence()
		
	numlines = seqlen / 80
	if (seqlen % 80): numlines += 1
	
	ruler = "    "
	for i in range (8):
		ruler += "+---------"
	
	## translate each requested reading frame
	if (translate):
		for i in range (len(frames)):
			if (frames[i]):
				trans[i] = gencode.translate_seq (seq, mode, i)
	
	## now let's print...
	for line in range (numlines):
		pos = line * 80
		print
		i = 2
		if (translate):
			while i >= 0: 
				if (frames[i]): 	
					buf = labels[i] + " "
#					buf = "   "
					if (line == 0):
						for j in range(offsets[i]):
							buf += " "						
					print buf,trans[i][pos-offsets[i] if pos-offsets[i] >=0 else pos:pos+80-offsets[i]]
				i -= 1
		print "5'-", seq[pos:80+pos]
		print ruler, pos + 80
		rev = reverse_complement (seq[pos:80+pos])
		rev = rev[::-1]
		print "3'-", rev

def reverse_complement (seq):
	seq = seq[::-1]
	rseq = ""

	for c in seq:
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
	return rseq

##########################
# main code
##########################

translate = False
cont = 1
genetic_code = GeneticCode ()
dna = DNA ("", "")

#
# True if given reading frame is on, else False
# frames 0-2 are forward, 3-5 are reverse
frames = [ True, True, True, False, False, False ] 


while (cont):
	print
	print "What would you like to do?"
	print "(1) input a new DNA sequence"
	print "(2) read DNA sequence from file"
	print "(3) find open reading frames"
	print "(4) toggle translation"
	print "(5) toggle individual reading frames"
	print "(6) supply new genetic code"
	print "(7) change single amino acid coding"
	print "(8) write sequence to file"
	print "(9) exit"

	command = raw_input ("Command: ")
	if not command.isdigit():
		print "Error: invalid command."
	else:
		command = int(command)
		if command == 1:  ## input DNA sequence
			print
			s = raw_input ("Enter a DNA sequence: ")
			dna = DNA ("", s)
			print_sequence (str(dna), genetic_code, translate, frames)
			print
			
		elif command == 2: ## read seq from file
			print
			f = raw_input ("File to open: ")
			with open (f, "r") as read_file:
				s = read_file.read ()
				print "Read: ", s
#			print "File closed." if read_file.close else "File still open!"
			dna = DNA ("", s)
			print_sequence (str(dna), genetic_code, translate, frames)
			print
				
		elif command == 3: ## find open reading frames
# 			print
# 			print "Not yet implemented."
# 			print
			seq = str(dna)
			start_codons = [ ]
			stop_codons  = [ ]
			ORFs = [ ]
			
			for i in range(len(seq)-2):
				if seq[i]+seq[i+1]+seq[i+2] == "ATG": # start codon found
					start_codons.append (i)
				elif seq[i]+seq[i+1]+seq[i+2] == "TGA":
					stop_codons.append (i)
				elif seq[i]+seq[i+1]+seq[i+2] == "TAA":
					stop_codons.append (i)
				elif seq[i]+seq[i+1]+seq[i+2] == "TAG":
					stop_codons.append (i)
				
			for i in start_codons:
				for j in stop_codons:
					if (j > i) and ((j - i) % 3 == 0):
						ORFs.append ([i,j])
						print "Found ORF starting at %d and ending at %d" % (i,j+2)

			## sort ORFs largest to smallest
			ORFs = sorted (ORFs, key = lambda orf: orf[1] - orf[0], reverse=True)
			for i in ORFs:				
				print
				print_sequence (seq[i[0]:i[1]+3], genetic_code, 1, [1,0,0,0,0,0])
			print			
			### TODO ADD LOGIC SO ORFS CONTAIN NO STOP CODONS					
					
										
 		elif command == 4: ## toggle translation
 			translate = not translate
 			print
			print_sequence (str(dna), genetic_code, translate, frames)
			print
			
		elif command == 5: ## toggle open reading frames
			labels = [ "F1", "F2", "F3", "R1", "R2", "R3" ]
			print
			print "Frames (X = on):",
			for i in range (len(frames)):
				state = "X" if frames[i] else " "
				print "%s(%d):[%s]" % (labels[i],i,state),
			print
			f = raw_input ("Select frame (0-5): ")
			if not f.isdigit() and f not in range (6):
				print "Error: invalid input."
			else:
				frame = int(f)
				frames[frame] = not frames[frame]
			print "Frames (X = on):",
			for i in range (len(frames)):
				state = "X" if frames[i] else " "
				print "%s(%d):[%s]" % (labels[i],i,state),
			print	
			print_sequence (str(dna), genetic_code, translate, frames)
			print
				
		elif command == 6: ## new genetic code
			print
			print "Your genetic code must be EXACTLY 64 characters long, one for each codon. Your"
			print "string should list each amino acid in single-letter code in the following order"
			print "order: AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA ... TTC, TTG, TTT. Use '*' "
			print "for stop codons."
			print
			print "The default is: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CYCLFLF"
			g = raw_input ("Enter a new genetic code or RETURN to leave unchanged: ")
			if g == "":
				print "Genetic code unchanged."
			elif g.isalpha() and (len(g) == 64):
				genetic_code = GeneticCode (g)
				print "Using new genetic code: ", g
			else: 
				print "Error: not a valid genetic code."
			print_sequence (str(dna), genetic_code, translate, frames)
			print
			
		elif command == 7: ## change codon
			print
			c = raw_input ("Codon to change: ")
			single_letter = raw_input ("New single-letter code: ")
			three_letter = raw_input ("New three-letter name: ")
			if ((len(c) != 3) or (len(three_letter) != 3) or (len(single_letter) != 1)):
				print "Error: invalid input."
			else:
				print
				oldcodon = genetic_code.translate_codon (c)
				string = "Confirm reassignment of " + c + " from " + oldcodon + " to " + single_letter + "? "
				confirm = raw_input (string)
				if (confirm == 'y'): 
					genetic_code.set_codon (c, single_letter, three_letter)
					print_sequence (str(dna), genetic_code, translate, frames)
				else:
					print "Codon assignment unchanged."
			print
			
		elif command == 8:
			print
			print "Not yet implemented."
			print
			
		elif command == 9:
			cont = 0
			print
			print "Good bye."
		else:
			print "Error: invalid command."
			
