# bioinformatics tools

debug = 0


##########################
# functions
##########################
def reverse_complement (seq, type='dna'):
	seq = seq[::-1]
	rseq = ""
	
	for c in seq:
		if c == 'A':
			if (type == 'dna'):
				rseq += 'T'
			elif (type == 'rna'):
				rseq += 'U'
		elif c == 'C':
			rseq += 'G'
		elif c == 'G':
			rseq += 'C'
		elif c == 'T' or c == 'U':
			rseq += 'A'
		else:
			print "Error: reverse_complement: invalid character:", c
	return rseq

def filter_sequence (seq, filter):
	s = seq.upper ()
	newseq = ""
	for c in s:
		if c in filter:
			newseq += c
	return newseq

########################
# helper function for translate_codon
########################
def _lookup (nt):
	if nt == 'A':
		return 0
	elif nt == 'C':
		return 1
	elif nt == 'G':
		return 2
	return 3	
	
def translate_codon (codon, gencode, mode = 0):
	### assumes a previously cleaned up sequence (upper case, no invalid chars)
	### if (mode) returns 3 letter translation; otherwise, single letter
	
	translation = ""
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
	
def translate_seq (seq, gencode, type = "dna", mode = 0, frame = 0):
	### assumes a previously cleaned up sequence
	### if (mode) returns 3 letter translation; otherwise, single letter
	### frame: 0-2 forward dir, 3-5 reverse dir
		
	translation = [ ]
	codons = [ ]
	offsets = [ 0,1,2,0,1,2 ]
	
	if frame not in range (6):
		print "Error: translate_seq: invalid frame", frame
		return translation
	else:
		if frame >= 3: 
			seq = reverse_complement (seq,type)

		offset = offsets[frame]
		ncodons = (len(seq)-offset) / 3		
		
		seq = seq[offset:]
		if (debug): print "translate_seq: working seq = ", seq
		for i in range(ncodons):			
			codon = seq[(i*3)]+seq[(i*3)+1]+seq[(i*3)+2]
			if (debug) : print "translate_seq: found codon",codon
			codons.append (codon) 
		for i in codons:
			codon = translate_codon (i, gencode, mode)
			translation.append (codon)
		translation = ''.join(translation)
		if (debug): print "translate_seq: translation =", translation
		return translation
	
def print_sequence (seq, gencode, translate = 0, mode = 1):
	trans = [ "","","","","","" ]
	offsets = [ 0, 1, 2, 0, 1, 2 ]
	seqlen = len (seq)

	numlines = seqlen / 80
	if (seqlen % 80): numlines += 1
	
	ruler = "    "
	for i in range (8):
		ruler += "+---------"
	
	## translate each requested reading frame
	if (translate):
		for i in range (len(frames)):
			if (frames[i]):
				trans[i] = translate_seq (seq, gencode, "dna", mode, i)
	
	## now let's print...
	for line in range (numlines):
		pos = line * 80
		print
		if (translate):
			for i in range (len(frames)): 
				if (frames[i]): 	
					buf = "   "
					if (line == 0):
						for j in range(offsets[i]):
							buf += " "						
					print buf,trans[i][pos-offsets[i] if pos-offsets[i] >=0 else pos:pos+80-offsets[i]]
		print "5'-", seq[pos:80+pos]
		print ruler, pos + 80
		rev = reverse_complement (seq[pos:80+pos])
		rev = rev[::-1]
		print "3'-", rev


##########################
# universal variables
##########################
valid_nucleotides = "ACGTU"
valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

##########################
# AAA, AAC, AAG, AAU, ACA, ACC, ACG, ACU, AGA ...
# A = 0, C = 1, G = 2, T/U = 3
# Index = (16 * first) + (4 * second) + third 
##########################
std_genetic_code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CYCLFLF"

three_letter_code = { 'A' : 'Ala', 'C' : 'Cys', 'D' : 'Asp', 'E' : 'Glu', 'F' : 'Phe', \
'G' : 'Gly', 'H' : 'His', 'I' : 'Ile', 'K' : 'Lys', 'L' : 'Leu', 'M' : 'Met', 'N' : 'Asn', \
'P' : 'Pro', 'Q' : 'Gln', 'R' : 'Arg', 'S' : 'Ser', 'T' : 'Thr', 'V' : 'Val', 'W' : 'Trp', \
'Y' : 'Tyr', '*' : '***' }

##########################
# code starts
##########################

translate = False
cont = 1
genetic_code = std_genetic_code

#
# True if given reading frame is on, else False
# frames 0-2 are forward, 3-5 are reverse
frames = [ True, False, False, False, False, False ] 


while (cont):
	print
	print "What would you like to do?"
	print "(1) input a new DNA sequence"
	print "(2) find open reading frames"
	print "(3) toggle translation"
	print "(4) toggle reading frame"
	print "(5) supply new genetic code"
	print "(9) exit"

	command = raw_input ("Command: ")
	if not command.isdigit():
		print "Error: invalid command."
	else:
		command = int(command)
		if command == 1:
			print
			seq = raw_input ("Enter a DNA sequence: ")
			if (debug): print "Filtering sequence..."
			seq = filter_sequence (seq, valid_nucleotides)
			print_sequence (seq, genetic_code, translate)
			
		elif command == 2:
			print "Not yet implemented."

 		elif command == 3:
 			translate = not translate
			print_sequence (seq, genetic_code, translate)

		elif command == 4:
			f = raw_input ("Select frame (0-5): ")
			if not f.isdigit() and f not in range (6):
				print "Error: invalid input."
			else:
				frame = int(f)
				frames[frame] = not frames[frame]
				print "Reading frame %d changed." % frame
				print_sequence (seq, genetic_code, translate)
				
		elif command == 5:
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
				genetic_code = g
				print "Using new genetic code: ", genetic_code
			else: 
				print "Error: not a valid genetic code."
			
		elif command == 9:
			cont = 0
			print
			print "Good bye."
		else:
			print "Error: invalid command."
			
			
# seq1 = "acgtacgtacgtacgt"
# seq2 = "abcdefghijklmnopqrstuvwxyz"
# 
# print "Filtering 'acgtacgtacgtacgt..."
# print "for DNA bases:", filter_sequence (seq1, valid_dna_bases)
# print "for amino acids:", filter_sequence (seq1, valid_amino_acids)
# seq1 = filter_sequence (seq1, valid_dna_bases)
# print
# 
# print "Filtering abcdefghijklmnopqrstuvwxyz..."
# print "for DNA bases:", filter_sequence (seq2, valid_dna_bases)
# print "for RNA bases:", filter_sequence (seq2, valid_rna_bases)
# print "for amino acids:", filter_sequence (seq2, valid_amino_acids)
# seq2 = filter_sequence (seq2, valid_amino_acids)
# print
# 
# print "The reverse complement of AACCTTGG is", reverse_complement ("AACCTTGG")
# print
# 
# print "AUG =", translate_codon ("AUG", genetic_code, 1)
# print "UUU =", translate_codon ("UUU", genetic_code, 1)
# print "CUG =", translate_codon ("CUG", genetic_code, 1)
# print "GGG =", translate_codon ("GGG", genetic_code, 1)
# print "GGGG =", translate_codon ("GGGG", genetic_code, 1)
# print
# 
# print "AUGUUUCUGGGGGGGTGACC =", translate_seq ("AUGUUUCUGGGGGGGTGACC", genetic_code)
# print "AUGUUUCUGGGGGGGTGACC =", translate_seq ("AUGUUUCUGGGGGGGTGACC", genetic_code, 1)
# print
# 
# print "Translationg ATGTTTCTGGGGGGGTGACC in all 6 reading frames..."
# print "F1:",translate_seq ("ATGTTTCTGGGGGGGTGACC", genetic_code, 1, 0)
# print "F2:",translate_seq ("ATGTTTCTGGGGGGGTGACC", genetic_code, 1, 1)
# print "F3:",translate_seq ("ATGTTTCTGGGGGGGTGACC", genetic_code, 1, 2)
# print "R1:",translate_seq ("ATGTTTCTGGGGGGGTGACC", genetic_code, 1, 3)
# print "R2:",translate_seq ("ATGTTTCTGGGGGGGTGACC", genetic_code, 1, 4)
# print "R3:",translate_seq ("ATGTTTCTGGGGGGGTGACC", genetic_code, 1, 5)
