## bioinformatics tools v3
from GeneticCode import *
from Ribosome import *
from BioPolymer import *
from ORFFinder import *
# import pdb

###############################
# functions
###############################

# given start = '0', prints ruler labeled 1 to 80 ...
def print_ruler (start, linelen=80):
	# pdb.set_trace()
	base = "---------+"
	ndecades = linelen / 10
	ruler = base * ndecades
	# add 3 spaces to acct for 5', 3'...
	print start+1, '\t   ' + ruler, (start + linelen)


def find_ORFs (seq, gencode):
	seqlen = len (seq)
	starts = []

	# find start and stop codons for this genetic code...
	start_codons = [ ]
	stop_codons = [ ]
	for i,c in enumerate (str(gencode)):
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
				# print_dna_sequence (seq[start:start+i+3], "ORF-"+str(nORFs), gencode)
				starts.remove (start)
				break;
			i += 3

	for ORF in ORFs:
		print "\nORF starts at pos", ORF[0], "and ends at pos", ORF[1]	
	return ORFs

	# anything left in 'starts' is potential ORF; print to end
	# for start in starts:
	# 	print "\nPotential incomplete ORF..."
	# 	print_dna_sequence (seq[start:], "ORF-"+str(nORFs), gencode)
	# 	nORFs += 1


def print_dna_sequence (dna, name, gencode, mode = 1, linelen=80):
### not yet printing reverse frames
###
	labels = [ "F1", "F2", "F3", "R1", "R2", "R3" ]
	trans = [ "","","","","","" ]
	offsets = [ 0, 1, 2, 0, 1, 2 ]

	seqlen = dna.length()
	seq = str(dna)
	ribosome = Ribosome(gencode)

	nlines = seqlen / linelen
	if ((seqlen % linelen) != 0): nlines += 1

	## translate each requested reading frame
	if (translate):
		for i in range (len(frames)):
			if (frames[i]):
				trans[i] = ribosome.translate_seq (seq[offsets[i]:], mode)
	
	start = 0
	print "\nSequence:", name
	for line in range (nlines):
		pos = line * 80
		print
		i = 2
		if (translate):
			while i >= 0: 
				if (frames[i]): 	
					buf = "\t" + labels[i] 
#					buf = "   "
					if (line == 0):
						for j in range(offsets[i]):
							buf += " "						
					print buf,trans[i][pos-offsets[i] if pos-offsets[i] >=0 else pos:pos+80-offsets[i]]
				i -= 1

		c = complement (seq[start:start + linelen])
		print '\t5\'-' + seq [start:start + linelen]
		print_ruler (start, linelen)
		print '\t3\'-' + c
		start += linelen

	print "SeqLen:", seqlen, "bases"

def print_protein_sequence (protein, name, mode = 1, linelen=80):

	seqlen = protein.length()
	seq = str(protein)

	nlines = seqlen / linelen
	if ((seqlen % linelen) != 0): nlines += 1

	start = 0
	print "\nSequence:", name
	for line in range (nlines):
		pos = line * 80
		print
		print '\tN- ' + seq [start:start + linelen]
		print_ruler (start, linelen)
		# print '\t3\'-' + c
		start += linelen

#
# complement() returns to complement of seq in 3'-5' direction
def complement (seq):
	cseq = ""
	for c in seq:
		if c == 'A' : cseq += 'T'
		elif c == 'C' : cseq += 'G'
		elif c == 'G' : cseq += 'C'
		elif c == 'T' : cseq += 'A'
	return cseq


##########################
# main code
##########################

translate = False
cont = 1
genetic_code = GeneticCode ()
# genetic_code.code()

dna = DNA ("")
name = ""

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
	print "(8) write sequence to file (* not implemented *)"
	print "(9) transcribe"
	print "(10) reverse complement"
	print "(11) exit"

	command = raw_input ("Command: ")
	if not command.isdigit():
		print "Error: invalid command."
	else:
		command = int(command)
		if command == 1:  ## input DNA sequence
			print
			name = raw_input ("Sequence name: ")
			s = raw_input ("Enter a DNA sequence: ")
			dna = DNA (s, name)
			print_dna_sequence (dna, name, genetic_code)
			print "MW:", dna.molwt()
			print
			
		elif command == 2: ## read seq from plain text file
			print
			f = raw_input ("File to open: ")
			with open (f, "r") as read_file:
				s = read_file.read ()
				print "Read: ", s
#			print "File closed." if read_file.close else "File still open!"
			dna = DNA (s)
			print_dna_sequence (dna, name, genetic_code)
			print
				
		elif command == 3: ## find open reading frames
			ORFs = [ ]
			finder = ORFFinder (genetic_code)
			ribosome = Ribosome (genetic_code)

			ORFs = finder.find(dna)
			nORFs = len (ORFs)
			if nORFs:
				print
				print "Found %d ORFs\n" % (nORFs)

				i = 0
				for ORF in ORFs:
					print "[%d]" % i, 
					print "length:", ((ORF[1] - ORF[0]) / 3) - 1, "amino acids"
					# print ribosome.translate_seq (str(dna)[ORF[0]:ORF[1]], 0)
					i += 1

				o = raw_input ("Select an ORF: ")
				ORF = ORFs[int(o)]
				p = Protein (ribosome.translate_seq (str(dna)[ORF[0]:ORF[1]], 0), 
					name + ".orf-" + o)
				print_protein_sequence (p, p.seqname())
				print 
				print "Length: ", p.length(), "amino acids"
				print "Mol Wt: ", p.molwt (), "Da"
			else:
				print "Found zero ORFs\n"
										
 		elif command == 4: ## toggle translation
 			translate = not translate
 			print
			print_dna_sequence (dna, name, genetic_code)
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
			print_dna_sequence (dna, name, genetic_code)
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
			print_dna_sequence (dna, name, genetic_code)
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
					print_dna_sequence (dna, name, genetic_code)
				else:
					print "Codon assignment unchanged."
			print
			
		elif command == 8:
			print
			print "Not yet implemented."
			print
			
		elif command == 9:
			rna = dna.transcribe()
			print str(rna)
			print

		elif command == 10:
			dna = dna.reverse_complement()
			print_dna_sequence (dna, dna.seqname(), genetic_code)
			print

		elif command == 11:
			cont = 0
			print
			print "Good bye."

		else:
			print "Error: invalid command."
			
