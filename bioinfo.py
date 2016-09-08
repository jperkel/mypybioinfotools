s = "ACCGATCAGAACGTATGGTTATCCGTACAGACCGATGAGAACGTACGGTTATCCGTACA" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATGAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAC" \
"GACCGATCAGAACGTACGGTTATCCGTACAGACCGATCAGAACGTACGGTTATCCGTAA" \
"TGAT"

dna_bases = "ACGT"
amino_acids = "ACDEFGHIKLMNPQRSTVWY*"
stop_codons = [ "TAG", "TAA", "TGA"]
start_codon = "ATG"


def _lookup (base):
	retval = -1
	if base == 'A': retval = 0
	elif base == 'C' : retval = 1
	elif base == 'G' : retval = 2
	elif base == 'T' : retval = 3
	return retval

def filter_seq (seq, type):
	s = ""
	if type == "DNA": f = dna_bases
	else: f = amino_acids

	seq = seq.upper()
	for c in seq:
		if (c in f):
			s += c
	return s 

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

#
# given start = '0', prints ruler labeled 1 to 80 ...
def print_ruler (start, linelen=80):
	base = "---------+"
	ndecades = linelen / 10
	ruler = base * ndecades
	# add 3 spaces to acct for 5', 3'...
	print start+1, '\t   ' + ruler, (start + linelen)

def complement (seq):
	cseq = ""
	for c in seq:
		if c == 'A' : cseq += 'T'
		elif c == 'C' : cseq += 'G'
		elif c == 'G' : cseq += 'C'
		elif c == 'T' : cseq += 'A'
	return cseq

def find_ORFs (seq):
	seqlen = len (seq)
	starts = []

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

print_pretty_seq (s, "My DNA Sequence")
find_ORFs ("ATGAATGAAAAAATGAAATAAAAAAAAAAAAA")
# find_ORFs (s)

# print
# print filter_seq ("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "DNA")
# print filter_seq ("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "protein")
# print filter_seq ("abcdefghijklmnopqrstuvwxyz", "DNA")
# print filter_seq ("abcdefghijklmnopqrstuvwxyz", "protein")