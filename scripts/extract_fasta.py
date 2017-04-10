# extract_BLAST_alignment.py
# TO RUN EXECUTE
# python extract_BLAST_alignment.py [input_fasta] [DB]

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import sys

query_fasta = sys.argv[1]
QUERY_GENE = query_fasta.split('/')[-1].split('.')[0]
DB = sys.argv[2]
blastx_cline = NcbiblastxCommandline(cmd='blastn',query=query_fasta, db=DB,
                                     qcov_hsp_perc=60,max_target_seqs=1000,outfmt=5,out=query_fasta+'blast_out.xml')
blastx_cline()
blast_records = NCBIXML.read(open(query_fasta+'blast_out.xml'))
output = open(QUERY_GENE+'_blast_hits.fasta', 'w')
q=SeqIO.read(query_fasta, 'fasta')
SeqIO.write(q,output, 'fasta')

for alignment in blast_records.alignments:
	for hsp in alignment.hsps:
		if (len(q.seq)==hsp.query_end) and (hsp.query_start==1):
			seq_str = hsp.sbjct
		elif (len(q.seq)==hsp.query_end) and (hsp.query_start>1):
			seq_str = "-"*len(str(q.seq[:hsp.query_start-1]))+str(hsp.sbjct)
		elif len(q.seq) > hsp.query_end:
			seq_str = "-"*len(str(q.seq[:hsp.query_start-1]))+str(hsp.sbjct)+"-"*len(q.seq[hsp.query_end:])
		record = SeqRecord(Seq(seq_str),id=QUERY_GENE, name=alignment.hit_def, description=alignment.hit_def)
		SeqIO.write(record,output, 'fasta')
		print(alignment.title)
		print("gaps ",hsp.gaps)
		print(hsp.query_start,hsp.query_end)
		print("---------------------")



output.close()
