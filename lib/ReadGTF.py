import pysam
import argparse

class GTFLine:
    def __init__(self,text):
        text_list = text.strip('\n').split('\t')
        self.seqname = text_list[0]
        self.source = text_list[1]
        self.feature = text_list[2]
        self.start = int(text_list[3])
        self.end = int(text_list[4])
        self.score = text_list[5]
        self.strand = text_list[6]
        self.frame = text_list[7]
        self.ID = text_list[8].split(';')[1]
        self.attribute = None

        if len(text_list) > 8:
            self.attribute = text_list[8]

    def __eq__(self,other):
        return self.seqname == other.seqname and self.start == other.start and \
            self.end == other.end

    def __cmp__(self,other):
        if self.seqname < other.seqname:
            return -1
        elif self.seqname > other.seqname:
            return 1
        else:
            if self.start - other.start != 0:
                return self.start - other.start
            else:
                return self.end - other.end

    def __str__(self):
        if self.attribute is not None:
            return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(self.seqname,
                self.source,self.feature,self.start,self.end,self.score,self.strand,
                self.frame,self.attribute)
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(self.seqname,
                self.source,self.feature,self.start,self.end,self.score,self.strand,
                self.frame)

    def __hash__(self):
        return hash((self.seqname,self.start,self.end))


def readGTF(gtf_file):
	f1 = open(gtf_file,'r').readlines()
	feature_list = []

	for line in f1:
		l2 = line.strip('\n')
		try:
			gtf_line = GTFLine(l2)
			if gtf_line.feature == 'exon':
				# print gtf_line.start, gtf_line.end, gtf_line.ID
				feature_list.append(gtf_line)
		except Exception as e:
			pass
	return feature_list
