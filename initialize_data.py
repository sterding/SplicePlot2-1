import pysam
import pandas
import cPickle as pickle
from lib.ReadGTF import readGTF
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
from collections import Counter
import argparse
import os


class VCFLine:
    def __init__(self,VCF_name,region):
        '''
            VCF_name is the name of the gzipped vcf file
            region is the position of the SNP, in format chr4:12345
        '''
        VCF_object = pysam.Tabixfile(VCF_name)

        VCF_header = list(VCF_object.header)
        last_header_line = VCF_header[len(VCF_header)-1]
        self.samples = last_header_line.split()[9:]

        region_list = region.split(':')

        self.contig = region_list[0]
        self.position = int(region_list[1])
        self.id = None
        self.ref = None
        self.alt = None
        self.genotype_calls = {}
        
        try:
            feature_iterator = VCF_object.fetch(region)
            for feature in feature_iterator:
                vcf_line_array = feature.strip('\n').split()
                
                contig_name = vcf_line_array[0]
                position = int(vcf_line_array[1])
                if contig_name == self.contig and position == self.position:
                    self.id = vcf_line_array[2]
                    self.ref = vcf_line_array[3]
                    self.alt = filter(lambda x: x != '.', vcf_line_array[4].split(','))
                    
                    genotype_calls_list = vcf_line_array[9:]

                    for i, indiv_id in enumerate(self.samples):
                        if '.' not in genotype_calls_list[i].split(':')[0]:
                            self.genotype_calls[indiv_id] = genotype_calls_list[i]
                    break

            if self.id == None:
                print "There is no variant at {0}".format(region)
                raise Exception

        except ValueError:
            print "{0} is not a valid SNP position for this VCF file".format(region)
            raise Exception
                

    def __getitem__(self,key):
        try:
            return self._determine_genotype_bases(self.genotype_calls[key])
        except KeyError:
            return None

    def __contains__(self,key):
        return key in self.genotype_calls

    def alleles_list(self):
        alleles = [self.ref]
        alleles.extend(self.alt)
        return alleles

    def genotypes_list(self):
        genotypes_list = []
        alleles_list = self.alleles_list()
        for i in range(len(alleles_list)):
            for j in range(i,len(alleles_list)):
                genotypes_list.append('{0}{1}'.format(alleles_list[i],alleles_list[j]))
        return genotypes_list

    def _determine_genotype_bases(self,genotype_string):
        calls = genotype_string.split(':')[0]
        calls_as_numeric = None
        if '|' in calls:
            calls_as_numeric = sorted(map(int,calls.split('|')))
        else:
            calls_as_numeric = sorted(map(int,calls.split('/')))

        genotype_string = ''.join(map(lambda x: self.alleles_list()[x],calls_as_numeric))
        return genotype_string
        
        

def average_read_depth_by_genotype(read_depth_dict,vcf_file_name,var_pos,circ_rna_data,circ_rna_junc):
    '''
        average_read_depth_by_genotype averages the ReadDepth objects in read_depth_dict based on the genotype
            in a .vcf file

        read_depth_dict is a dict of ReadDepth objects, where the keys are the individual IDs and the values are
            the corresponding ReadDepth objects

        vcf_file_name is the location of the .vcf file

        var_pos is the position of the SNP, with format chr1:12345

        return values:
            A dict containing the average read depths. The keys are the genotypes, and the values are the
                corresponding ReadDepth objects
            
            A dict which maps indiv IDs to genotype. The keys are the indiv IDs, and the values are the
                corresponding genotypes

            A list of genotypes which occur in the dataset, where homozygous reference is first

    '''

    try:
    	vcf_line = VCFLine(vcf_file_name,var_pos)
        possible_genotypes_bucket_counts = {}
        average_read_depths_dict = {}
        genotypes_in_data = set()
        genotype_by_id = {}
        for indiv_id, read_depth_object in read_depth_dict.items():
            if indiv_id in vcf_line:
                indiv_genotype = vcf_line[indiv_id] #GG GT TT
                genotype_by_id[indiv_id] = indiv_genotype
                genotypes_in_data.add(indiv_genotype)

                if indiv_genotype not in possible_genotypes_bucket_counts:
                    possible_genotypes_bucket_counts[indiv_genotype] = 1
                    average_read_depths_dict[indiv_genotype] = read_depth_object
                else:
                    possible_genotypes_bucket_counts[indiv_genotype] += 1
                    average_read_depths_dict[indiv_genotype] = average_read_depths_dict[indiv_genotype] + read_depth_object

        genotype_info = average_read_depths_dict.keys()
        # print genotype_info

        # print len(genotype_by_id), len(circ_rna_data) # 84 97
        # print circ_rna_data, circ_rna_junc
        if circ_rna_data != None and circ_rna_junc != None:
	        for idxx in range(len(circ_rna_data)):
	            circ_rna_data_ = circ_rna_data[idxx]
	            circ_rna_junc_ = circ_rna_junc[idxx]
	            # print circ_rna_junc_, circ_rna_data_
	            if circ_rna_data_ is not None and circ_rna_junc_ is not None:
	                for genotype in genotype_info:
	                    temp_count = 0
	                    for idx in genotype_by_id.keys():
	                        # print genotype_by_id[idx], genotype
	                        if genotype_by_id[idx] == genotype:
	                            if idx in circ_rna_data_.keys():
	                                # print genotype_by_id[idx], genotype, idx, circ_rna_data_[idx]
	                                temp_count += circ_rna_data_[idx]

	                    temp_junction_circ = {}
	                    temp_junction_circ[circ_rna_junc_] = temp_count
	                    # print temp_junction_circ
	                    # average_read_depths_dict[genotype].circRNA[idxx] = temp_junction_circ
	                    average_read_depths_dict[genotype].circRNA.update(temp_junction_circ)
	                    # print average_read_depths_dict[genotype].circRNA
	        print average_read_depths_dict

        for genotype, counts in possible_genotypes_bucket_counts.items():
            # print genotype, counts #GG 76 TT 1 GT 7
            average_read_depths_dict[genotype].divide_by_constant(counts)
            # print average_read_depths_dict[genotype]

        # normalize read depths
        # print 'Normalizing average read depths...'
        # for genotype in genotype_info:
        #     dist = average_read_depths_dict[genotype].wiggle
        #     dist_norm = [float(i)*1.5/max(dist) for i in dist]
        #     average_read_depths_dict[genotype].wiggle = dist_norm

        filtered_genotypes_list = filter(lambda x: x in genotypes_in_data, vcf_line.genotypes_list())
        return average_read_depths_dict, genotype_by_id, filtered_genotypes_list, possible_genotypes_bucket_counts

    except IOError:
        print 'There is no .vcf file at {0}'.format(vcf_file_name)
        raise Exception


def map_indiv_id_to_bam_name(id_map_file):

    '''
        map_indiv_id_to_bam_name maps the individual ids in the vcf file to the corresponding .bam files

        id_map_file is the name of the file containing the id to bam file correspondences

        return value:
            a dict that maps individual id to .bam file location

            a list of all the .bam file locations to process
    '''

    bam_to_id = {}
    bam_list = []

    try:
        lines = open(id_map_file,'r').readlines()
        for line in lines:
            line = line.strip('\n').split()

            file_path = line[1]
            if os.path.exists(file_path):
                bam_to_id[line[1]] = line[0]
                bam_list.append(line[1])
            else:
                print 'Skipping {0}. Invalid file path'.format(file_path)

        return bam_to_id, bam_list

    except IOError:
        print 'There is no mapping file at {0}'.format(id_map_file)
        raise Exception


class Exon:
    def __init__(self,chrm,low,high,strand):
        self.chrm = chrm
        self.low = low
        self.high = high
        self.strand = strand

    @classmethod
    def create_from_gtf(cls,line):
        '''
            create_from_gtf creates an Exon object from a single line of a gtf file

            line is a single line of text (a string) from a gtf file
        '''
        info = line.strip('\n').split()
        return cls(info[0],int(info[3]),int(info[4]),info[6])

    def determine_proportion_covered(self,read_depth):
        if read_depth.chrm != self.chrm or self.low < read_depth.low or self.high > read_depth.high:
            return 0
        
        # determine the top and bottom indices
        bottom_index = self.low - read_depth.low
        top_index = bottom_index + (self.high - self.low)

        bases_covered = 0
        for item in range(bottom_index, top_index + 1):
            if read_depth.wiggle[item] > 0:
                bases_covered += 1.0

        return bases_covered / (top_index + 1 - bottom_index)

    def determine_average_coverage(self,read_depth):
        return self.determine_proportion_covered(read_depth) / (self.high - self.low + 1.0)

    def __str__(self):
        return '{0}:{1}-{2}{3}'.format(self.chrm,self.low,self.high,self.strand)


class EvaluatedExon:
    def __init__(self,exon,value):
        self.exon = exon
        self.value = value
        self.length = exon.high - exon.low + 1

    def __cmp__(self,other):
        if self.value - other.value != 0:
            return self.value - other.value

        return self.length - other.length

    def __str__(self):
        return '{0},{1}'.format(self.exon.__str__(), self.value)


def get_splice_range_coordinates(range_name):
    try:
        chrom = range_name.split(':')[0]
        start_idx, end_idx = map(int, range_name.split(':')[1].split('-'))

        return chrom, start_idx, end_idx
    except:
        print '{0} is not a valid range name'.format(range_name)
        raise Exception

                
def determine_exons_and_coordinates(gtf_file, chrom, start_idx, end_idx):

    try:
        gtf_list = readGTF(gtf_file)
        relevant_exons_iterator = []

        for gtf_line in gtf_list:
            # print gtf_line.start, gtf_line.end
            if gtf_line.start >= start_idx and gtf_line.end <= end_idx:
                relevant_exons_iterator.append(gtf_line)
            elif gtf_line.start >= start_idx and gtf_line.start <= end_idx:
                gtf_line.end = end_idx
                relevant_exons_iterator.append(gtf_line)
            elif gtf_line.end >= start_idx and gtf_line.end <= end_idx:
                gtf_line.start = start_idx
                relevant_exons_iterator.append(gtf_line)

        filtered_exons_list = []
        min_coordinate = float('inf')
        max_coordinate = float('-inf')
        all_exons_idx = []
        splice_junc_coordinate_list = []

        for exon in relevant_exons_iterator:
            print "exon:{}-{}".format(exon.start, exon.end)
            filtered_exons_list.append(exon)
            all_exons_idx.append(exon.start)
            all_exons_idx.append(exon.end)

        min_coordinate = min(all_exons_idx)
        max_coordinate = max(all_exons_idx)
        shared_site = all_exons_idx[1]
        other_sites = all_exons_idx[2::2]

        if min_coordinate == float('inf') or max_coordinate == float('-inf'):
            print 'The given range coordinates do not correspond to exons in the annotation'
            raise Exception

        # print min_coordinate, max_coordinate, shared_site, other_sites
        return min_coordinate, max_coordinate, shared_site, other_sites, filtered_exons_list
    except IOError:
        print 'There is no gtf file at {0}'.format(gtf_file)
        raise Exception 


def initialize_read_depths_and_determine_exons(junction_name,gtf_file_name,bam_list,bam_to_id_dict):

    '''
        initialize_read_depths_and_determine_exons creates a dictionary of read depths and assembles a possible set of exons

        junction_name is the name of the junction being examined. It is a string with the 
            format chr1:17055-17915,chr1:17055-17606,chr1:17055-17233,
            where the numbers represent the genomic coordinates of the splice sites

        gtf_file_name is a string representing the path to the gtf file containing known exons

        bam_list is a list of strings containing the file paths to the bam files

        return values:
            A dictionary of ReadDepth objects, where the key is the file path (as a string) to the bam file,
                and where the value is a ReadDepth object corresponding to the bam file

            A mRNAsObject, which represents the set of possible mRNA segments determined from the junctions
    '''

    chrom, start_idx, end_idx = get_splice_range_coordinates(junction_name)

    minimum_coordinate,maximum_coordinate,shared_site,other_sites,filtered_exons_list = determine_exons_and_coordinates(gtf_file_name,chrom,start_idx,end_idx)

    # print minimum_coordinate, maximum_coordinate
    # print start_idx, end_idx
    minimum_coordinate = start_idx
    maximum_coordinate = end_idx

    read_depth_dict = {}
    for bam_file in bam_list:
        current_read_depth = ReadDepth.determine_depth(bam_file,chrom,minimum_coordinate,maximum_coordinate)
        print bam_file
        read_depth_dict[bam_to_id_dict[bam_file]] = current_read_depth
    
    mRNAs = []
    resize_lower_bound = minimum_coordinate
    resize_upper_bound = maximum_coordinate

    list_id = []
    for i in range(len(filtered_exons_list)):
        # print filtered_exons_list[i].ID
        list_id.append(filtered_exons_list[i].ID)
    id_count = Counter(list_id)

    exon = filtered_exons_list[0]

    for key in id_count:
        exon_splice = []
        for idx in range(len(filtered_exons_list)):
            temp_exon = [filtered_exons_list[idx].start, filtered_exons_list[idx].end]
            if filtered_exons_list[idx].ID == key:
                exon_splice.append(temp_exon)

        mRNAs.append(exon_splice)
   
    # resize each of the read depth objects so that their lengths correspond to the possible mRNAs
    for key in read_depth_dict:
        read_depth_dict[key].shrink(resize_lower_bound,resize_upper_bound)

    mRNAs_info = mRNAsObject(exon.seqname,exon.strand,resize_lower_bound,resize_upper_bound,mRNAs)
    return read_depth_dict, mRNAs_info


def create_data_frame(read_depth_dict,junction_name,var_pos,genotype_lookup_dict,filtered_genotypes_list):

    '''
        create_data_frame creates a pandas.DataFrame object in the format required by the hive
            and structure plotting functions

        read_depth_dict is a dict, where the keys are the bam file paths and the values are ReadDepth objects

        junction_name is the name of the junction, in the usual format

        genotype_lookup_dict is a dict, where the key is the indiv ID (as specified in .vcf file)
            and where the value is the individual's genotype

        var_pos is the position of the SNP, in the format chr1:12345

        return values:
            A pandas.DataFrame. Each row represents an individual, and the row name is the indiv id.
                The first column contains the genotypes of each individual. All remaining columns contain
                the splicing ratios corresponding to a splice junction

    '''

    precursor_dict = {}
    data_frame_index = []
    junctions_list = junction_name.split(',')
    for key, value in read_depth_dict.items():

        indiv_genotype = genotype_lookup_dict[key]

        data_frame_index.append(key)

        # add the genotype to precursor_dict
        if var_pos not in precursor_dict:
            precursor_dict[var_pos] = []
        precursor_dict[var_pos].append(indiv_genotype)

        # calculate the splicing ratios
        total_junc_read_count = 0
        for junction in junctions_list:
            if junction in value.junctions_dict:
                total_junc_read_count += value.junctions_dict[junction]

        # add splicing ratios to precursor_dict
        for junction in junctions_list:
            write_value = 0.0

            if junction in value.junctions_dict:
                write_value = value.junctions_dict[junction] * 1.0 / total_junc_read_count

            if junction not in precursor_dict:
                precursor_dict[junction] = []

            precursor_dict[junction].append(write_value)

    # rearrange the data frame so that genotype is the first column
    df = pandas.DataFrame(precursor_dict,index=data_frame_index)
    new_col_order = [var_pos]
    new_col_order.extend(junctions_list)
    return df.reindex(columns=new_col_order)


def read_circ_rna(data_file, junction_name):
    chrom, start_idx, end_idx = get_splice_range_coordinates(junction_name)
    # print chrom, start_idx, end_idx

    circRNA_all = {}
    junctions_all = {}
    try:
        lines = open(data_file,'r').readlines()
        genotype = lines[0].strip('\n').split()[1:]
        # print genotype
        index = 0
        for idx_ in range(1, len(lines)):
            circRNA = {}
            temp_junction = None
            read_count = lines[idx_].strip('\n').split()[1:]
            junction = lines[idx_].strip('\n').split()[0]
            for idx in range(len(genotype)):
                temp_genotype = genotype[idx].split('_')[1]
                circRNA[temp_genotype] = int(read_count[idx])

            circ_start_idx = junction.split('_')[1]
            circ_end_idx = junction.split('_')[2]
            # print circ_start_idx, circ_end_idx
            if int(circ_start_idx) >= int(start_idx) and int(circ_end_idx) <= int(end_idx):
                temp_junction = junction.split('_')[0]+':'+junction.split('_')[1]+'-'+junction.split('_')[2]
                circRNA_all[index] = circRNA
                junctions_all[index] = temp_junction
                index = index + 1
        # print circRNA_all,junctions_all
        return circRNA_all, junctions_all

    except IOError:
        print 'There is no circle RNA file at {0}'.format(data_file)
        raise Exception


def calculate_average_expression_and_data_frame(var_pos,junction_name,vcf,annotation,map_file,circ_rna_file):
    #add read circle RNA data
    if circ_rna_file is not None:
        circ_rna_data, circ_rna_junc = read_circ_rna(circ_rna_file, junction_name)
    else:
        circ_rna_data = None
        circ_rna_junc = None

    # print circ_rna_data, circ_rna_junc

    bam_to_id_dict, bam_list = map_indiv_id_to_bam_name(map_file)
    
    read_depths_dict, mRNA_info_object = initialize_read_depths_and_determine_exons(junction_name,
        annotation,bam_list,bam_to_id_dict)

    new_read_depths_dict = read_depths_dict
    # print new_read_depths_dict
    genotype_averages_dict, genotype_by_id, filtered_genotypes_list, genotype_dict_counts = average_read_depth_by_genotype(new_read_depths_dict,vcf,var_pos,circ_rna_data,circ_rna_junc)
    # print genotype_averages_dict
    data_frame = create_data_frame(new_read_depths_dict,junction_name,var_pos,genotype_by_id,filtered_genotypes_list)
    # print data_frame

    return genotype_averages_dict, data_frame, mRNA_info_object, filtered_genotypes_list, genotype_dict_counts


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Initialize and pickle alternative splice junction data')

    parser.add_argument('varpos',type=str,help='string describing position of SNP. should have format chr_name:base_number')
    parser.add_argument('junc',type=str,help="string representing junction. should have format chr_name:lower_base-upper_base,chr_name:lower_base-upper_base, where lower_base and upper_base represent possible intronic regions")
    parser.add_argument('--vcf',type=str,required=True,help='location of the vcf file')
    parser.add_argument('--gtf',type=str,required=True,help='location of the gtf file')
    parser.add_argument('--mf',type=str,required=True,help='location of the map file')
    parser.add_argument('--circ',type=str,required=False,help='location of the circle-RNA file')
    parser.add_argument('--output',type=str,required=False,default=None,help='location of output pickle file. Optional parameter')

    args = parser.parse_args()
    try:
        if args.circ is not None:
            genotype_averages_dict, data_frame, mRNA_info_object, filtered_genotypes_list, genotype_dict_counts = calculate_average_expression_and_data_frame(args.varpos,args.junc,args.vcf,args.gtf,args.mf, args.circ)
        else:
            genotype_averages_dict, data_frame, mRNA_info_object, filtered_genotypes_list, genotype_dict_counts = calculate_average_expression_and_data_frame(args.varpos,args.junc,args.vcf,args.gtf,args.mf, None)

        output_file_path = '{0}/pickle_files/'.format(os.path.dirname(os.path.abspath(__file__)))
        if args.output is not None:
            output_file_path = args.output

        stem = output_file_path
        tail = '{0}@{1}.p'.format(args.varpos,args.junc)
        
        if output_file_path[len(output_file_path)-2:] == '.p':
            stem, tail = os.path.split(output_file_path)

        try:
            os.makedirs(stem)
        except OSError:
            if os.path.isdir(stem):
                pass
            else:
                print 'Cannot create directory {0}'.format(stem)
                raise Exception

        if stem != '' and stem[len(stem)-1] != '/':
            stem = stem + '/'

        pickle_file = open('{0}{1}'.format(stem,tail),'wb')
        pickle.dump(args.varpos,pickle_file)
        pickle.dump(args.junc,pickle_file)
        pickle.dump(genotype_averages_dict,pickle_file)
        pickle.dump(mRNA_info_object,pickle_file)
        pickle.dump(data_frame,pickle_file)
        pickle.dump(filtered_genotypes_list,pickle_file)
        pickle.dump(genotype_dict_counts, pickle_file)

        pickle_file.close()
        print 'Done!'
    except:
        print 'Failed'
