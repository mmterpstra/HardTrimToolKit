#!/usr/bin/env python3
import argparse
import pysam
import sys
import gzip
import pprint
import time
from datetime import timedelta
import logging
from hardtrimtoolkit import version
version = version.__version__
def main():
    hardtrim_adapters_cli()
def hardtrim_adapters_cli():
    start_time = time.time()
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='# %(asctime)s ##INFO## %(message)s')
    logging.info("Started "+sys.argv[0]+" version "+version+".")
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description=' Takes a readname sorted bam file and uses the adapter trimming tag to change the soft trims into hard trims. Stores original optionally in the TR and TQ tags.',
        epilog='Example useage python3 marktrimming.py --input input.bam --output output.bam --fastq qname_R1.fastq.gz --fastq qname_R2.fastq.gz')
    parser.add_argument('--output',help='Output bam. Should support "-" as /dev/stdout', required=True)
    parser.add_argument('--input', help='Input bam. Should support "-" as /dev/stdin', required=True )
    parser.add_argument('--ubam-out',help='Write uncompressed output bam. Helpful for speeding up on stream processing',  type=bool, required=False, default=False )
    parser.add_argument('--ubam-in',help='Read uncompressed input bam.', type=bool, required=False, default=False )
    parser.add_argument('--trimming-tag',help='Tag guide in trimmings. Default is XT for gatk compatibility', type=str, required=False, default="XT")
    parser.add_argument('--original-alignment-tag',help='Tag to move original alignment to. Default is OA for gatk compatibility', type=str, required=False, default="OA")
    parser.add_argument('--version', action='version', version='%(prog)s '+version)


    args = parser.parse_args()
    #main subroutine
    count,trimmed_count_first_of_pair,trimmed_count_second_of_pair = hardtrim_adapters(args)
    
    logging.info("Ended script.")
    logging.info(" Run time " + str(timedelta(seconds=time.time()- start_time)) + ".")
    logging.info(" Records processsed "+str(count) + ".")
    logging.info(" Trimmed first of pair "+ str(trimmed_count_first_of_pair) + ".")
    logging.info(" Trimmed second of pair "+ str(trimmed_count_second_of_pair) + ".")

def hardtrim_adapters(args):
    print("To be added")
    #from https://samtools.github.io/hts-specs/SAMtags.pdf part 1.4 (on 20250220): 
    #OA:Z:(RNAME ,POS ,strand ,CIGAR,MAPQ,NM ;)+ The original alignment information of the record
    #prior to realignment or unalignment by a subsequent tool. Each original alignment entry contains
    #the following six field values from the original record, generally in their textual SAM representations,
    #separated by commas (‘,’) and terminated by a semicolon (‘;’): RNAME, which must be explicit
    #(unlike RNEXT, ‘=’ may not be used here); 1-based POS; ‘+’ or ‘-’, indicating forward/reverse strand
    #respectively (as per bit 0x10 of FLAG); CIGAR; MAPQ; NM tag value, which may be omitted (though
    #the preceding comma must be retained).
    #In the presence of an existing OA tag, a subsequent tool may append another original alignment entry
    #after the semicolon, adding to—rather than replacing—the existing OA information.
    #The OA field is designed to provide record-level information that can be useful for understanding the
    #provenance of the information in a record. It is not designed to provide a complete history of the
    #template alignment information. In particular, realignments resulting in the the removal of Secondary
    #or Supplementary records will cause the loss of all tags associated with those records, and may also
    #leave the SA tag in an invalid state.
    #OC:Z:cigar Original CIGAR, usually before realignment. Deprecated in favour of the more general OA.
    #OP:i:pos Original 1-based POS, usually before realignment. Deprecated in favour of the more general OA.
    #OQ:Z:qualities Original base quality, usually before recalibration. Same encoding as QUAL.

    #psuedocode:
    #read rn sorted bam for same name reads then until other name spotted and save that to buffer.
        #store og data like mate cigar, cigar and so on
        #hardtrim read,qual,cigar based on trimming tag(XT) and update alignment position if needed and remove trailing indels on cigar.


    return(0,0,0)

def revert_hardtrim_adapters_cli():
    start_time = time.time()
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='# %(asctime)s ##INFO## %(message)s')
    logging.info("Started "+sys.argv[0]+" version "+version+".")
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description=' Takes a readname sorted bam file and uses the adapter trimming tag to change the soft trims into hard trims. Stores original optionally in the TR and TQ tags.',
        epilog='Example useage python3 marktrimming.py --input input.bam --output output.bam --fastq qname_R1.fastq.gz --fastq qname_R2.fastq.gz')
    parser.add_argument('--output',help='Output bam. Should support "-" as /dev/stdout', required=True)
    parser.add_argument('--input', help='Input bam. Should support "-" as /dev/stdin', required=True )
    parser.add_argument('--ubam-out',help='Write uncompressed output bam. Helpful for speeding up on stream processing',  type=bool, required=False, default=False )
    parser.add_argument('--ubam-in',help='Read uncompressed input bam.', type=bool, required=False, default=False )
    parser.add_argument('--trimming-tag',help='Tag to add in trimmings. Default is XT for gatk compatibility', type=str, required=False, default="XT")
    parser.add_argument('--original-alignment-tag',help='Tag to move original alignment from. Default is OA for gatk compatibility', type=str, required=False, default="OA")
    parser.add_argument('--version', action='version', version='%(prog)s '+version)


    args = parser.parse_args()
    #main subroutine
    count,trimmed_count_first_of_pair,trimmed_count_second_of_pair = revert_hardtrim_adapters(args)
    
    logging.info("Ended script.")
    logging.info(" Run time " + str(timedelta(seconds=time.time()- start_time)) + ".")
    logging.info(" Records processsed "+str(count) + ".")
    logging.info(" Trimmed first of pair "+ str(trimmed_count_first_of_pair) + ".")
    logging.info(" Trimmed second of pair "+ str(trimmed_count_second_of_pair) + ".")



def revert_hardtrim_adapters(args):
    return(0,0,0)