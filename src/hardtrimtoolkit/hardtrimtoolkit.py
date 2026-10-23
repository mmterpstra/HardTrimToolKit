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
    logging.basicConfig(level=logging.INFO, format='# %(asctime)s # %(filename)s.%(funcName)s # %(levelname)s # %(message)s')
    logging.info("Started "+sys.argv[0]+" version "+version+".")
    logging.info("WARNING experimental version also outside of the fields position/cigar/seq and qual it does not update correctly.")
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description='Takes a readname sorted bam file and uses the adapter trimming tag to change the soft trims into hard trims.' \
        ' Stores original optionally in the XO tag for reversion by the reverthardtrim tool. Should be used preferably' \
        ' before sorting the alignment because coordinate sorting is made undone. ' \
        'Also it doesnt update all the info TLEN so preferably you use the qname sorted file to fix the data with picard tools.',
        epilog='Example useage hardtrimadapters --input input.bam --output output.bam --fastq qname_R1.fastq.gz --fastq qname_R2.fastq.gz')
    parser.add_argument('--output',help='Output bam. Should support "-" as /dev/stdout', required=True)
    parser.add_argument('--input', help='Input bam. Should support "-" as /dev/stdin', required=True )
    parser.add_argument('--ubam-out',help='Write uncompressed output bam. Helpful for speeding up on stream processing',
         action='store_true', required=False, default=False )
    parser.add_argument('--ubam-in',help='Read uncompressed input bam.',  action='store_true', required=False, default=False )
    parser.add_argument('--trimming-tag',help='Tag guide in trimmings. Default is XT for gatk compatibility', type=str, required=False, default="XT")
    parser.add_argument('--original-alignment-tag',help='Tag to move original alignment to. Default is XO. Not the same format as OA.' \
        ' format as : "chr,pos,cigar,seq,qual;"', type=str, required=False, default="XO")
    parser.add_argument('--omit-hardtrim',help='Option to remove in the hardtrim in the cigar.e.g. ' \
        'trim 101M + XT 10 -> 91M instead of 91M10H', action='store_true', required=False, default=False )
    parser.add_argument('--version', action='version', version='%(prog)s '+version)


    args = parser.parse_args()
    #main subroutine
    count,trimmed_count_first_of_pair,trimmed_count_second_of_pair = hardtrim_adapters(args)
    
    logging.info("Ended script.")
    logging.info("Run time " + str(timedelta(seconds=time.time()- start_time)) + ".")
    logging.info("Records processsed "+str(count) + ".")
    logging.info("Trimmed first of pair "+ str(trimmed_count_first_of_pair) + ".")
    logging.info("Trimmed second of pair "+ str(trimmed_count_second_of_pair) + ".")

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
    count = 0
    trimmed_count_first_of_pair = 0
    trimmed_count_second_of_pair = 0
    readmode = 'rb'
    
    #assume aligned unsorted file
    inbam = pysam.AlignmentFile(args.input, readmode, check_sq=True, require_index=False)
    outheader = inbam.header.to_dict()
    outheader['PG'] = [{'ID': '1','PN':'hardtrimadapters','DS':'Hard removal of trimmed adapter sequence in aligned bams','VN':'version'+version,'CL':' '.join(sys.argv)}]
    #outheader['HD']['SO'] processing changes based on this value coordinate, queryname or unsorted
    writemode = 'wb'
    if args.ubam_out:
        writemode = 'wbu'
    outfile = pysam.AlignmentFile(args.output, writemode, header=outheader, check_sq=True, require_index=False)
    saveoriginal = True
    recordbuffer = []
    
    def validate_sam_tag(tag:str):
        if not(len(tag)==2):
            raise ValueError("Tag length '"+tag+"' is longer or shorter then 2 chars! This breaks sam spec.")

    validate_sam_tag(args.trimming_tag)
    validate_sam_tag(args.original_alignment_tag)

    def add_hardtrim_bamrecords(
        bamrecords: list[pysam.AlignedSegment],
        trimming_tag: str='XT',
        original_alignment_tag: str='XO',
        append_hardtrim: bool=False
    ):
        for bamrecord in bamrecords:
            #update tlen? or let it be fixed in postprocessing?
            if(bamrecord.has_tag("XT")):
                #print(file=sys.stderr, bamrecord)
                #Does trim on XT field in bamrecord
                trim_aligned_segment_by_query_offset(read=bamrecord, trim=bamrecord.get_tag(trimming_tag))

    for recordno, bamrecord in enumerate(inbam):
        if recordno % 1000000 == 0: 
            logging.info(" Records processsed "+str(recordno) + ".")
        if len(recordbuffer) == 0 :
            if not(bamrecord.is_paired):
                #process
                recordbuffer.append(bamrecord)
                add_hardtrim_bamrecords(
                    bamrecords = recordbuffer,
                    trimming_tag = args.trimming_tag,
                    original_alignment_tag=args.original_alignment_tag,
                    append_hardtrim =  not(args.omit_hardtrim))
                for record in recordbuffer:
                    record.set_tag(args.trimming_tag,None)
                    outfile.write(record)
                recordbuffer = []

            else:
                recordbuffer.append(bamrecord)
            
        elif recordbuffer[0].query_name != bamrecord.query_name:
            #process paired reads assuming not paired if possible
            #process
            add_hardtrim_bamrecords(
                bamrecords = recordbuffer,
                trimming_tag = args.trimming_tag,
                original_alignment_tag=args.original_alignment_tag,
                append_hardtrim = not(args.omit_hardtrim))
            for record in recordbuffer:
                record.set_tag(args.trimming_tag,None)
                if record.is_read1:
                    trimmed_count_first_of_pair=trimmed_count_first_of_pair+1
                elif record.is_read2:
                    trimmed_count_second_of_pair=trimmed_count_second_of_pair+1

                outfile.write(record)
            #push new
            recordbuffer = [bamrecord]
        else:
            recordbuffer.append(bamrecord)
    if len(recordbuffer) != 0 :
        add_hardtrim_bamrecords(
            bamrecords = recordbuffer,
            trimming_tag = args.trimming_tag,
            original_alignment_tag=args.original_alignment_tag,
            append_hardtrim =  not(args.omit_hardtrim))
        for record in recordbuffer:
            record.set_tag(args.trimming_tag,None)
            if record.is_read1:
                trimmed_count_first_of_pair=trimmed_count_first_of_pair+1
            elif record.is_read2:
                trimmed_count_second_of_pair=trimmed_count_second_of_pair+1
            outfile.write(record)
    count=recordno
    inbam.close()
    outfile.close()
    return(count,trimmed_count_first_of_pair,trimmed_count_second_of_pair)

def revert_hardtrim_adapters_cli():
    start_time = time.time()
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='# %(asctime)s ##INFO## %(message)s')
    logging.info("Started "+sys.argv[0]+" version "+version+".")
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description=' Uses the XO tag to revert mappings to pre-hardtrimmed bams.',
        epilog='Example useage revert_hardtrim_adapters --input input.bam --output output.bam ')
    parser.add_argument('--output',help='Output bam. Should support "-" as /dev/stdout', required=True)
    parser.add_argument('--input', help='Input bam assumes assumes sorted by readname. Should support "-" as /dev/stdin', required=True )
    parser.add_argument('--ubam-out',help='Write uncompressed output bam. Helpful for speeding up on stream processing',  type=bool, required=False, default=False )
    parser.add_argument('--ubam-in',help='Read uncompressed input bam.', type=bool, required=False, default=False )
    parser.add_argument('--trimming-tag',help='Tag to add in trimmings. Default is XT for gatk compatibility', type=str, required=False, default="XT")
    parser.add_argument('--original-alignment-tag',help='Tag to move original alignment from. Default is XO', type=str, required=False, default="XO")
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
    #pseudocode
    #handle switches
    #iterate bam file:
    #   update tags ao with the original alignment results
    #   hardtrim the read and update position and sam info to match
    readmode = 'rb'
    inbam = pysam.AlignmentFile(args.input, readmode, check_sq=False, require_index=False)
    outheader = inbam.header.to_dict()
    if not 'PG' in outheader: 
        outheader['PG'] = [{'ID': '1','PN':'reverthardtrimadapters','DS':'External trimming tool integration','VN':'version'+version,'CL':' '.join(sys.argv)}]
    else:
        outheader['PG'].append({'ID': '1','PN':'reverthardtrimadapters','DS':'Reverts your bam to input','VN':'version'+version,'CL':' '.join(sys.argv)})

    def validate_sam_tag(tag:str):
        if not(len(tag)==2):
            raise ValueError("Tag length '"+tag+"' is longer or shorter then 2 chars! This breaks sam spec.")
    
    validate_sam_tag(args.trimming_tag)
    validate_sam_tag(args.original_alignment_tag)
    trimmed_count_first_of_pair = 0
    trimmed_count_second_of_pair = 0
    writemode = 'wb'
    if args.ubam_out:
        writemode = 'wbu'
    outfile = pysam.AlignmentFile(args.output, writemode, header=outheader, check_sq=False, require_index=False)
    for lineno, bamrecord in enumerate(inbam):
        if lineno % 1000000 == 0: 
            logging.info(" Records processsed "+str(lineno) + ".")
        print(pprint(bamrecord))
        revert_hardtrim(bamrecord,args.trimming_tag,args.original_alignment_tag)
        exit(1)
        outfile.write(bamrecord)
        count = count + 1
    outfile.close()
    return count,trimmed_count_first_of_pair,trimmed_count_second_of_pair

def trim_aligned_segment_by_ref_offset(read: pysam.AlignedSegment, trim: int, trimmed_seq_tag: str) -> pysam.AlignedSegment:
    """
    WIP not needed atm
    Re‑implementation of the Perl TrimSam() logic using pysam.AlignedSegment.
    This trims `trim` bases from either the 5′ or 3′ end depending on orientation and offset. 
    The offset is on reference space not on query space but on reference space. 
    Reference space is calced by intersection of [bcs]am record with bed file.  

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to modify
    trim : int
        Number of aligned bases to trim (overlap amount)

    Returns
    -------
    pysam.AlignedSegment
        Modified read with updated CIGAR, sequence, qualities, and position.
    """

    # ----------------------------------------------------------------------
    # 1. Determine trimming direction
    # ----------------------------------------------------------------------
    # Perl logic:
    #   - If reverse + unpaired  → trim 5'
    #   - If reverse & pair==1   → trim 5'
    #   - If forward & pair==2   → trim 5'
    #   Otherwise trim 3'
    # ----------------------------------------------------------------------

    is_reverse = read.is_reverse
    is_paired = read.is_paired
    pe_end = None
    if is_paired:
        if read.is_read1:
            pe_end = 1
        elif read.is_read2:
            pe_end = 2

    process_reverse = False
    if (
        (is_reverse and not is_paired) or
        (is_reverse and is_paired and pe_end == 1) or
        (not is_reverse and is_paired and pe_end == 2)
    ):
        # Trim 5′ end (left side in reference direction)
        process_reverse = True

    # ----------------------------------------------------------------------
    # 2. Extract original cigar, seq and qual
    # ----------------------------------------------------------------------
    cigar = read.cigartuples            # list of (op, length)
    seq = read.query_sequence
    qual = read.query_qualities

    original_record = ','.join(
            [
                str(read.reference_name),
                str(read.reference_start),
                read.cigarstring,
                ''.join(read.query_sequence),
                pysam.array_to_qualitystring(read.query_qualities)+';'
            ]
    )
    # op codes from pysam:
    # 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7=X,8==
    # We will reconstruct a new cigar list.
    new_cigar = []
    hardclip_count = trim

    # Work on a mutable copy
    cigar_work = cigar.copy()

    # ----------------------------------------------------------------------
    # 3. Walk through cigar ops from correct side until trim exhausted
    # ----------------------------------------------------------------------
    def take(op_tuple):
        """Pop from left or right depending on trimming direction."""
        if process_reverse:
            return cigar_work.pop(0)
        else:
            return cigar_work.pop()

    # Helpers for sequence clipping
    def clip_left(seq:str, qual:str,n:int):
        seq = seq[n:]
        qual = qual[n:]
        return seq, qual

    def clip_right(seq:str, qual:str,n:int):
        seq = seq[:-n]
        qual = qual[:-n]
        return seq, qual
        
    while cigar_work and trim >= 0:
        op, length = take(cigar_work)
        #| Operator   | Operation                    | Int   | Comment                     | 
        #| -------    | ----------                   | ---   | -------                     | 
        #| M          | pysam.CMATCH       | 0     | Seq match                   | 
        #| I          | pysam.CINS         | 1     | Insertion                   | 
        #| D          | pysam.CDEL         | 2     | Deletion                    |
        #| N          | pysam.CREF_SKIP    | 3     | Reference skip(splicing)    |
        #| S          | pysam.CSOFT_CLIP   | 4     | Not aligned seq soft clip   | 
        #| H          | pysam.CHARD_CLIP   | 5     | Not aligned seq hard clip   | 
        #| P          | pysam.CPAD         | 6     | Silent indel for padding    | 
        #| =          | pysam.CEQUAL       | 7     | Seq match with seq == ref   |
        #| X          | pysam.CDIFF        | 8     | Seq match with seq != ref   |
        #| B          | pysam.CBACK        | 9     | backward compat ?           |

        if op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF):  # M, =, X → aligned bases
            if trim >= length:
                # consume whole block
                if process_reverse:
                    (seq, qual)=clip_left(seq, qual, length)
                    read.reference_start += length
                else:
                    (seq, qual)=clip_right(seq, qual, length)
                trim -= length
                # record zero-length op
                amount_new = 0
                op_new = op
            else:
                # partial trimming inside this block
                amount_new = length - trim
                op_new = op
                if process_reverse: 
                    (seq, qual)=clip_left(seq, qual, trim)
                    read.reference_start += trim
                else:
                    (seq, qual)=clip_right(seq, qual, trim)
                trim = -1  # done
        elif op in (pysam.CINS, pysam.CSOFT_CLIP):  # I or S → affects read seq
            hardclip_count += length
            amount_new = 0
            op_new = op
            if process_reverse:
                (seq, qual)=clip_left(seq, qual, length)
            else:
                (seq, qual)=clip_right(seq, qual, length)
        elif op in (pysam.CDEL, pysam.CREF_SKIP):  # D or N → reference-only
            amount_new = 0
            op_new = op
            if process_reverse:
                read.reference_start += length
            if trim >= length:
                trim -= length
            else:
                trim -= length  # allow negative
        elif op in (5,6):
            raise ValueError("Hardclip/Pad in original CIGAR not supported")

        # Append reduced op to new cigar
        if process_reverse:
            new_cigar.append((op_new, amount_new))
        else:
            new_cigar.insert(0, (op_new, amount_new))

    # ----------------------------------------------------------------------
    # 4. Append remaining cigar ops
    # ----------------------------------------------------------------------
    if process_reverse:
        new_cigar.extend(cigar_work)
        new_cigar.insert(0, (pysam.CHARD_CLIP, hardclip_count))  # leading H
    else:
        new_cigar = cigar_work + new_cigar
        new_cigar.append((pysam.CHARD_CLIP, hardclip_count))     # trailing H

    # ----------------------------------------------------------------------
    # 5. Clean zero‑length ops and rebuild
    # ----------------------------------------------------------------------
    final_cigar = [(op, ln) for (op, ln) in new_cigar if ln > 0]

    read.cigartuples = final_cigar
    read.query_sequence = seq
    read.query_qualities = qual

    return read


def trim_aligned_segment_by_query_offset(read: pysam.AlignedSegment, trim: int, original_record_tag: str="XO", append_hardtrim: bool=True) -> pysam.AlignedSegment:
    """
    Re‑implementation of the perl trim by bed script. 
    ( github.com/mmterpstra/pipeline-util/src/trimByBed.pl ) TrimSam() logic using pysam.AlignedSegment.
    This trims `trim` bases from either the 5′ or 3′ end depending on orientation and offset. 
    The offset is on query space.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to modify
    trim : int
        Number of aligned bases to trim (overlap amount)

    Returns
    -------
    pysam.AlignedSegment
        Modified read with updated CIGAR, sequence, qualities, and position.
    """

    # ----------------------------------------------------------------------
    # 1. Determine trimming direction
    # ----------------------------------------------------------------------
    # Perl logic:
    #   - If reverse + unpaired  → trim 5'
    #   - If reverse & pair==1   → trim 5'
    #   - If forward & pair==2   → trim 5'
    #   Otherwise trim 3'
    # ----------------------------------------------------------------------

    is_reverse = read.is_reverse
    is_paired = read.is_paired
    pe_end = None
    if is_paired:
        if read.is_read1:
            pe_end = 1
        elif read.is_read2:
            pe_end = 2

    process_reverse = False
    if (
        (is_reverse and not is_paired) or
        (is_reverse and is_paired and pe_end == 1) or
        (not is_reverse and is_paired and pe_end == 2)
    ):
        # Trim 5′ end (left side in reference direction)
        process_reverse = True

    # ----------------------------------------------------------------------
    # 2. Extract original cigar, seq and qual
    # ----------------------------------------------------------------------
    cigar = read.cigartuples            # list of (op, length)
    seq = read.query_sequence
    qual = read.query_qualities
    position = read.reference_start

    original_record = ','.join(
        [
            str(read.reference_name),
            str(read.reference_start),
            read.cigarstring,
            ''.join(read.query_sequence),
            pysam.array_to_qualitystring(read.query_qualities)+';'
        ]
    )
    # op codes from pysam:
    # 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7=X,8==
    # We will reconstruct a new cigar list.
    new_cigar = []
    hardclip_count = trim

    # Work on a mutable copy
    cigar_work = cigar.copy()

    # ----------------------------------------------------------------------
    # 3. Walk through cigar ops from correct side until trim exhausted
    # ----------------------------------------------------------------------
    def take(op_tuple):
        """Pop from left or right depending on trimming direction."""
        if process_reverse:
            return cigar_work.pop(0)
        else:
            return cigar_work.pop()

    # Helpers for sequence trimming
    def clip_left(seq:str, qual:str,n:int):
        seq = seq[n:]
        qual = qual[n:]
        return seq, qual

    def clip_right(seq:str, qual:str,n:int):
        seq = seq[:-n]
        qual = qual[:-n]
        return seq, qual
        
    while cigar_work and trim >= 0:
        op, length = take(cigar_work)
        #| Operator   | Operation                    | Int   | Comment                     | 
        #| -------    | ----------                   | ---   | -------                     | 
        #| M          | pysam.CMATCH       | 0     | Seq match                   | 
        #| I          | pysam.CINS         | 1     | Insertion                   | 
        #| D          | pysam.CDEL         | 2     | Deletion                    |
        #| N          | pysam.CREF_SKIP    | 3     | Reference skip(splicing)    |
        #| S          | pysam.CSOFT_CLIP   | 4     | Not aligned seq soft clip   | 
        #| H          | pysam.CHARD_CLIP   | 5     | Not aligned seq hard clip   | 
        #| P          | pysam.CPAD         | 6     | Silent indel for padding    | 
        #| =          | pysam.CEQUAL       | 7     | Seq match with seq == ref   |
        #| X          | pysam.CDIFF        | 8     | Seq match with seq != ref   |
        #| B          | pysam.CBACK        | 9     | backward compat ?           |

        if op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF):  # M, =, X → aligned bases
            if trim >= length:
                # consume whole block
                if process_reverse:
                    (seq, qual)=clip_left(seq, qual, length)
                    read.reference_start += length
                else:
                    (seq, qual)=clip_right(seq, qual, length)
                trim -= length
                # record zero-length op
                amount_new = 0
                op_new = op
            else:
                # partial trimming inside this block
                amount_new = length - trim
                op_new = op
                if process_reverse: 
                    (seq, qual)=clip_left(seq, qual, trim)
                    read.reference_start += trim
                else:
                    (seq, qual)=clip_right(seq, qual, trim)
                trim = -1  # done
        elif op in (pysam.CINS, pysam.CSOFT_CLIP):  # I or S → affects read seq
            hardclip_count += length
            amount_new = 0
            op_new = op
            if process_reverse:
                (seq, qual)=clip_left(seq, qual, length)
            else:
                (seq, qual)=clip_right(seq, qual, length)
            if trim >= length:
                trim -= length
            else:
                trim -= length  # allow negative
        elif op in (pysam.CDEL, pysam.CREF_SKIP):  # D or N → reference-only
            amount_new = 0
            op_new = op
            if process_reverse:
                read.reference_start += length
            #This doesnt matter for sequence space trimming
            #if trim >= length:
            #    trim -= length
            #else:
            #    trim -= length  # allow negative
        elif op in (pysam.CBACK):
            raise ValueError("Hardclip/Pad in original CIGAR not supported")

        # Append reduced op to new cigar
        if process_reverse:
            new_cigar.append((op_new, amount_new))
        else:
            new_cigar.insert(0, (op_new, amount_new))

    # ----------------------------------------------------------------------
    # 4. Append remaining cigar ops
    # ----------------------------------------------------------------------
    if append_hardtrim:
        if process_reverse:
            new_cigar.extend(cigar_work)
            new_cigar.insert(0, (pysam.CHARD_CLIP, hardclip_count))  # leading H
        else:
            new_cigar = cigar_work + new_cigar
            new_cigar.append((pysam.CHARD_CLIP, hardclip_count))     # trailing H

    # ----------------------------------------------------------------------
    # 5. Clean zero‑length ops and rebuild
    # ----------------------------------------------------------------------
    final_cigar = [(op, ln) for (op, ln) in new_cigar if ln > 0]

    read.cigartuples = final_cigar
    read.query_sequence = seq
    read.query_qualities = qual
    read.set_tag(original_record_tag, original_record, value_type="Z")
    
    return read

def revert_hardtrim(
    read: pysam.AlignedSegment,
    trimming_tag: str="XT",
    original_record_tag: str="XO"
):
    """Restore the original alignment stored in `original_record_tag`.

    Expected XO format:
        query_name,reference_start,cigar,sequence,qualities;

    note the , split and the ending ; for record safety.
    """
   

    if read.has_tag(original_record_tag):

        (
            original_record_name,
            original_reference_start,
            original_cigar_string,
            original_query_sequence,
            original_qualities_and_more
        ) = read.get_tag(original_record_tag).split(',',maxsplit=4)
        original_qualities=original_qualities_and_more[:len(original_query_sequence)]
        original_more = original_qualities_and_more[len(original_query_sequence):]
        
        #checks for format sanity
        if not(original_record_name == read.query_name):
            print(read,sys.stderr)
            raise ValueError("Record names should be identical something broke or is incorrectly handled. Affliced record name :'"+original_record_name+"'")
        if not(original_more.startswith(";")):
            raise ValueError("Quality string is not followed by ';' looks somebody broke the spec")
        if len(original_more) > 1 :
            #while not 
            raise ValueError("Processing multiple original_record_tags (tag'"+original_record_tag+"',info:'"+read.get_tag(original_record_tag)+"') not supported.")
        
        trimming_length = len(original_query_sequence) - read.query_length
        
        read.reference_start = int(original_reference_start)
        read.cigarstring = original_cigar_string
        read.query_sequence = original_query_sequence
        read.query_qualities = pysam.qualitystring_to_array(original_qualities)
        
        read.set_tag(trimming_tag,trimming_length)
        if read.has_tag(original_record_tag):
            read.set_tag(original_record_tag,None)

    return read