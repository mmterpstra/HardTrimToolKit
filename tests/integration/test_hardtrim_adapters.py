import pysam
import sys
import os
from hardtrimtoolkit import trim_aligned_segment_by_query_offset
from hardtrimtoolkit.hardtrimtoolkit import hardtrim_adapters
from argparse import Namespace

#@pytest.mark.integration
def test_trim_aligned_segment_by_query_offset(caplog):
    caplog.set_level("DEBUG")
    #this is a space for testing the results of trim_aligned_segment_by_query_offset()
    def get_test_read_1():
        header_dict = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                {'LN': 1584, 'SN': 'chr2'}],
            'RG':[{"ID":'Sample', 'BC':'NNNNNN'}],
        }
        header=pysam.AlignmentHeader.from_dict(header_dict)
    
        # basic test
        read = pysam.AlignedSegment(header)
        read.query_name = "read_28833_29006_6945"
        read.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        read.flag = 99
        read.reference_id = 0
        read.reference_start = 32 #this is zero based and sam is 1 based
        read.mapping_quality = 20
        read.cigar = ((pysam.CMATCH,10), (pysam.CDEL,1), (pysam.CMATCH,25))
        read.next_reference_id = 0
        read.next_reference_start=199
        read.template_length=167
        read.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
        read.tags = (("NM", 1),
            ("RG", "Sample"),("XT", 8))
        return read 
    
    read= get_test_read_1()
    trim_aligned_segment_by_query_offset(read, 8, 'XO')
    assert("This is the read being checked:" + read.to_string())
    assert(read.to_string() == "\t".join([
        'read_28833_29006_6945',
        '99',
        'chr1',
        '33',
        '20',
        '10M1D17M8H',
        '=',
        '200',
        '167',
        'AGCTTAGCTAGCTACCTATATCTTGGT',
        '<<<<<<<<<<<<<<<<<<<<<:<9/,&',
        'NM:i:1',
        'RG:Z:Sample',
        'XT:i:8',
        'XO:Z:chr1,32,10M1D25M,AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG,<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<;']))


    def get_test_read_2():
        header_dict = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                {'LN': 1584, 'SN': 'chr2'}],
            'RG':[{"ID":'Sample', 'BC':'NNNNNN'}],
        }
        header=pysam.AlignmentHeader.from_dict(header_dict)
    
        # basic test
        read = pysam.AlignedSegment(header)
        read.query_name = "read_28833_29006_6945"
        read.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        read.flag = 16
        read.reference_id = 0
        read.reference_start = 32 #this is zero based and sam is 1 based
        read.mapping_quality = 20
        read.cigar = ((pysam.CMATCH,10), (pysam.CDEL,1), (pysam.CMATCH,25))
        read.next_reference_id = 0
        read.next_reference_start=199
        read.template_length=167
        read.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
        read.tags = (("NM", 1),
            ("RG", "Sample"),
            ("XT", 8))
        return read 

    read= get_test_read_2()
    trim_aligned_segment_by_query_offset(read, 8, 'XO')
    assert("This is the read being checked:" + read.to_string())
    assert(read.to_string() == "\t".join([
        'read_28833_29006_6945',
        '16',
        'chr1',
        '41',
        '20',
        '8H2M1D25M',
        '=',
        '200',
        '167',
        'TAGCTACCTATATCTTGGTCTTGGCCG',
        '<<<<<<<<<<<<<:<9/,&,22;;<<<',
        'NM:i:1',
        'RG:Z:Sample',
        'XT:i:8',
        'XO:Z:chr1,32,10M1D25M,AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG,<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<;']))

    read = get_test_read_2()
    assert(read.to_string() == 
        "\t".join(
            [
                'read_28833_29006_6945',
                '16',
                'chr1',
                '33',#32(orignial zero based pos)+1(make it one based) +9 (trim) +1(Deletion) = 45?
                '20',
                '10M1D25M',
                '=',
                '200',
                '167',
                'AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG',
                '<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<',
                'NM:i:1',
                'RG:Z:Sample',
                'XT:i:8'
            ]
        )
    )
    
    trim_aligned_segment_by_query_offset(read, 10, 'XO')
    assert("This is the read being checked:" + read.to_string())
    assert(read.to_string() == "\t".join([
        'read_28833_29006_6945',
        '16',
        'chr1',
        '44',#32(original pos)+1 (zero based conversion) +10 (trim)   +1(Deletion) = 44?
        '20',
        '10H25M',
        '=',
        '200',
        '167',
        'GCTACCTATATCTTGGTCTTGGCCG',
        '<<<<<<<<<<<:<9/,&,22;;<<<',
        'NM:i:1',
        'RG:Z:Sample',
        'XT:i:8',
        'XO:Z:chr1,32,10M1D25M,AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG,<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<;']))
    #assert any("hardtrim" in r.message for r in caplog.records)


#@pytest.mark.integration
def test_hardtrim_cli(caplog):
    caplog.set_level("DEBUG")
    
    def write_testbam(bam):
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                    {'LN': 1584, 'SN': 'chr2'}],
            'RG': [{"ID":'Sample', 'BC':'NNNNNN'}]
        }
        with pysam.AlignmentFile(bam, "wb", header=header) as outf:
            a = pysam.AlignedSegment()
            a.query_name = "read_28833_29006_6945"
            a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            a.flag = 99
            a.reference_id = 0
            a.reference_start = 32
            a.mapping_quality = 20
            a.cigar = ((0,10), (2,1), (0,25))
            a.next_reference_id = 0
            a.next_reference_start=199
            a.template_length=167
            a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
            a.tags = (("NM", 1),
                ("RG", "Sample"),
                ("XT",8)
            )
            outf.write(a)
            b = pysam.AlignedSegment()
            b.query_name = "read_28833_29006_6945"
            b.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            b.flag = 16
            b.reference_id = 0
            b.reference_start = 32 #this is zero based and sam is 1 based
            b.mapping_quality = 20
            b.cigar = ((pysam.CMATCH,10), (pysam.CDEL,1), (pysam.CMATCH,25))
            b.next_reference_id = 0
            b.next_reference_start=199
            b.template_length=167
            b.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
            b.tags = (("NM", 1),
                ("RG", "Sample"),
                ("XT", 8))
            outf.write(b)

    testinputbam="tests/integration/test_hardtrim_adapters.bam"
    testoutbam="tests/integration/test_hardtrim_adapters_output.bam"
    write_testbam(testinputbam)
    args = Namespace(
        input=testinputbam,
        output=testoutbam,
        ubam_out=False,
        ubam_in=False,
        trimming_tag="XT",
        original_alignment_tag="XO",
        append_hardtrim=True,
    )
    hardtrim_adapters(args)

    test_result = pysam.samtools.view(testoutbam, split_lines=True)

    expected_result = [
        'read_28833_29006_6945\t99\tchr1\t33\t20\t10M1D17M8H\t=\t200\t167\tAGCTTAGCTAGCTACCTATATCTTGGT\t<<<<<<<<<<<<<<<<<<<<<:<9/,&\tNM:i:1\tRG:Z:Sample\tXO:Z:chr1,32,10M1D25M,AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG,<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<;',
        #minor nitpick the logic regarding the two reads doesn't make sense from a sequencing perspective. it just shows two reads being trimmed.
        'read_28833_29006_6945\t16\tchr1\t41\t20\t8H2M1D25M\t=\t200\t167\tTAGCTACCTATATCTTGGTCTTGGCCG\t<<<<<<<<<<<<<:<9/,&,22;;<<<\tNM:i:1\tRG:Z:Sample\tXO:Z:chr1,32,10M1D25M,AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG,<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<;'
    ]

    assert len(test_result) == len(expected_result)

    for res, exp in zip(test_result, expected_result):
        assert res==exp

    os.remove(testinputbam)
    os.remove(testoutbam)