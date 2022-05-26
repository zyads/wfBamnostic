"""
Assemble and sort some COVID reads...
"""

# import subprocess
from pathlib import Path
from typing import List, Optional, Union

import bamnostic as bs
from latch import small_task, workflow
from latch.types import LatchDir, LatchFile

# Bamnostic is meant to be a reduced drop-in replacement for pysam. As such it has much the same API as pysam
# with regard to BAM-related operations.
# Note: the pileup() method is not supported at this time.


# User input BAM file.
@small_task
def bam_tsk(
    bam_file: Optional[LatchFile],
    output_dir: Optional[LatchDir],
    reads_n: Optional[int],
    complex_read_name: Optional[str],
    chr: Optional[str],
    start: int,
    end: int,
) -> LatchFile:
    bam = bs.AlignmentFile(bam_file, "rb")  # filepath_or_object, mode='rb

    header = bam.header  # Data validation
    head = bam.head(n=reads_n)  # Returns the first n alignments

    for complex_read in bam:
        if complex_read.read_name == complex_read_name:
            returnableread = complex_read
            break

    # serial access
    for i, read in enumerate(bam.fetch(chr, start, end)):
        if i >= 3:
            break
        returnableread2 = read

    # make output .out file with output variables in it
    out_file = Path(output_dir) / f"{complex_read_name}.out"
    with open(out_file, "w") as f:
        # write header, head, read, and serial access to output file
        f.write(f"Header: {header}\n")
        f.write(f"Head: {head}\n")
        f.write(f"Read: {returnableread}\n")
        f.write(f"Serial Access: {returnableread2}\n")

    return LatchFile(str(out_file), "latch:///" + str(out_file))


@workflow
def bam_wf(
    bam_file: Optional[LatchFile] = None,
    output_dir: Optional[LatchDir] = None,
    reads_n: Optional[int] = 20,
    complex_read_name: Optional[str] = None,
    chr: Optional[str] = None,
    start: int = 0,
    end: int = 1,
) -> LatchFile:
    """

     __metadata__:
             display_name: BAM file parser and random access tools
             author:
                 name:
                 email:
                 github:
             repository: https://github.com/mills-lab/bamnostic
             license:
                 id: MIT
    Args:

             bamfile:

               __metadata__:
                 display_name: In

             output_dir:


               __metadata__:
                 display_name: Out

             reads_n:

               __metadata__:
                 display_name: Reads

             complex_read_name:

               __metadata__:
                 display_name: Complex Read Name

             chr:

               __metadata__:
                 display_name: Chromosome

             start:

               __metadata__:
                 display_name: Start

             end:

               __metadata__:
                 display_name: End


    """

    return bam_tsk(
        bam_file=bam_file,
        output_dir=output_dir,
        reads_n=reads_n,
        complex_read_name=complex_read_name,
        chr=chr,
        start=start,
        end=end,
    )


# Keep wrapping parameters
# Print log to print error messages if chromosome not found

# Make new task for different nodes

# CIGAR string:
# M: match
# I: insertion
# D: deletion
# N: skipped
# S: soft clipping
# H: hard clipping
# P: padding
# =: sequence match
# X: sequence mismatch
#
# QUAL string:
# ASCII representation of Phred quality score
#
# TAGS:
#     XS: The number of soft clipped bases
#     XH: The number of hard clipped bases
#     XQ: The number of bases with quality below the minimum quality threshold
#     XA: The number of aligned bases
#     XC: The number of aligned bases that are clipped
#     XG: The number of aligned bases that are gap
#     XN: The number of aligned bases that are N
#     XO: The number of aligned bases that are not primary
#     XQ: The number of aligned bases that are not primary and not aligned
#     XR: The number of aligned bases that are not primary and not aligned and not clipped
#     XU: The number of aligned bases that are not primary and not aligned and not clipped and not gap
#     XW: The number of aligned bases that are not primary and not aligned and not clipped and not gap and not N
#     XT: The number of aligned bases that are not primary and not aligned and not clipped and not gap and not N and not quality below the minimum quality threshold
#     X0: The number of aligned bases that are not primary and not aligned and not clipped and not gap and not N and not quality below the minimum quality threshold and not aligned
#     X1: The number of aligned bases that are not primary and not aligned and not clipped and not gap and not N and not quality below the minimum quality threshold and not aligned and not clipped

# Remember that BAM file is 1 indexed so must need to account for that???

# Write the reads to a file.
# with open(output_dir / 'complex_read.txt', 'w') as f:
# returnableread.save(output_dir / (complex_read_name + ".bam"))
# bam.close() #Remember to close the file, otherwise it will be locked?
