"""
Assemble and sort some COVID reads...
"""

# import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir

from typing import List, Union, Optional
import bamnostic as bs


# Bamnostic is meant to be a reduced drop-in replacement for pysam. As such it has much the same API as pysam
# with regard to BAM-related operations.
# Note: the pileup() method is not supported at this time.


#User input BAM file.
@small_task
def bam_tsk(
    bam_file: Optional[LatchFile],
    output_dir: Optional[LatchDir],

    reads_n: Optional[int],

    complex_read_name : str,

    chr: str,
    start: int,
    end: int,

) -> LatchFile: 
    bam = bs.AlignmentFile(bam_file, 'rb') # filepath_or_object, mode='rb

    header = bam.header #Data validation
    head = bam.head(n=reads_n) #Returns the first n alignments

    for complex_read in bam:
        if complex_read.read_name == complex_read_name:
            returnableread = complex_read
            break

    #serial access
    for i, read in enumerate(bam.fetch(chr, start, end)):
        if i >= 3:
            break
        returnableread2 = read

    #make output .out file with output variables in it
    out_file = Path(output_dir) / f"{complex_read_name}.out"
    with open(out_file, 'w') as f:
        # write header, head, read, and serial access to output file
        f.write(f"Header: {header}\n")
        f.write(f"Head: {head}\n")
        f.write(f"Read: {returnableread}\n")
        f.write(f"Serial Access: {returnableread2}\n")

    return LatchFile(str(out_file), "latch:///" + str(out_file))


# @small_task
# def guide_counter_task(
#     reads: LatchFile,
#     output_name: str,
# ) -> (LatchFile, LatchFile, LatchFile):

#     _guide_counter_cmd = [
#         "guide-counter",
#         "count",
#         "--input",
#         str(Path(reads).resolve()),
#         "--library",
#         "/root/brunello.csv",
#         "--output",
#         output_name,
#     ]

#     counts = Path(f"/root/{output_name}.counts.txt")
#     extended_counts = Path(f"/root/{output_name}.extended-counts.txt")
#     stats = Path(f"/root/{output_name}.stats.txt")

#     subprocess.run(_guide_counter_cmd)

#     return (
#         LatchFile(str(counts), f"latch://{counts}"),
#         LatchFile(str(extended_counts), f"latch://{extended_counts}"),
#         LatchFile(str(stats), f"latch://{stats}"),
#     )



@workflow
def bam_wf(
    bam_file: Optional[LatchFile] = None,
    output_dir: Optional[LatchDir] = None,

    reads_n: Optional[int] = 20,
    complex_read_name : str = None,

    chr: str = None,
    start: int = 0,
    end: int = 1,

) -> LatchFile:
    """
     __metadata__:
    #         display_name: BAM file parser and random access tools
    #         author:
    #             name:
    #             email:
    #             github:
    #         repository: https://github.com/mills-lab/bamnostic
    #         license:
    #             id: MIT
    Args:

    #         bamfile:

    #           __metadata__:
    #             display_name: In

    #         output_dir:


    #           __metadata__:
    #             display_name: Out

    #         reads_n:

    #           __metadata__:
    #             display_name: Reads

    #         complex_read_name:

    #           __metadata__:
    #             display_name: Complex Read Name

    #         chr:
    
    #           __metadata__:
    #             display_name: Chromosome

    #         start:

    #           __metadata__:
    #             display_name: Start

    #         end:

    #           __metadata__:
    #             display_name: End


    """

    return bam_tsk(
        bam_file = bam_file,
        output_dir = output_dir,

        reads_n = reads_n,
        complex_read_name = complex_read_name,

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



#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

# @small_task
# def bactopia_tsk(
#     fastq_one: Optional[LatchFile],
#     fastq_two: Optional[LatchFile],
#     input_dir: Optional[LatchDir],
#     output_dir: LatchDir,
#     sample_name: List[Union[str, int]],
#     genome_size: Optional[int],
#     species: Optional[Species],
#     species_genome_size: Optional[SpeciesGenomeSize],
#     ask_merlin: bool = False,
#     coverage: int = 100,
#     hybrid: bool = False,
#     skip_logs: bool = False,
#     skip_fastq_check: bool = False,
#     skip_qc: bool = False,
#     skip_error_correction: bool = False,
#     no_miniasm: bool = False,
#     skip_pseudogene_correction: bool = False,
#     skip_adj_correction: bool = False,
#     skip_prodigal_tf: bool = False,
#     rawproduct: bool = False,
#     centre: str = "Bactopia",
#     min_contig_length: int = 500,
#     amr_plus: bool = False,
# ) -> LatchDir:
#     # example opening a LatchFile
#     with open(Path(fastq_one), "w") as f:
#         lines = f.readlines() 
    
#         # ... Logic Here ...
    
#         local_output_dir = Path("/root/outputs")
  
#         # example returning Flyte directory
#         return LatchDir(
#             str(local_output_dir.resolve()),
#             remote_directory=_fmt_dir(output_dir.remote_source), # remote_source is a LatchDir? What does this do?
#         )


# @workflow
# def bactopia_wf(
#     output_dir: LatchDir,
#     sample_name: List[Union[str, int]] = "sample1",
#     fastq_one: Optional[LatchFile] = None,
#     fastq_two: Optional[LatchFile] = None,
#     input_dir: Optional[LatchDir] = None,
#     genome_size: Optional[int] = None,
#     species: Species = Species.none,
#     species_genome_size: SpeciesGenomeSize = SpeciesGenomeSize.mash,
#     ask_merlin: bool = False,
#     coverage: int = 100,
#     hybrid: bool = False,
#     skip_logs: bool = False,
#     skip_fastq_check: bool = False,
#     skip_qc: bool = False,
#     skip_error_correction: bool = False,
#     no_miniasm: bool = False,
#     skip_pseudogene_correction: bool = False,
#     skip_adj_correction: bool = False,
#     skip_prodigal_tf: bool = False,
#     rawproduct: bool = False,
#     amr_plus: bool = False,
#     centre: str = "Bactopia",
#     min_contig_length: int = 500,
# ) -> LatchDir:

#     return bactopia_tsk(
#         fastq_one=fastq_one,
#         fastq_two=fastq_two,
#         input_dir=input_dir,
#         output_dir=output_dir,
#         sample_name=sample_name,
#         genome_size=genome_size,
#         species=species,
#         species_genome_size=species_genome_size,
#         ask_merlin=ask_merlin,
#         coverage=coverage,
#         hybrid=hybrid,
#         skip_logs=skip_logs,
#         skip_fastq_check=skip_fastq_check,
#         skip_qc=skip_qc,
#         skip_error_correction=skip_error_correction,
#         no_miniasm=no_miniasm,
#         skip_pseudogene_correction=skip_pseudogene_correction,
#         skip_adj_correction=skip_adj_correction,
#         skip_prodigal_tf=skip_prodigal_tf,
#         rawproduct=rawproduct,
#         centre=centre,
#         min_contig_length=min_contig_length,
#         amr_plus=amr_plus,
#     )

#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????





# #Input bam file to be assembled.
# @small_task
# def assemble_bam(bam_file: LatchFile) -> LatchFile:
#     """
#     Assemble a BAM file.
#     """
#     bam_file = Path(bam_file)
#     bam_file_name = bam_file.stem
#     bam_file_dir = bam_file.parent
#     bam_file_assembled = bam_file_dir / f"{bam_file_name}_assembled.bam"
#     bam_file_assembled.touch()
#     subprocess.run(
#         [   "bamtools",
#             "assemble",
#             "-i",
#             bam_file,
#             "-o",
#             bam_file_assembled,
#         ],
#         check=True,
#     )
#     return bam_file_assembled
    
# #Sort the assembled BAM file.
# @small_task
# def sort_bam(bam_file: LatchFile) -> LatchFile:
#     """
#     Sort a BAM file.
#     """
#     bam_file = Path(bam_file)
#     bam_file_name = bam_file.stem
#     bam_file_dir = bam_file.parent
#     bam_file_sorted = bam_file_dir / f"{bam_file_name}_sorted.bam"
#     bam_file_sorted.touch()
#     subprocess.run(
#         [   "bamtools",
#             "sort",
#             "-i",
#             bam_file,
#             "-o",
#             bam_file_sorted,
#         ],
#         check=True,
#     )
#     return bam_file_sorted

#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

# @small_task
# def assembly_task(read1: LatchFile, read2: LatchFile) -> LatchFile:

#     # A reference to our output.
#     sam_file = Path("covid_assembly.sam").resolve()

#     _bowtie2_cmd = [
#         "bowtie2/bowtie2",
#         "--local",
#         "-x",
#         "wuhan",
#         "-1",
#         read1.local_path,
#         "-2",
#         read2.local_path,
#         "--very-sensitive-local",
#         "-S",
#         str(sam_file),
#     ]

#     subprocess.run(_bowtie2_cmd)

#     return LatchFile(str(sam_file), "latch:///covid_assembly.sam")


# @small_task
# def sort_bam_task(sam: LatchFile) -> LatchFile:

#     bam_file = Path("covid_sorted.bam").resolve()

#     _samtools_sort_cmd = [
#         "samtools",
#         "sort",
#         "-o",
#         str(bam_file),
#         "-O",
#         "bam",
#         sam.local_path,
#     ]

#     subprocess.run(_samtools_sort_cmd)

#     return LatchFile(str(bam_file), "latch:///covid_sorted.bam")





# @workflow
# def assemble_and_sort(read1: LatchFile, read2: LatchFile) -> LatchFile:
#     """Description...

#     markdown header
#     ----

#     Write some documentation about your workflow in
#     markdown here:

#     > Regular markdown constructs work as expected.

#     # Heading

#     * content1
#     * content2

#     __metadata__:
#         display_name: Assemble and Sort FastQ Files
#         author:
#             name:
#             email:
#             github:
#         repository:
#         license:
#             id: MIT

#     Args:

#         read1:
#           Paired-end read 1 file to be assembled.

#           __metadata__:
#             display_name: Read1

#         read2:
#           Paired-end read 2 file to be assembled.

#           __metadata__:
#             display_name: Read2
#     """
    # sam = assembly_task(read1=read1, read2=read2)
    # return sort_bam_task(sam=sam)