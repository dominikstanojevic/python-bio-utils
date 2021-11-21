from ont_fast5_api.fast5_file import Fast5File
import numpy as np

from dataclasses import dataclass

from typing import *

from .alignment import Strand


@dataclass
class AlignmentStats:
    n_matches: int
    n_mismatches: int
    n_insertions: int
    n_deletions: int


@dataclass
class AlignmentInfo:
    ctg: str
    ref_start: int
    ref_end: int
    n_clipped_start: int
    n_clipped_end: int
    strand: Strand
    stats: AlignmentStats


def get_event_table(file: Fast5File) -> Optional[np.ndarray]:
    '''
    Retrieves event table from the latest analysis.

    This function retrieves event table from the:
    "/Analyses/RawGenomeCorrected_###/BaseCalled_template/Events" where ### is
    the index of the latest analysis.

    Args:
        file: Single FAST5 file
    Return:
        Event table, None if RawGenomeCorrected is not present
    '''
    corrected_analysis = file.get_latest_analysis('RawGenomeCorrected')
    if corrected_analysis is None:
        return

    return file.get_analysis_dataset(corrected_analysis, 'Events')


def get_alignment_info(file: Fast5File) -> Optional[np.ndarray]:
    '''
    Retrieves alignment info attributes used for tombo alignment

    Retrieves alignment attributes from:
    "/Analyses/RawGenomeCorrected_###/BaseCalled_template/Alignment" where ###
    is the index of the latest analysis.

    Args:
        file: Single FAST5 file
    Return:
        AlignmentInfo object representing alignment information, None if 
        RawGenomeCorrected is not present
    '''
    corrected_analysis = file.get_latest_analysis('RawGenomeCorrected')
    if corrected_analysis is None:
        return

    attrs = file.get_analysis_attributes(f'{corrected_analysis}/BaseCalled_template/Alignment')
    ctg = attrs['mapped_chrom']
    rstart = attrs['mapped_start']
    rend = attrs['mapped_end']
    clipped_start = attrs['clipped_bases_start']
    clipped_end = attrs['clipped_bases_end']
    strand = attrs['mapped_strand']

    n_matches = attrs['num_matches']
    n_mismatches = attrs['num_mismatches']
    n_ins = attrs['num_insertions']
    n_del = attrs['num_deletions']
    stats = AlignmentStats(n_matches, n_mismatches, n_ins, n_del)

    return AlignmentInfo(ctg, rstart, rend, clipped_start, clipped_end, strand, 
                         stats)
