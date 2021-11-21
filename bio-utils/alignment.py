from mappy import Aligner, Alignment
import numpy as np

from enum import Enum
from pathlib import Path

from typing import *


class Strand(str, Enum):
    FORWARD = '+',
    REVERSE = '-'


def get_aligner(reference: Union[Path, str],
                preset: str='map-ont', 
                **kwargs) -> Optional[Aligner]:
    '''
    Returns aligner with for the given reference and preset.

    Full list of arguments can be found here: 
    `Link text <https://pypi.org/project/mappy/>`. This function also checks
    if the correct reference file was given, returning NoneType otherwise. Mappy
    will not give raise any error if the wrong reference file was provided.

    Args:
        reference: Reference file in FASTA fromat. Can be gzipped.
        preset: Preset for the aligner. Default is map-ont.
        kwargs: Keyword arguments for the mappy.Aligner class.
    Returns:
        Mappy aligner for the given reference, NoneType otherwise
    '''
    aligner = Aligner(reference, preset=preset, **kwargs)
    
    seq_names = aligner.seq_names
    if seq_names is None or len(seq_names) == 0:
        return

    return aligner

def align(query: str,
          aligner: Aligner,
          best: bool=True) -> Union[List[Alignment], Alignment, None]:
    '''
    Aligns the query sequence to the reference.

    This function aligns the given sequence to the reference indexed by aligner.
    Alignment in mappy returns k best alignments (parameter for aligner). Even
    if k=1, map function will return all primary alignments (possibly more than)
    one. If best flag is set, this function will return only the best primary 
    alignment.

    Args:
        query: The query sequence that is going to be aligned
        aligner: aligner that indexed a reference sequence
        best: Flag indicating that only the best primary alignment should be
              returned
    Returns:
        All alignments if best flag is not set, best primary alignment if the 
        flag is set, None if there is no alignments
    '''
    alignments = aligner.map(query)

    if len(alignments) == 0:
        return
    if best:
        return alignments[0]

    return alignments


def reference_to_query(alignment: Alignment) -> List[int]:
    '''
    Maps reference indices to the query indices using cigar string.

    This function maps every reference base in the alignment to the base(s) in
    the query sequence. Elements of the list are starting query indices for every
    reference base. If insertion occured, reference will span multiple query bases
    (map[pos+1] - map[pos] > 1). If deletion occured there will be no mapped query
    bases (map[pos+1] - map[pos] = 0). For match or mismatch map[pos+1] - map[pos] = 1.
    This mapping was generated from cigar string. Last element in this mapping is
    the exclusive end.

    Args:
        alignment: mappy.Alignment object representing alignment
    Returns:
        Reference to query base mapping
    '''
    ref_len = alignment.r_en - alignment.r_st
    cigar = alignment.cigar if alignment.strand == 1 else reversed(alignment.cigar)
    rpos, qpos = 0, alignment.q_st

    ref_to_query = np.empty((ref_len + 1,), dtype=np.uint32)
    for l, op in cigar:
        if op == 0 or op == 7 or op == 8:  # Match or mismatch
            for i in range(l):
                ref_to_query[rpos + i] = qpos + i
            rpos += l
            qpos += l
        elif op == 1:  # Insertion
            qpos += l
        elif op == 2:  # Deletion
            for i in range(l):
                ref_to_query[rpos + i] = qpos
            rpos += l
        else:
            raise TypeError('Invalid cigar operation.')

    ref_to_query[rpos] = qpos  # Add the last one (excluded end)

    return ref_to_query