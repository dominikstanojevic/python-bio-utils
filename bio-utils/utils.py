from Bio import SeqIO

from pathlib import Path
import re

from typing import *


def get_files(path: Union[Path, str],  
              extension: str,
              recursive: bool=False) -> Generator[Path, None, None]:
    '''
    Returns generator over files in the root path that match the given
    extension.

    If path is a file, this function checks if the path ends with the given
    extension. If recursive flag is set to True, subfolders are searched for 
    files that ends with the give extension.

    Args:
        path: Root path or path to the file.
        extension: Extension to be searched for.
        recursive: Flag to indicate that the search should be performed 
                   recursively.
    Returns:
        All files which path ends with the give extension.
    '''
    if isinstance(path, str):
        path = Path(path)

    if path.is_file():
        if path.name.endswith(extension):
            yield path
    else:
        if recursive:
            files = path.glob(f'**/*{extension}')
        else:
            files = path.glob(f'*{extension}')

        yield from files


MotifPositions = dict[str, Tuple[Set[int], Set[int]]]
def build_reference_idx(path: Union[Path, str], 
                        motif: str,
                        rel_idx: int) -> MotifPositions:
    '''
    Generates position for specific motif.

    This function generates position dictionary for specific motif and relative 
    index (specific base in the motif). Positions are stored separately for forward 
    and reverse strand. Keys to the dictionary are seuqnce names that are 
    present in the FASTA file.

    Args:
        path: Path to the FASTA file.
        motif: Motif to be searched for in the sequences.
        rel_idx: Index relative to the start of the motif.
    Returns:
        Dictionary of [seq_name, positions] elements in which positions are 
        defined as tuple of position on forward and on reverse strand
    '''
    positions = OrderedDict()

    for record in SeqIO.parse(path, 'fasta'):
        contig = record.name
        length = len(record.seq)

        fwd_pos = {m.start() + rel_idx for m in re.finditer(motif, str(record.seq), re.I)}

        def pos_for_rev(i):
            return length - (i + rel_idx) - 1
        rev_pos = {pos_for_rev(m.start()) 
                   for m in re.finditer(motif, str(record.seq.reverse_complement()), re.I)}

        positions[contig] = (fwd_pos, rev_pos)

    return positions