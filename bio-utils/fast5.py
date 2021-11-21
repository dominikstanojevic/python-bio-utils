from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.fast5_interface import get_fast5_file
from ont_fast5_api.fast5_read import Fast5Read
import numpy as np

from pathlib import Path

from typing import *


def get_reads(path: Union[str, Path], 
              mode='r') -> Generator[Union[Fast5Read, Fast5File], None, None]:
    '''Retrieves all reads stored in FAST5 file.

    For single FAST5 file it retrieves itself. For muti FAST5 file it retrieves
    all the reads stored in the file.

    Args:
        path: Path to the FAST5 file.
        mode: Mode for FAST5 file opening.
    
    Returns:
        Generator of all reads present in the FAST5 file, or the FAST5 file 
        itself.
    '''
    if isinstance(path, Path):
        path = str(path)

    with get_fast5_file(path, mode) as f5:
        yield from f5.get_reads()


def get_raw_signal(read: Fast5Read, continuous: bool=True) -> np.ndarray:
    '''
    Retrives signal from the fast5 file.

    Args:
        read: Fast5 read
        continuous: Flag indicating if returned signal should be discrete or
                    continuous
    Returns:
        Signal for the given read
    '''
    return read.get_raw_data(scale=continuous)


def get_fastq(read: Fast5Read) -> Optional[str]:
    '''
    Retrieves fastq string from the basecall analysis group.

    The string is retrieved from the latest analysis.
    Path to the dataset: "/Analyses/Basecall_1D_###/BaseCalled_template/Fastq"
    where ### is the index of the latest analysis.

    Args:
        read: FAST5 read
    Returns:
        FASTQ string from the latest analysis, None if not present
    '''
    basecall = read.get_latest_analysis('Basecall_1D')
    if basecall is None:
        return

    fastq = read.get_analysis_dataset(basecall, 'BaseCalled_template/Fastq')
    return fastq


def get_block_stride(read: Fast5Read) -> Optional[int]:
    '''
    Retrieves block size for the latest basecall.

    Path to the "block_stride" attribute:
    "/Analyses/Basecall_1D_###/Summary/basecall_1d_template" where ### is the 
    index of the latest analysis.

    Args:
        read: FAST5 read
    Returns:
        Block stride, None if analysis is not present
    '''
    bc_analysis = read.get_latest_analysis('Basecall_1D')
    if bc_analysis is None:
        return

    bc_summary = read.get_summary_data(bc_analysis)
    return bc_summary['basecall_1d_template']['block_stride']

def get_raw_start_index(read: Fast5Read) -> Optional[int]:
    '''
    Retrieves starting index for the first basecalled base.

    The starting index is written as an "first_sample_template" attribute in
    "/Analyses/Segmentation_###/Summary/segmentation" where ### is the index of
    the latest analysis.

    Args:
        read: FAST5 read
    Returns:
        The starting raw index for the latest analysis, None if analysis is not
        present.    
    '''
    segmentation = read.get_latest_analysis('Segmentation')
    if segmentation is None:
        return

    summary = read.get_summary_data(segmentation)
    return summary['segmentation']['first_sample_template']


def get_move_table(read: Fast5Read) -> Optional[np.ndarray]:
    '''
    Returns move table.

    Retrieves move table from "/Analyses/Basecall_1D_###/BaseCalled_template/Move"
    where ### is the latest analysis index. If no analysis is present or move
    table was not written during basecalling, it returns None.

    Args:
        read: FAST5 read
    Returns:
        Move table, None if not present
    '''
    bc_analysis = read.get_latest_analysis('Basecall_1D')
    if bc_analysis is None:
        return

    return read.get_analysis_dataset(bc_analysis, 'BaseCalled_template/Move')


def get_offset_scale(read: Fast5Read) -> Tuple[float, float]:
    '''
    Retrieves offset and scale for the given read.

    This function offset and scale used for converting discrete levels to 
    signal values in pA. Equation: pA = (levels + offset) * scale

    Args:
        read: FAST5 read
    Returns:
        Tuple of offset and scale for the given read
    '''
    channel_info = read.get_channel_info()

    digitisation = channel_info['digitisation']
    rng = channel_info['range']
    scale = rng / digitisation

    return channel_info['offset'], scale


def sequence_to_signal(move_table: np.ndarray,
                       raw_start_idx: int,
                       block_stride: int) -> np.ndarray:
    '''Maps basecalled sequence to the signal.

    This function maps indices of basecalled sequence to the points corresponding
    to the each base. Mapped element is the starting signal point for every base.
    Signal points are processed in blocks - bases are not called for individual
    points, only for one or more blocks.

    Args:
        move_table: Move table
        raw_start_idx: Index for the first signal point
        block_stride: Block size
    Returns:
        Sequence-to-signal mapping
    '''
    move_table = np.append(move_table, 1)
    return move_table.nonzero()[0] * block_stride + raw_start_idx


def parse_fastq(fastq: str) -> Tuple[str, np.ndarray]:
    '''Parses FASTQ string and returns both sequence and qualities.

    This function parses FASTQ string (possibly obtained from FAST5 file) and
    returns genomic sequence and corresponding base qualities (Phred).

    Args:
        fastq: FASTQ string
    Returns:
        Genomic sequence and corresponding qualities.
    '''
    data = fastq.strip().split('\n')

    sequence = data[1]
    qualities = np.ndarray([ord(c) - 33 for c in data[3]], dtype=np.uint8)

    return sequence, qualities


def normalize_mad(signal: np.ndarray,
                  scale_factor: Optional[float]=1.4826) -> np.ndarray:
    '''Performs signal normalization using median absolute deviation.

    This function takes signal and performs signal normalization by shifting
    values using median and scaling using median absolute deviation. 
    Optional constant scale factor can be used to estimate standard deviation 
    (s = factor * MAD). Default factor is factor used for normally distributed
    data (1.4826).

    Args:
        signal: Signal to be normalized.
        scale_factor: Scale factor (default: 1.4826).
    '''
    if scale_factor is None:
        scale_factor = 1

    med = np.median(signal)
    shifted = signal - med
    mad = np.median(np.abs(shifted))

    return shifted / (scale_factor * mad)

