#!/usr/bin/env python3

import os
import sys
import typing
import json
import gzip
import struct
import multiprocessing as mp

from glob import iglob
from collections import OrderedDict

import numpy as np
from scipy.spatial.distance import pdist, squareform

"""
    python requirements:
        numpy
        scipy
        flexx

    system requirements:
        bgzip

    optional requirements:
        pypy

    system requirements pypy:
        libatlas-base-dev libblas3 liblapack3 liblapack-dev libblas-dev gfortran
"""


"""
    1 thread  - no alignment
    real    220m11.411s
    user    436m20.040s
    sys       0m41.524s

    4 threads - no alignment
    real     80m24.037s
    user    373m12.301s
    sys       0m40.070s

    6 threads - no alignment
    real     81m47.819s
    user    498m14.096s
    sys       1m27.114s

    8 threads - no alignment
    real     57m19.298s
    user    398m 2.372s
    sys       3m 0.994s

    ===========
    6 threads - no alignment - 20 bins
    real     6m12.305s
    user    39m 2.290s
    sys      0m 4.687s

    6 threads - w/ alignment - 20 bins
    real     7m18.622s
    user    45m43.042s
    sys      0m 5.512s

    ==========
    6 threads - w/ alignment
    real     99m40.466s
    user    611m 1.784s
    sys       1m46.982s

"""


"""
    w/o alignment
    1023M 360_merged_2.50.vcf.gz
    5.3M 360_merged_2.50.vcf.gz.250000.000000.SL2.50ch00.npz
    23M 360_merged_2.50.vcf.gz.250000.000001.SL2.50ch01.npz
    12M 360_merged_2.50.vcf.gz.250000.000002.SL2.50ch02.npz
    15M 360_merged_2.50.vcf.gz.250000.000003.SL2.50ch03.npz
    18M 360_merged_2.50.vcf.gz.250000.000004.SL2.50ch04.npz
    16M 360_merged_2.50.vcf.gz.250000.000005.SL2.50ch05.npz
    11M 360_merged_2.50.vcf.gz.250000.000006.SL2.50ch06.npz
    14M 360_merged_2.50.vcf.gz.250000.000007.SL2.50ch07.npz
    15M 360_merged_2.50.vcf.gz.250000.000008.SL2.50ch08.npz
    15M 360_merged_2.50.vcf.gz.250000.000009.SL2.50ch09.npz
    16M 360_merged_2.50.vcf.gz.250000.000010.SL2.50ch10.npz
    15M 360_merged_2.50.vcf.gz.250000.000011.SL2.50ch11.npz
    19M 360_merged_2.50.vcf.gz.250000.000012.SL2.50ch12.npz
    4.0K 360_merged_2.50.vcf.gz.250000.npz
    496K 360_merged_2.50.vcf.gz.csi
    2.8M 360_merged_2.50.vcf.gz.gzi
    8.0K 360_merged_2.50.vcf.gz.gzj

    w/ alignment
    5.7M 360_merged_2.50.vcf.gz.250000.000000.SL2.50ch00.npz
    25M 360_merged_2.50.vcf.gz.250000.000001.SL2.50ch01.npz
    13M 360_merged_2.50.vcf.gz.250000.000002.SL2.50ch02.npz
    16M 360_merged_2.50.vcf.gz.250000.000003.SL2.50ch03.npz
    19M 360_merged_2.50.vcf.gz.250000.000004.SL2.50ch04.npz
    18M 360_merged_2.50.vcf.gz.250000.000005.SL2.50ch05.npz
    12M 360_merged_2.50.vcf.gz.250000.000006.SL2.50ch06.npz
    16M 360_merged_2.50.vcf.gz.250000.000007.SL2.50ch07.npz
    16M 360_merged_2.50.vcf.gz.250000.000008.SL2.50ch08.npz
    17M 360_merged_2.50.vcf.gz.250000.000009.SL2.50ch09.npz
    17M 360_merged_2.50.vcf.gz.250000.000010.SL2.50ch10.npz
    16M 360_merged_2.50.vcf.gz.250000.000011.SL2.50ch11.npz
    20M 360_merged_2.50.vcf.gz.250000.000012.SL2.50ch12.npz
    4.0K 360_merged_2.50.vcf.gz.250000.npz
"""

DEBUG                         = False
DEBUG_MAX_BIN                 = 10
DEFAULT_BIN_SIZE              = 250_000
DEFAULT_SAVE_ALIGNMENT        = True
DEFAULT_METRIC                = 'jaccard'
DEFAULT_DISTANCE_TYPE_MATRIX  = np.float32
DEFAULT_COUNTER_TYPE_MATRIX   = np.uint16
DEFAULT_COUNTER_TYPE_PAIRWISE = np.uint32
DEFAULT_POSITIONS_TYPE        = np.uint32
DEFAULT_THREADS               = mp.cpu_count()

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
METRIC_RAW_NAME = 'RAW'
METRIC_VALIDS                 = [
    'braycurtis' , 'canberra'  , 'chebyshev'    , 'cityblock'     ,
    'correlation', 'cosine'    , 'dice'         , 'euclidean'     ,
    'hamming'    , 'jaccard'   , 'jensenshannon', 'kulsinski'     ,
    'mahalanobis', 'matching'  , 'minkowski'    , 'rogerstanimoto',
    'russellrao' , 'seuclidean', 'sokalmichener', 'sokalsneath'   ,
    'sqeuclidean', 'yule'      ,
    METRIC_RAW_NAME
]


MatrixType           = typing.OrderedDict[typing.Tuple[str,str],int]
ChromosomeMatrixType = typing.OrderedDict[int, typing.List[int]]
BinSnpsType          = typing.OrderedDict[int, int]
BinPairwiseCountType = typing.OrderedDict[int, typing.List[int]]
BinPositionType      = typing.OrderedDict[int, typing.List[int]]
BinPositionTypeInt   = typing.OrderedDict[int, typing.List[typing.List[int]]]
BinAlignmentType     = typing.OrderedDict[int, typing.List[str]]
BinAlignmentTypeInt  = typing.OrderedDict[int, typing.List[typing.List[str]]]
TriangleIndexType    = typing.OrderedDict[typing.Tuple[int,int], int]
IUPACType            = typing.Dict[typing.FrozenSet[str], str]
SampleNamesType      = typing.List[str]
ChromosomeNamesType  = typing.List[str]
ChromosomesType      = typing.OrderedDict[str, "Chromosome"]


class Chromosome():
    """
        import reader
        chrom = reader.Chromosome(vcf_name="../data/360_merged_2.50.vcf.gz", bin_width=250_000, metric="RAW", chromosome_order=0, chromosome_name="SL2.50ch00")
        chrom.filename
        chrom.exists
        chrom.load()
        chrom.matrixNp.shape
        # (87, 64980)
        chrom.matrixNp
        # array([
        # [ 58,  18,  18, ...,   0,   0,   0],
        # [146,  44,  44, ...,   0,   0,   0],
        # [ 92,  18,  16, ...,   0,   0,   0],
        # ...,
        # [ 12,  28,  24, ...,   0,   0,   0],
        # [ 10,  14,  10, ...,   0,   0,   0],
        # [ 68,  78,  64, ...,   0,   0,   0]], dtype=uint16)
        chrom.matrixNp[0,:]
        # array([58, 18, 24, ..., 12,  0,  0], dtype=uint16)
        chrom.matrixNp.shape
        # (87, 64980)
        chrom.sample_count
        # 361
        chrom.matrixNp.sum(axis=1)
        # array([14135926, 23254563, 19677883, 15874453,  8422700, 17916464,
        #        13698021, 16151376,  7881003, 16836865, 19910092, 17204299,
        #        16263240,  7178839, 14871309,  9396853,  5977816, 10596943,
        #        15257031, 20200738, 16195559, 11704718,  8208319, 18487166,
        #        33227260, 19990183,  9529424,  8905346,  6964558,  4755606,
        #         7747795,   511011,  4415965,  5857671,  6814355,  6168579,
        #         1453627,  1744796,  4442980,  2088842,  5274180,  2468293,
        #         1546854,  2285387,  1936388,  1377816,  2966110,  2667990,
        #         4795700,  3008045,  4955744,  4958694,  6057739,  5587246,
        #         8511794,  7176901,  4027879,  4136236,  3431414,  4764401,
        #         4960412,  4165828,  5062155,  3305293,  5670359,  5311667,
        #         6011455,  6375255,  5160923,  4773837,  5889080,  5730958,
        #         6715606,  6523644,  5566348,  6531633,  6160864,  5562834,
        #         5614820,  3962912,  4594549,  3875737,  4870790,  4368941,
        #         4920943,  4689962,  5360525], dtype=uint64)
        m = reader.triangleToMatrix(chrom.sample_count, chrom.matrixNp[0,:])
        m
        #
        #array([
        # [ 0, 58, 18, ..., 16, 32,  0],
        # [58,  0, 18, ..., 16, 34,  0],
        # [18, 18,  0, ..., 14, 16,  0],
        # ...,
        # [16, 16, 14, ...,  0, 12,  0],
        # [32, 34, 16, ..., 12,  0,  0],
        # [ 0,  0,  0, ...,  0,  0,  0]], dtype=uint16)
        #
        m.shape
        #
        #(361, 361)
        #
        chrom.matrix_bin_matrix_dist(1)
    """

    def __init__(self,
            vcf_name                      : str,
            bin_width                     : int,
            metric                        : str,

            chromosome_order              : int,
            chromosome_name               : str,

            type_matrix_counter           : np.dtype        = DEFAULT_COUNTER_TYPE_MATRIX,
            type_matrix_distance          : np.dtype        = DEFAULT_DISTANCE_TYPE_MATRIX,
            type_pairwise_counter         : np.dtype        = DEFAULT_COUNTER_TYPE_PAIRWISE,
            type_positions                : np.dtype        = DEFAULT_POSITIONS_TYPE
        ):

        self.vcf_name                     : str             = vcf_name
        self.bin_width                    : int             = bin_width
        self.chromosome_order             : int             = chromosome_order
        self.chromosome_name              : str             = chromosome_name
        self.metric                       : str             = metric

        self.type_matrix_counter          : np.dtype        = type_matrix_counter
        self.type_matrix_distance         : np.dtype        = type_matrix_distance
        self.type_pairwise_counter        : np.dtype        = type_pairwise_counter
        self.type_positions               : np.dtype        = type_positions

        self.type_matrix_counter_max_val  : int             = np.iinfo(self.type_matrix_counter  ).max
        self.type_matrix_distance_max_val : int             = np.finfo(self.type_matrix_distance ).maxexp
        self.type_pairwise_counter_max_val: int             = np.iinfo(self.type_pairwise_counter).max
        self.type_positions_max_val       : int             = np.iinfo(self.type_positions       ).max

        self.matrix_size                  : int             = None
        self.bin_max                      : int             = None
        self.bin_min                      : int             = None
        self.bin_count                    : int             = None

        self.bin_snps_min                 : int             = None
        self.bin_snps_max                 : int             = None

        self.chromosome_snps              : int             = None
        self.chromosome_first_position    : int             = None
        self.chromosome_last_position     : int             = None

        self.sample_names                 : SampleNamesType = None
        self.sample_count                 : int             = None

        self.matrixNp                     : np.ndarray      = None
        self.binsnpNp                     : np.ndarray      = None
        self.pairwiNp                     : np.ndarray      = None
        self.alignmentNp                  : np.ndarray      = None
        self.positionNp                   : np.ndarray      = None

        self._is_loaded                   : bool            = False

        assert metric in METRIC_VALIDS

    @property
    def filename(self) -> str:
        basefolder = self.vcf_name + "_ib"
        dirname    = os.path.join(basefolder, str(self.bin_width), self.metric)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        return f"{dirname}{os.path.sep}ib_{self.chromosome_order:06d}.{self.chromosome_name}.npz"

    @property
    def exists(self) -> bool:
        return os.path.exists(self.filename)

    @property
    def is_loaded(self) -> bool:
        return self._is_loaded

    @property
    def matrix(self) -> np.ndarray:
        return self.matrixNp

    @property
    def matrix_dtype(self) -> np.dtype:
        return self.matrixNp.dtype

    def _get_infos(self):
        names  = [
            "bin_count"                  , "bin_min"                      , "bin_max"                      , "bin_width"               ,
            "bin_snps_min"               , "bin_snps_max"                 ,
            "chromosome_snps"            , "chromosome_order"             , "chromosome_first_position"    , "chromosome_last_position",
            "matrix_size"                , "sample_count"                 , 
            "type_matrix_counter_max_val", "type_matrix_distance_max_val" , "type_pairwise_counter_max_val", "type_positions_max_val"
        ]
        vals = [getattr(self,n) for n in names]

        return (names, vals)

    def _get_meta(self):
        type_matrix_counter_name   = self.type_matrix_counter.__name__
        type_matrix_distance_name  = self.type_matrix_distance.__name__
        type_pairwise_counter_name = self.type_pairwise_counter.__name__
        type_positions_name        = self.type_positions.__name__

        names  = ["vcf_name"      , "metric"   , "chromosome_name"   ]
        vals   = [getattr(self,n) for n in names]

        names += ["type_matrix_counter_name", "type_matrix_distance_name", "type_pairwise_counter_name", "type_positions_name"]
        vals  += [type_matrix_counter_name  , type_matrix_distance_name  , type_pairwise_counter_name  , type_positions_name  ]

        return names, vals

    def todict(self):
        info_names,info_values     = self._get_infos()
        meta_names,meta_values     = self._get_meta()

        data =      {k:v for k,v in zip(info_names, info_values)}
        data.update({k:v for k,v in zip(meta_names, meta_values)})
        
        data['sample_names'] = self.sample_names
        data['countMatrix '] = self.matrixNp.tolist()
        data['countTotals']  = self.binsnpNp.tolist()
        data['countPairw']   = self.pairwiNp.tolist()
        data['alignments']   = self.alignmentNp.tolist()
        data['positions']    = self.positionNp.tolist()

        return data

    def _sample_pos(self, sample_name: str):
        return self.sample_names.index(sample_name)

    def __repr__(self):
        return str(self)

    def __str__(self):
        res = []
        for k in [
            "vcf_name",
            "chromosome_name",
            "chromosome_snps",
            "chromosome_order",
            "chromosome_first_position",
            "chromosome_last_position",
            "bin_count",
            "bin_width",
            "bin_min",
            "bin_max",
            "bin_snps_min",
            "bin_snps_max",
            "matrix_size",
            "metric",
            "sample_names",
            "sample_count",
            "type_matrix_counter",
            "type_matrix_distance",
            "type_pairwise_counter",
            "type_positions",
            "type_matrix_counter_max_val",
            "type_matrix_distance_max_val",
            "type_pairwise_counter_max_val",
            "type_positions_max_val",
            "matrixNp",
            "binsnpNp",
            "pairwiNp",
            "alignmentNp",
            "positionNp"
        ]:

            v = getattr(self, k)
            s = None
            if   isinstance(v, int):
                s = f"{v:,d}"
            elif isinstance(v, str):
                s = f"{v:s}"
            elif isinstance(v, list):
                s = f"{len(v):,d}"
            elif isinstance(v, np.ndarray):
                s = f"{str(v.shape):s}"
            else:
                s = f"{str(v):s}"
            res.append(f"  {k:.<30s}{s:.>30s}")
        return "\n".join(res)

    @staticmethod
    def _processVcf_read_header(vcf_name: str) -> typing.Tuple[SampleNamesType, int, int, TriangleIndexType]:
        with openFile(vcf_name, 'rt') as fhd:
            for line in fhd:
                line = line.strip()
                
                if len(line) <= 2:
                    continue

                if line[:2] == "##":
                    # print("header", line)
                    continue
                
                if line[0] == "#":
                    # print("title", line)
                    cols         = line[1:].split("\t")
                    sample_names = cols[9:]
                    sample_count = len(sample_names)
                    # matrix_size  = calculateMatrixSize(sample_count)
                    _, matrix_size, indexes = triangleToIndex(sample_count)

                    print( "sample_names", ",".join(sample_names))
                    print(f"num samples  {sample_count:12,d}")
                    print(f"matrix_size  {matrix_size :12,d}")
                    print(f"indexes      {len(indexes):12,d}")
                    return sample_names, sample_count, matrix_size, indexes

                else:
                    raise ValueError("data before header error", line)

    @staticmethod
    def _processVcf_read_chrom(
            vcf_name             : str,
            bin_width            : int,
            metric               : str,

            chromosome_order     : int,
            chromosome_name      : str,

            matrix_size          : int,
            indexes              : TriangleIndexType,
            
            sample_names         : SampleNamesType,
            sample_count         : int,

            save_alignment       : bool,
            IUPAC                : IUPACType,
            diff_matrix          : MatrixType,

            type_matrix_counter  : np.dtype,
            type_pairwise_counter: np.dtype,
            type_positions       : np.dtype
        ) -> typing.Tuple[int, int]:

        print(f"reading {chromosome_name}")

        bgzip = BGzip(vcf_name)

        bin_snps                 : BinSnpsType          = OrderedDict()
        bin_pairwise_count       : BinPairwiseCountType = OrderedDict()
        bin_alignment            : BinAlignmentTypeInt  = OrderedDict()
        bin_positions            : BinPositionTypeInt   = OrderedDict()
        chromosome_matrix        : MatrixType           = OrderedDict()
        chromosome_first_position: int = 0
        chromosome_last_position : int = 0
        chromosome_snps          : int = 0

        lastBinNum = None
        for line in bgzip.get_chromosome(chromosome_name):
            # print(line)
            cols       = line.split("\t")
            chrom      = cols[0]
            ref        = cols[3]
            alts       = cols[4].split(",")
            opts       = [ref] + alts
            pos        = int(cols[1])
            samples    = cols[9:]
            binNum     = pos // bin_width
            assert len(samples) == sample_count
            
            if len(ref) != 1:
                print('H')
                continue

            if any([len(a) > 1 for a in alts]):
                print('h')
                continue

            if binNum not in chromosome_matrix:
                if DEBUG:
                    if len(chromosome_matrix) > DEBUG_MAX_BIN:
                        break

                print(f"  New bin: {chrom} :: {pos:12,d} => {binNum:12,d}")

                if lastBinNum is not None:
                    chromosome_matrix[ lastBinNum] = np.array(chromosome_matrix [lastBinNum], type_matrix_counter  )
                    bin_pairwise_count[lastBinNum] = np.array(bin_pairwise_count[lastBinNum], type_pairwise_counter)
                    # bin_snps          [lastBinNum] = np.array(bin_snps          [lastBinNum], type_pairwise_counter)

                    if save_alignment:
                        bin_alignment[lastBinNum]  = ["".join(b) for b in bin_alignment[lastBinNum]]
                        # bin_positions[lastBinNum]  = [p for p in bin_positions  [lastBinNum] if p != -1]
                        bin_positions[lastBinNum]  = np.array(bin_positions     [lastBinNum], type_positions       )

                lastBinNum                 = binNum
                chromosome_matrix[binNum]  = [0] * matrix_size
                bin_pairwise_count[binNum] = [0] * sample_count
                bin_snps[binNum]           =  0

                bin_alignment[binNum]      = None
                if save_alignment:
                    bin_alignment[binNum]  = [[] for _ in range(sample_count)]
                    bin_positions[binNum]  = []

            samples                   = [s.split(";")[0]            for s in samples]
            samples                   = [s if len(s) == 3 else None for s in samples]
            # samples = [s.replace("|", "").replace("/", "") if s is not None else None for s in samples]
            samples                   = [tuple([int(i) for i in s.replace("/", "|").split("|")]) if s is not None else None for s in samples]
            vals                      = chromosome_matrix[binNum]
            paiw                      = bin_pairwise_count[binNum]
            chromosome_snps          += 1
            bin_snps[binNum]         += 1
            chromosome_last_position  = pos

            aling                     = None
            if save_alignment:
                aling                 = bin_alignment[binNum]
                bin_positions[binNum].append(pos)

            if chromosome_first_position == 0:
                chromosome_first_position = pos

            for sample1num in range(sample_count):
                sample1 = samples[sample1num]

                if sample1 is None:
                    if save_alignment:
                        aling[sample1num].append('N')
                    continue
                elif save_alignment:
                    alts1 = frozenset([opts[s] for s in sample1])
                    nuc   = IUPAC[alts1]
                    aling[sample1num].append(nuc)

                for sample2num in range(sample1num+1,sample_count):
                    sample2 = samples[sample2num]

                    if sample2 is None:
                        continue

                    k = (sample1, sample2) if sample1 <= sample2 else (sample2, sample1)

                    value = diff_matrix.get(k, None)
                    if value is None:
                        print(line)
                        raise ValueError("multiallelic", k)

                    pairind           = indexes[(sample1num,sample2num)]
                    vals[pairind   ] += value
                    paiw[sample1num] += value
                    paiw[sample2num] += value

        print(f"cleaning {chromosome_name}")

        chromosome_matrix[ lastBinNum] = np.array(chromosome_matrix [lastBinNum], type_matrix_counter  )
        bin_pairwise_count[lastBinNum] = np.array(bin_pairwise_count[lastBinNum], type_pairwise_counter)
        # bin_snps          [lastBinNum] = np.array(bin_snps          [lastBinNum], type_pairwise_counter)

        if save_alignment:
            bin_alignment[lastBinNum]  = ["".join(b) for b in bin_alignment[lastBinNum]]
            # bin_positions[lastBinNum]  = [p for p in bin_positions  [lastBinNum] if p != -1]
            bin_positions[lastBinNum]  = np.array(bin_positions     [lastBinNum], type_positions       )
            position_max_size          = max([p.shape[0] for p in bin_positions.values()])

            for binNum, binval in bin_positions.items():
                binz                   = np.zeros(position_max_size, type_positions)
                binz[:binval.shape[0]] = binval
                bin_positions[binNum]  = binz

        chromosome_bin_count, chromosome_snps = Chromosome._processVcf_save_chrom_data(
            vcf_name                     = vcf_name,
            bin_width                    = bin_width,
            metric                       = metric,

            chromosome_order             = chromosome_order,
            chromosome_name              = chromosome_name,

            chromosome_snps              = chromosome_snps,
            chromosome_matrix            = chromosome_matrix,
            chromosome_first_position    = chromosome_first_position,
            chromosome_last_position     = chromosome_last_position,

            bin_alignment                = bin_alignment if save_alignment else None,
            bin_positions                = bin_positions if save_alignment else None,
            bin_snps                     = bin_snps,
            bin_pairwise_count           = bin_pairwise_count,

            sample_names                 = sample_names,
            sample_count                 = sample_count,
            matrix_size                  = matrix_size,

            type_matrix_counter          = type_matrix_counter,
            type_pairwise_counter        = type_pairwise_counter,
            type_positions               = type_positions
        )

        # chromosome_bins = chromosome.bin_max
        # chromosome_bins = len(chromosome_matrix)

        print(f"returning {chromosome_name}")

        return chromosome_bin_count, chromosome_snps

    @staticmethod
    def _processVcf_save_chrom_data(
            vcf_name                 : str,
            bin_width                : int,
            metric                   : str,

            chromosome_order         : int,
            chromosome_name          : str,

            chromosome_snps          : int,
            chromosome_matrix        : ChromosomeMatrixType,
            chromosome_first_position: int,
            chromosome_last_position : int,

            bin_alignment            : BinAlignmentType,
            bin_positions            : BinPositionType,
            bin_snps                 : BinSnpsType,
            bin_pairwise_count       : BinPairwiseCountType, 
            
            sample_names             : SampleNamesType,
            sample_count             : int,
            matrix_size              : int,

            type_matrix_counter      : np.dtype,
            type_pairwise_counter    : np.dtype,
            type_positions           : np.dtype
        ) -> typing.Tuple[int, int]:

        print(f"creating {chromosome_name}")

        # self.chromosome_names.append(chromosome_name)
        # self.chromosome_count += 1

        chromosome = Chromosome(
            vcf_name                  = vcf_name,
            bin_width                 = bin_width,
            metric                    = METRIC_RAW_NAME,
            
            chromosome_order          = chromosome_order,
            chromosome_name           = chromosome_name,

            type_matrix_counter       = type_matrix_counter,
            type_pairwise_counter     = type_pairwise_counter,
            type_positions            = type_positions
        )

        chromosome.addFromVcf(
            chromosome_snps           = chromosome_snps,
            chromosome_matrix         = chromosome_matrix,
            chromosome_first_position = chromosome_first_position,
            chromosome_last_position  = chromosome_last_position,
            
            bin_alignment             = bin_alignment,
            bin_positions             = bin_positions,
            bin_snps                  = bin_snps,
            bin_pairwise_count        = bin_pairwise_count,
            
            matrix_size               = matrix_size,
            sample_names              = sample_names
        )

        print(f"saving {chromosome_name}")

        chromosome.save()

        Chromosome._calculateDistance(
            vcf_name         = vcf_name,
            bin_width        = bin_width,
            metric           = metric,
            chromosome_order = chromosome_order,
            chromosome_name  = chromosome_name,
            chromosome       = chromosome
        )

        return chromosome.bin_count, chromosome.chromosome_snps

    @staticmethod
    def _calculateDistance(
            vcf_name              : str,
            bin_width             : int,
            metric                : str,
            chromosome_order      : int,
            chromosome_name       : str,
            chromosome            : "Chromosome" = None
        ):

        print(f"converting {chromosome_name} to {metric}")

        if chromosome is None:
            chromosome = Chromosome(
                vcf_name              = vcf_name,
                bin_width             = bin_width,
                metric                = METRIC_RAW_NAME,
                chromosome_order      = chromosome_order,
                chromosome_name       = chromosome_name,
            )
            chromosome.load()

        assert chromosome.metric == METRIC_RAW_NAME
        assert chromosome.exists

        print(f"saving {chromosome_name}")
        chromosome.save(metric=metric)

        return chromosome.bin_max, chromosome.chromosome_snps

    def addFromVcf(self,
            chromosome_snps          : int,
            chromosome_matrix        : ChromosomeMatrixType,
            chromosome_first_position: int,
            chromosome_last_position : int,

            bin_alignment            : BinAlignmentType,
            bin_positions            : BinPositionType,
            bin_snps                 : BinSnpsType,
            bin_pairwise_count       : BinPairwiseCountType,

            matrix_size              : int,
            sample_names             : SampleNamesType
        ):

        print(f"adding chromosome data: {self.chromosome_name}")

        self.matrix_size                  = matrix_size

        bin_names                         = list(chromosome_matrix.keys())

        self.bin_min                      = min(bin_names)
        self.bin_max                      = max(bin_names)
        self.bin_count                    = len(bin_names)

        self.chromosome_snps              = chromosome_snps
        self.sample_count                 = len(sample_names)
        self.sample_names                 = sample_names
        self.chromosome_first_position    = chromosome_first_position
        self.chromosome_last_position     = chromosome_last_position

        self.matrixNp                     = np.zeros((self.bin_max, self.matrix_size ), self.type_matrix_counter  )
        self.pairwiNp                     = np.zeros((self.bin_max, self.sample_count), self.type_pairwise_counter)
        self.binsnpNp                     = np.zeros( self.bin_max,                     self.type_pairwise_counter)

        self.bin_snps_min                 = min([v for v in bin_snps.values()])
        self.bin_snps_max                 = max([v for v in bin_snps.values()])

        self.alignmentNp                  = None
        self.positionNp                   = None
        if bin_alignment:
            self.alignmentNp              = np.zeros((self.bin_max, self.sample_count), np.unicode_)
            self.positionNp               = np.zeros((self.bin_max, self.bin_snps_max), self.type_positions)
        else:
            self.alignmentNp              = np.zeros(0, np.unicode_        )
            self.positionNp               = np.zeros(0, self.type_positions)

        for binNum in range(self.bin_max):
            # print("  binNum", binNum)
            chromosome_matrix_bin  = chromosome_matrix .get(binNum, np.zeros(self.matrix_size , self.type_matrix_counter  ))
            bin_pairwise_count_bin = bin_pairwise_count.get(binNum, np.zeros(self.sample_count, self.type_pairwise_counter))
            bin_snps_bin           = bin_snps          .get(binNum, 0)
            
            bin_alignment_bin      = None
            bin_positions_bin      = None
            if bin_alignment:
                # print("   creating alignment")
                bin_alignment_bin  = bin_alignment     .get(binNum, np.zeros(self.sample_count, np.unicode_        ))
                bin_positions_bin  = bin_positions     .get(binNum, np.zeros(self.bin_snps_max, self.type_positions))
            
            # print("bin_pairwise_count_bin", bin_pairwise_count_bin)
            assert not any([v > self.type_matrix_counter_max_val   for v in chromosome_matrix_bin ]), f"value went over the maximum value ({self.type_matrix_counter_max_val  }) for container {self.type_matrix_counter  }: {[v for v in chromosome_matrix_bin  if v > self.type_matrix_counter_max_val  ]}"
            assert not any([v > self.type_pairwise_counter_max_val for v in bin_pairwise_count_bin]), f"value went over the maximum value ({self.type_pairwise_counter_max_val}) for container {self.type_pairwise_counter}: {[v for v in bin_pairwise_count_bin if v > self.type_pairwise_counter_max_val]}"
            
            # binData             = [maxVal if v > maxVal else v for v in binData]
            self.matrixNp[binNum,:]    = chromosome_matrix_bin
            self.pairwiNp[binNum,:]    = bin_pairwise_count_bin
            self.binsnpNp[binNum  ]    = bin_snps_bin
            if bin_alignment:
                # print("   adding alignment")
                self.alignmentNp[binNum,:] = bin_alignment_bin
                self.positionNp [binNum,:] = bin_positions_bin
        
        self._is_loaded = True

        print(f"chromosome data added: {self.chromosome_name}")

    def save(self, metric: str = None):
        if self.metric != METRIC_RAW_NAME:
            if metric is not None:
                raise ValueError("Can't convert from one distance to another. Open RAW file for that")

        assert self.is_loaded, "chromosome not loaded"

        if metric is not None:
            assert metric in METRIC_VALIDS
            print(f"converting chromosome {self.chromosome_name} to {metric}")
            matrix_dist                      = self.convert_matrix_to_distance(metric)
            self.metric                      = metric
            self.matrixNp                    = matrix_dist

        print(f"{'saving numpy array:':.<32s}{self.filename:.>30s}")
        print(self)

        info_names, info_vals      = self._get_infos()
        meta_names, meta_vals      = self._get_meta()

        sample_namesNp             = np.array(self.sample_names, np.unicode_)
        info_namesNp               = np.array(info_names, np.unicode_)
        info_valuesNp              = np.array(info_vals , np.int64   )
        meta_namesNp               = np.array(meta_names, np.unicode_)
        meta_valuesNp              = np.array(meta_vals , np.unicode_)

        np.savez_compressed(self.filename,
            countMatrix  = self.matrixNp,
            countTotals  = self.binsnpNp,
            countPairw   = self.pairwiNp,
            alignments   = self.alignmentNp,
            positions    = self.positionNp,
            sample_names = sample_namesNp,
            info_names   = info_namesNp,
            info_values  = info_valuesNp,
            meta_names   = meta_namesNp,
            meta_values  = meta_valuesNp
        )

    def load(self):
        print(f"{'loading numpy array:':.<32s}{self.filename:.>30s}")
        data                               = np.load(self.filename, mmap_mode='r', allow_pickle=False)

        self.matrixNp                      = data['countMatrix']
        self.binsnpNp                      = data['countTotals']
        self.pairwiNp                      = data['countPairw']
        self.alignmentNp                   = data['alignments']
        self.positionNp                    = data['positions']
        
        sample_namesNp                     = data['sample_names']
        
        info_namesNp                       = data['info_names']
        info_valuesNp                      = data['info_values']

        meta_namesNp                       = data['meta_names']
        meta_valuesNp                      = data['meta_values']

        sample_names                       = sample_namesNp.tolist()

        info_names                         = info_namesNp.tolist()
        info_values                        = info_valuesNp.tolist()
        info_values                        = [int(v) for v in info_values]
        info_dict                          = {info_names[k]: info_values[k]  for k in range(len(info_names))}
        # print(info_dict)

        meta_names                         = meta_namesNp.tolist()
        meta_values                        = meta_valuesNp.tolist()
        meta_dict                          = {meta_names[k]: meta_values[k]  for k in range(len(meta_names))}
        # print(meta_dict)

        vcf_name                           = meta_dict["vcf_name"]
        assert os.path.basename(vcf_name) == os.path.basename(self.vcf_name)
        # self.vcf_name                      = vcf_name

        self.matrix_size                   = info_dict["matrix_size"]
        assert self.matrix_size == self.matrixNp.shape[1]

        self.bin_count                     = info_dict["bin_count"]
        self.bin_min                       = info_dict["bin_min"]
        self.bin_max                       = info_dict["bin_max"]
        self.bin_snps_min                  = info_dict["bin_snps_min"]
        self.bin_snps_max                  = info_dict["bin_snps_max"]
        assert self.bin_max == self.matrixNp.shape[0]
        assert self.bin_max == self.pairwiNp.shape[0]

        bin_width                          = info_dict["bin_width"]
        assert bin_width == self.bin_width
        self.bin_width = bin_width

        self.chromosome_snps               = info_dict["chromosome_snps"]
        chromosome_name                    = meta_dict["chromosome_name"]
        chromosome_order                   = info_dict["chromosome_order"]
        assert self.chromosome_name  == chromosome_name
        assert self.chromosome_order == chromosome_order
        self.chromosome_name               = chromosome_name
        self.chromosome_order              = chromosome_order

        self.chromosome_first_position     = info_dict["chromosome_first_position"]
        self.chromosome_last_position      = info_dict["chromosome_last_position"]
        assert (self.chromosome_first_position // self.bin_width) == self.bin_min
        assert (self.chromosome_last_position  // self.bin_width) == self.bin_max

        metric                             = meta_dict["metric"]
        assert self.metric == metric
        assert metric in METRIC_VALIDS
        self.metric = metric

        type_matrix_counter_name           = meta_dict["type_matrix_counter_name"]
        type_matrix_distance_name          = meta_dict["type_matrix_distance_name"]
        type_pairwise_counter_name         = meta_dict["type_pairwise_counter_name"]
        type_positions_name                = meta_dict["type_positions_name"]

        self.type_matrix_counter           = getattr(np, type_matrix_counter_name)
        self.type_matrix_distance          = getattr(np, type_matrix_distance_name)
        self.type_pairwise_counter         = getattr(np, type_pairwise_counter_name)
        self.type_positions                = getattr(np, type_positions_name)

        assert self.matrixNp.dtype   in (self.type_matrix_counter, self.type_matrix_distance)
        assert self.pairwiNp.dtype   == self.type_pairwise_counter
        assert self.positionNp.dtype == self.type_positions

        type_matrix_counter_max_val        = info_dict["type_matrix_counter_max_val"]
        type_matrix_distance_max_val       = info_dict["type_matrix_distance_max_val"]
        type_pairwise_counter_max_val      = info_dict["type_pairwise_counter_max_val"]
        type_positions_max_val             = info_dict["type_positions_max_val"]

        assert np.iinfo(self.type_matrix_counter  ).max    == type_matrix_counter_max_val
        assert np.finfo(self.type_matrix_distance ).maxexp == type_matrix_distance_max_val
        assert np.iinfo(self.type_pairwise_counter).max    == type_pairwise_counter_max_val
        assert np.iinfo(self.type_positions       ).max    == type_positions_max_val

        self.type_matrix_counter_max_val   = type_matrix_counter_max_val
        self.type_matrix_distance_max_val  = type_matrix_distance_max_val
        self.type_pairwise_counter_max_val = type_pairwise_counter_max_val
        self.type_positions_max_val        = type_positions_max_val
        
        self.sample_count                  = info_dict["sample_count"]
        self.sample_names                  = sample_names
        assert len(self.sample_names) == self.sample_count
        assert self.sample_count      == self.pairwiNp.shape[1]
        if self.alignmentNp.shape[0] != 0:
            assert self.sample_count  == self.alignmentNp.shape[1]
        
        self._is_loaded = True

        print(self)

    def convert_matrix_to_distance(self, metric: str):
        assert metric in METRIC_VALIDS
        assert self.is_loaded

        print(f"chromosome {self.chromosome_name} - calculating distance {metric}")
        dist = np.zeros((self.matrixNp.shape[0], self.matrixNp.shape[1]), dtype=self.type_matrix_distance)

        for binNum in range(self.bin_max):
            # print(f"chromosome {self.chromosome_name} - calculating distance {metric} - bin {binNum}")
            dist[binNum,:] = matrixDistance(self.matrixNp[binNum,:], metric=metric, dtype=self.type_matrix_distance)
        
        return dist

    def matrix_sample(self, sample_name: str, metric=None) -> np.ndarray:
        """
            import reader
            chrom = reader.Chromosome("../data/360_merged_2.50.vcf.gz", 250_000, 0, "SL2.50ch00")
            chrom.load()
            chrom.matrix_sample('TS-1').shape
        """

        if metric is None:
            metric == self.metric

        data = np.zeros((self.bin_max, self.sample_count), dtype=self.matrixNp.dtype)
        
        for binNum in range(self.bin_max):
            if metric == self.metric:
                data[binNum,:] = self.matrix_bin_sample(binNum, sample_name)

            else:
                data[binNum,:] = self.matrix_bin_dist_sample(binNum, sample_name, metric=metric)
        
        return data

    def matrix_bin(self, binNum: int) -> np.ndarray:
        assert binNum <  self.bin_max
        assert binNum >= 0
        return self.matrix[binNum,:]

    def matrix_bin_square(self, binNum: int) -> np.ndarray:
        matrix_bin = self.matrix_bin(binNum)
        mat        = triangleToMatrix(matrix_bin)
        matsum     = mat.sum(axis=0)
        assert sum(matsum - self.pairwiNp[binNum]) == 0
        # np.fill_diagonal(mat, self.pairwiNp[binNum])
        return mat

    def matrix_bin_sample(self, binNum: int, sample_name: str) -> np.ndarray:
        sample_pos = self._sample_pos(sample_name)
        return self.matrix_bin_square(binNum)[sample_pos,:]

    def matrix_bin_dist(self, binNum: int, metric: str = DEFAULT_METRIC, dtype: np.dtype = DEFAULT_DISTANCE_TYPE_MATRIX) -> np.ndarray:
        """
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html

            'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine',
            'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon', 'kulsinski',
            'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
            'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'.
        
            https://stackoverflow.com/questions/35758612/most-efficient-way-to-construct-similarity-matrix
            https://www.dabblingbadger.com/blog/2020/2/27/implementing-euclidean-distance-matrix-calculations-from-scratch-in-python

            import numpy as np
            from scipy.spatial.distance import pdist, squareform

            import reader
            chrom = reader.Chromosome("../data/360_merged_2.50.vcf.gz", 250_000, 0, "SL2.50ch00")
            chrom.load()

            metric=DEFAULT_METRIC
            matrix_bin_matrix = chrom.matrix_bin_matrix(0)[:5,:5]
            matrix_bin_matrix
            1.0/(1.0 + squareform(pdist(matrix_bin_matrix, metric=metric)))
        """

        assert self.metric == METRIC_RAW_NAME, "can only calculate distance from raw tables"

        matrix_bin = self.matrix_bin(binNum)

        if metric == METRIC_RAW_NAME:
            return matrix_bin
        else:
            return matrixDistance(matrix_bin, metric=metric, dtype=dtype)

    def matrix_bin_dist_square(self, binNum: int, metric: str = DEFAULT_METRIC) -> np.ndarray:
        return squareform(self.matrix_bin_dist(binNum, metric=metric))

    def matrix_bin_dist_sample(self, binNum: int, sample_name: str, metric: str = DEFAULT_METRIC) -> np.ndarray:
        sample_pos                    = self._sample_pos(sample_name)
        matrix_bin_matrix_dist        = self.matrix_bin_dist_square(binNum, metric=metric)
        matrix_bin_matrix_dist_sample = matrix_bin_matrix_dist[sample_pos,:]
        return matrix_bin_matrix_dist_sample

    def matrix_bin_dist_sample_square(self, binNum: int, sample_name: str, metric: str = DEFAULT_METRIC) -> np.ndarray:
        return squareform(self.matrix_bin_dist_sample(binNum, sample_name, metric=metric))



class Genome():
    def __init__(self,
            vcf_name              : str,
            bin_width             : int                 = DEFAULT_BIN_SIZE,
            metric                : str                 = DEFAULT_METRIC,

            diff_matrix           : MatrixType          = None,
            IUPAC                 : IUPACType           = None,

            save_alignment        : bool                = DEFAULT_SAVE_ALIGNMENT,
            type_matrix_counter   : np.dtype            = DEFAULT_COUNTER_TYPE_MATRIX,
            type_matrix_distance  : np.dtype            = DEFAULT_DISTANCE_TYPE_MATRIX,
            type_pairwise_counter : np.dtype            = DEFAULT_COUNTER_TYPE_PAIRWISE,
            type_positions        : np.dtype            = DEFAULT_POSITIONS_TYPE
        ):

        self.vcf_name             : str                 = vcf_name
        self.bin_width            : int                 = bin_width
        self.metric               : str                 = metric

        self._diff_matrix         : MatrixType          = diff_matrix
        self._IUPAC               : IUPACType           = IUPAC

        self._save_alignment      : bool                = save_alignment
        self.type_matrix_counter  : np.dtype            = type_matrix_counter
        self.type_matrix_distance : np.dtype            = type_matrix_distance
        self.type_pairwise_counter: np.dtype            = type_pairwise_counter
        self.type_positions       : np.dtype            = type_positions

        self.sample_names         : SampleNamesType     = None
        self.sample_count         : int                 = None

        self.chromosome_names     : ChromosomeNamesType = None
        self.chromosome_count     : int                 = None

        self.genome_bins          : int                 = None
        self.genome_snps          : int                 = None

        self._chromosomes         : ChromosomesType     = None

        assert self.metric in METRIC_VALIDS

    @property
    def filename(self) -> str:
        basefolder =  self.vcf_name + "_ib"
        dirname    = os.path.join(basefolder, str(self.bin_width))
        return f"{dirname}{os.path.sep}{os.path.basename(self.vcf_name)}.npz"

    @property
    def exists(self) -> bool:
        return os.path.exists(self.filename)

    @property
    def loaded(self) -> bool:
        if not self.exists:
            # print(f"genome loaded error - file does not exists {self.filename}", file=sys.stderr)
            return False
        
        if self._chromosomes is None:
            # print("genome loaded error - chromosomes is None", file=sys.stderr)
            return False

        # if len(self._chromosomes) == 0:
        #     return False

        return True

    @property
    def complete(self) -> bool:
        if not self.loaded:
            # print("genome complete error - not loaded", file=sys.stderr)
            return False
        
        for chromosome in self._chromosomes.values():
            if not chromosome.exists:
                # print(f"genome complete error - chromosome does not exists {chromosome.filename}", file=sys.stderr)
                return False

        return True

    def _processVcf(self, threads=DEFAULT_THREADS):
        self._chromosomes     = OrderedDict()
        self.chromosome_names = []
        self.chromosome_count = 0
        self.genome_bins      = 0
        self.genome_snps      = 0

        if self._IUPAC is None:
            self._IUPAC = genIUPAC()

        if self._diff_matrix is None:
            self._diff_matrix = genDiffMatrix()

        assert os.path.exists(self.vcf_name)

        sample_names, sample_count, matrix_size, indexes = Chromosome._processVcf_read_header(self.vcf_name)
        self.sample_names     = sample_names
        self.sample_count     = sample_count

        bgzip                 = BGzip(self.vcf_name)

        self.chromosome_names = bgzip.chromosomes
        self.chromosome_count = len(self.chromosome_names)

        mp.set_start_method('spawn')

        results = []
        with mp.Pool(processes=threads) as pool:
            for chromosome_order, chromosome_name in enumerate(self.chromosome_names):
                chromosome = Chromosome(
                    vcf_name              = self.vcf_name,
                    bin_width             = self.bin_width,
                    metric                = self.metric,
                    chromosome_order      = chromosome_order,
                    chromosome_name       = chromosome_name,
                )

                if chromosome.exists:
                    print(f"chromosome {chromosome_name} already exists")
                    continue

                chromosome = Chromosome(
                    vcf_name              = self.vcf_name,
                    bin_width             = self.bin_width,
                    metric                = METRIC_RAW_NAME,
                    chromosome_order      = chromosome_order,
                    chromosome_name       = chromosome_name,
                )

                if chromosome.exists:
                    print(f"chromosome {chromosome_name} already exists as raw - calculating distance")
                    res = pool.apply_async(
                        Chromosome._calculateDistance,
                        [],
                        {
                            "vcf_name"              : self.vcf_name,
                            "bin_width"             : self.bin_width,
                            "metric"                : self.metric,
                            "chromosome_order"      : chromosome_order,
                            "chromosome_name"       : chromosome_name,
                        }
                    )
                else:
                    print(f"reading chromosome {chromosome_name} from vcf")
                    res = pool.apply_async(
                        Chromosome._processVcf_read_chrom,
                        [],
                        {
                            "vcf_name"              : self.vcf_name,
                            "bin_width"             : self.bin_width,
                            "metric"                : self.metric,

                            "chromosome_order"      : chromosome_order,
                            "chromosome_name"       : chromosome_name,

                            "matrix_size"           : matrix_size,
                            "indexes"               : indexes,

                            "sample_names"          : self.sample_names,
                            "sample_count"          : self.sample_count,

                            "save_alignment"        : self._save_alignment,
                            "IUPAC"                 : self._IUPAC,
                            "diff_matrix"           : self._diff_matrix,

                            "type_matrix_counter"   : self.type_matrix_counter,
                            "type_pairwise_counter" : self.type_pairwise_counter,
                            "type_positions"        : self.type_positions
                        }
                    )
                results.append([False, res, chromosome_order, chromosome_name])
            
            while not all([r[0] for r in results]):
                for resnum, (finished, res, chromosome_order, chromosome_name) in enumerate(results):
                    if finished:
                        continue
                    
                    # if not res.ready():
                    #     continue

                    try:
                        chromosome_bin_count, chromosome_snps = res.get(timeout=1)
                    except mp.TimeoutError:
                        print(f"waiting for {chromosome_name}")
                        continue

                    print(f"got results from {chromosome_name}")

                    results[resnum][0] = True
                    self.genome_bins += chromosome_bin_count
                    self.genome_snps += chromosome_snps

        for resnum, (finished, res, chromosome_order, chromosome_name) in enumerate(results):
            print(f"loading {chromosome_name}")
            
            assert chromosome_name == self.chromosome_names[resnum]
            
            chromosome = Chromosome(
                vcf_name              = self.vcf_name,
                bin_width             = self.bin_width,
                metric                = self.metric,
                
                chromosome_order      = chromosome_order,
                chromosome_name       = chromosome_name,

                type_matrix_counter   = self.type_matrix_counter,
                type_pairwise_counter = self.type_pairwise_counter,
                type_positions        = self.type_positions
            )
            
            if not chromosome.exists:
                raise IOError(f"chromosome database does not exists: {chromosome.filename}")
            else:
                print(f"loading chromosome {chromosome_name} {chromosome.filename}")
            
            chromosome.load()
            
            self._chromosomes[chromosome_name] = chromosome
        
        print("all chromosomes loaded")

    def __repr__(self):
        return str(self)

    def __str__(self):
        res = []
        for k in [
            "vcf_name",
            "bin_width",
            "metric",
            "chromosome_names",
            "chromosome_count",
            "genome_bins",
            "genome_snps",
            "sample_names",
            "sample_count",
            "type_matrix_counter",
            "type_pairwise_counter",
            "type_positions",
        ]:

            v = getattr(self, k)
            s = None
            if   isinstance(v, int):
                s = f"{v:,d}"
            elif isinstance(v, str):
                s = f"{v:s}"
            elif isinstance(v, list):
                s = f"{len(v):,d}"
            elif isinstance(v, np.ndarray):
                s = f"{str(v.shape):s}"
            else:
                s = f"{str(v):s}"
            res.append(f"  {k:.<30s}{s:.>30s}")
        return "\n".join(res)

    def _load_db(self, preload=False):
        if self.complete:
            return

        print(f"{'loading numpy array:':.<32s}{self.filename:.>30s}")
        # filename           = f"{self.vcf_name}.{self.chromosome_order:06d}.{self.chromosome_name}.npz"
        data                       = np.load(self.filename, mmap_mode='r', allow_pickle=False)
        
        self.sample_names          = data["sample_names"].tolist()
        self.chromosome_names      = data["chromosome_names"].tolist()

        info_namesNp               = data['info_names']
        info_valuesNp              = data['info_values']

        meta_namesNp               = data['meta_names']
        meta_valuesNp              = data['meta_values']

        info_names                 = info_namesNp.tolist()
        info_values                = info_valuesNp.tolist()
        info_values                = [int(v) for v in info_values]
        info_dict                  = {info_names[k]: info_values[k]  for k in range(len(info_names))}
        # print(info_dict)

        meta_names                 = meta_namesNp.tolist()
        meta_values                = meta_valuesNp.tolist()
        meta_dict                  = {meta_names[k]: meta_values[k]  for k in range(len(meta_names))}
        # print(meta_dict)

        vcf_name                   = meta_dict["vcf_name"]
        assert os.path.basename(vcf_name) == os.path.basename(self.vcf_name)

        type_matrix_counter_name   = meta_dict["type_matrix_counter_name"]
        type_matrix_distance_name  = meta_dict["type_matrix_distance_name"]
        type_pairwise_counter_name = meta_dict["type_pairwise_counter_name"]
        type_positions_name        = meta_dict["type_positions_name"]

        self.type_matrix_counter   = getattr(np, type_matrix_counter_name)
        self.type_matrix_distance  = getattr(np, type_matrix_distance_name)
        self.type_pairwise_counter = getattr(np, type_pairwise_counter_name)
        self.type_positions        = getattr(np, type_positions_name)

        self.bin_width             = info_dict["bin_width"]
        self.chromosome_count      = info_dict["chromosome_count"]
        self.sample_count          = info_dict["sample_count"]
        self.genome_bins           = info_dict["genome_bins"]
        self.genome_snps           = info_dict["genome_snps"]

        assert self.chromosome_count == len(self.chromosome_names)
        assert self.sample_count     == len(self.sample_names)

        print(self)

        self._chromosomes = OrderedDict()

        if preload:
            chromosome_bins   = 0
            chromosome_snps   = 0
            
            for chromosome_name in self.chromosome_names:
                chromosome       = self.get_chromosome(chromosome_name)

                chromosome_bins += chromosome.bin_count
                chromosome_snps += chromosome.chromosome_snps

            assert self.genome_bins == chromosome_bins, f"self.genome_bins {self.genome_bins} == {chromosome_bins} chromosome_bins"
            assert self.genome_snps == chromosome_snps, f"self.genome_snps {self.genome_snps} == {chromosome_snps} chromosome_snps"

    def _get_infos(self):
        names  = ["bin_width"   , "chromosome_count"   , "sample_count"   , "genome_bins"   , "genome_snps"   ]
        vals =   [getattr(self,n) for n in names]

        return (names, vals)

    def _get_meta(self):
        type_matrix_counter_name   = self.type_matrix_counter.__name__
        type_matrix_distance_name  = self.type_matrix_distance.__name__
        type_pairwise_counter_name = self.type_pairwise_counter.__name__
        type_positions_name        = self.type_positions.__name__

        names  = ["vcf_name", "metric"]
        vals   = [getattr(self,n) for n in names]

        names += [
            "type_matrix_counter_name"  , "type_matrix_distance_name",
            "type_pairwise_counter_name", "type_positions_name"
        ]
        vals  += [
            type_matrix_counter_name  , type_matrix_distance_name,
            type_pairwise_counter_name, type_positions_name
        ]

        return names, vals

    def todict(self):
        info_names,info_values     = self._get_infos()
        meta_names,meta_values     = self._get_meta()

        data =      {k:v for k,v in zip(info_names, info_values)}
        data.update({k:v for k,v in zip(meta_names, meta_values)})
        
        data['sample_names'    ] = self.sample_names
        data['chromosome_names'] = self.chromosome_names

        return data

    def save(self):
        print(f"{'saving numpy array:':.<32s}{self.filename:.>30s}")
        print(self)

        sample_namesNp             = np.array(self.sample_names    , np.unicode_)
        chromosome_namesNp         = np.array(self.chromosome_names, np.unicode_)

        info_names,info_values     = self._get_infos()
        meta_names,meta_values     = self._get_meta()

        info_namesNp               = np.array(info_names , np.unicode_)
        info_valuesNp              = np.array(info_values, np.int64   )

        meta_namesNp               = np.array(meta_names , np.unicode_)
        meta_valuesNp              = np.array(meta_values, np.unicode_)

        np.savez_compressed(self.filename,
            sample_names     = sample_namesNp,
            chromosome_names = chromosome_namesNp,
            info_names       = info_namesNp,
            info_values      = info_valuesNp,
            meta_names       = meta_namesNp,
            meta_values      = meta_valuesNp
        )

    def load(self, preload=False, threads=DEFAULT_THREADS):
        if self.exists:
            print(f"database exists {self.filename}")
            self._load_db(preload=preload)

        else:
            print(f"database does not exists {self.filename}. reading vcf")
            self._processVcf(threads=threads)

            self.save()

            if not self.complete:
                raise IOError("not complete. not able to load the database")

    def get_chromosome(self, chromosome_name: str) -> Chromosome:
        # chromosome_position = self.chromosome_names.index(chromosome_name)
        assert chromosome_name in self.chromosome_names, f"invalid chromosome name {chromosome_name}: {','.join(self.chromosome_names)}"
        
        if chromosome_name not in self._chromosomes:
            chromosome_order = self.chromosome_names.index(chromosome_name)
            chromosome       = Chromosome(
                vcf_name              = self.vcf_name,
                bin_width             = self.bin_width,
                metric                = self.metric,

                chromosome_order      = chromosome_order,
                chromosome_name       = chromosome_name,

                type_matrix_counter   = self.type_matrix_counter,
                type_matrix_distance  = self.type_matrix_distance,
                type_pairwise_counter = self.type_pairwise_counter,
                type_positions        = self.type_positions 
            )
            
            if not chromosome.exists:
                raise IOError(f"chromosome database file {chromosome.filename} does not exists")

            chromosome.load()
            
            assert os.path.basename(chromosome.vcf_name) == os.path.basename(self.vcf_name)
            assert chromosome.bin_width                  == self.bin_width
            assert chromosome.chromosome_order           == chromosome_order
            assert chromosome.chromosome_name            == chromosome_name
            assert chromosome.metric                     == self.metric
            assert chromosome.sample_names               == self.sample_names
            assert chromosome.type_matrix_counter        == self.type_matrix_counter
            assert chromosome.type_matrix_distance       == self.type_matrix_distance
            assert chromosome.type_pairwise_counter      == self.type_pairwise_counter
            assert chromosome.type_positions             == self.type_positions

            self._chromosomes[chromosome_name] = chromosome

        return self._chromosomes[chromosome_name]



class Genomes():
    def __init__(self, folder_name: str, verbose=False):
        self.folder_name     : str        = folder_name
        self._genomes                     = None

        self.curr_genome_name: str        = None
        self.curr_genome_binw: int        = None
        self.curr_genome_metr: str        = None
        self.curr_genome     : Genome     = None
        self.curr_chrom_name : str        = None
        self.curr_chrom      : Chromosome = None

        self.update(verbose=verbose)

        if DEBUG:
            print("genomes    ", self.genomes)
            # print("genome_info", self.genomes[0], self.genome_info(self.genomes[0]))
            genomes     = self.genomes
            genome_name = genomes[0]
            print("genome_name", genome_name)
            bin_widths  = self.bin_widths(genome_name)
            print("bin_widths ", bin_widths)
            bin_width    = bin_widths[0]
            # print("bin_width_info", self.bin_width_info(genome_name, bin_width)
            metrics     = self.metrics(genome_name, bin_width)
            print("metrics    ", metrics)
            metric      = metrics[0]
            print("metric     ", metric)
            # print("metric_info", self.metric_info(genome_name, bin_width, metric)
            chroms = self.chromosome_names(genome_name, bin_width, metric)
            print("chroms     ", chroms)

            self.load_genome(genome_name, bin_width, metric)
            print("chromosomes", self.chromosomes)
            chromosome = self.chromosomes[0]
            print("chromosome " , chromosome)
            self.load_chromosome(genome_name, bin_width, metric, chromosome)

    @property
    def genomes(self) -> typing.List[str]:
        return list(self._genomes.keys())

    @property
    def genome(self) -> Genome:
        return self.curr_genome

    @property
    def chromosomes(self) -> ChromosomeNamesType:
        if self.genome is None:
            return None
        return self.genome.chromosome_names

    @property
    def chromosome(self) -> Chromosome:
        return self.curr_chrom

    def genome_info(self, genome_name: str):
        assert genome_name in self.genomes
        return self._genomes[genome_name]

    def bin_widths(self, genome_name: str) -> typing.List[str]:
        return list(self.genome_info(genome_name)["bin_widths"].keys())

    def bin_width_info(self, genome_name: str, bin_width: int):
        assert bin_width   in self.bin_widths(genome_name)
        return self.genome_info(genome_name)["bin_widths"][bin_width]

    def metrics(self, genome_name: str, bin_width: int) -> typing.List[str]:
        return list(self.bin_width_info(genome_name, bin_width)["metrics"].keys())

    def metric_info(self, genome_name: str, bin_width: int, metric: str) -> typing.Dict[str, typing.Any]:
        assert metric      in self.metrics(genome_name, bin_width)
        return self.bin_width_info(genome_name, bin_width)["metrics"][metric]

    def chromosome_names(self, genome_name: str, bin_width: int, metric: str) -> typing.List[typing.Tuple[int, str]]:
        metric_info      = self.metric_info(genome_name, bin_width, metric)
        chromosome_names = [(m["chromosome_pos"], m["chromosome_name"]) for m in metric_info]
        return chromosome_names

    def update(self, verbose=False):
        self._genomes = Genomes.listProjects(self.folder_name, verbose=verbose)

    def load_genome(self, genome_name: str, bin_width: int, metric: str) -> Genome:
        if not (
                genome_name == self.curr_genome_name and
                bin_width   == self.curr_genome_binw and
                metric      == self.curr_genome_metr
        ):
            self.curr_genome_name: str        = None
            self.curr_genome_binw: int        = None
            self.curr_genome_metr: str        = None
            self.curr_genome     : Genome     = None
            self.curr_chrom_name : str        = None
            self.curr_chrom      : Chromosome = None

            # projects[dbname]["bin_widths"][bin_width]["metrics"][metric]
            if genome_name   not in self.genomes:
                raise ValueError(f"no such database {genome_name}. {','.join(self.genomes)}")

            if bin_width not in self.bin_widths(genome_name):
                raise ValueError(f"no such bin width {bin_width} for database {genome_name}: {self.bin_widths(genome_name)}")

            if metric    not in self.metrics(genome_name, bin_width):
                raise ValueError(f"no such metric {metric} for database {genome_name} bin width {bin_width}: {self.metrics(genome_name, bin_width)}")

            # print(self.genomes[genome_name])
            project_path = self.genome_info(genome_name)["projectpath"]
            print("loading project_path", project_path)

            genome = Genome(
                vcf_name  = project_path,
                bin_width = bin_width,
                metric    = metric
            )

            print("genome.filename", genome.filename)
            
            assert genome.exists
            
            # print(genome)
            
            genome.load()

            assert genome.loaded
            assert genome.complete

            self.curr_genome_name: str        = genome_name
            self.curr_genome_binw: int        = bin_width
            self.curr_genome_metr: str        = metric
            self.curr_genome     : Genome     = genome
            self.curr_chrom_name : str        = None
            self.curr_chrom      : Chromosome = None

        return self.curr_genome

    def load_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str) -> Chromosome:
        genome = self.load_genome(genome_name, bin_width, metric)
        
        if chromosome_name != self.curr_chrom_name:
            chromosome           = genome.get_chromosome(chromosome_name)
            self.curr_chrom_name = chromosome_name
            self.curr_chrom      = chromosome
        
        return self.curr_chrom

    @staticmethod
    def listProjects(folder_name, verbose=False):
        basepath  = os.path.abspath(folder_name)
        filenames = list(iglob(os.path.join(basepath,"**/*.npz"), recursive=True))
        projects  = OrderedDict()

        for filepath in filenames:
            filedir     = filepath[len(basepath)+1:]
            filename    = os.path.basename(filedir)
            filefolder  = filedir[:-1*len(filename)-1]
            dbname      = filefolder.split(os.path.sep)[0][:-3]

            for ext in ['.vcf.bgz', '.vcf.gz', '.vcf']:
                dbname = dbname[:-1*len(ext)] if dbname.endswith(ext) else dbname
            dbname = dbname.replace("_", " ")

            datafolder, bin_width, metric   = None, None, None
            chromosome_pos, chromosome_name = None, None
            folder_parts                    = filefolder.strip(os.path.sep).split(os.path.sep)
            
            if filename.startswith("ib_"): # data
                try:
                    projectfolder                   = os.path.join(*folder_parts[:-2])
                    assert projectfolder.endswith("_ib")
                    projectfolder                   = projectfolder[:-3]
                    datafolder, bin_width, metric   = folder_parts[-3:]
                    fileparts                       = filename[3:-4].split('.')
                    chromosome_pos, chromosome_name = fileparts[0], ".".join(fileparts[1:])
                    chromosome_pos = int(chromosome_pos)
                except:
                    print(f"invalid folder {filefolder}")
                    continue

            else: # root
                try:
                    projectfolder        = os.path.join(*folder_parts[:-1])
                    assert projectfolder.endswith("_ib")
                    projectfolder        = projectfolder[:-3]
                    datafolder, bin_width = filefolder.strip(os.path.sep).split(os.path.sep)[-2:]
                except:
                    print(f"invalid folder {filefolder}")
                    continue

            projectpath = os.path.join(basepath, projectfolder)
            bin_width    = int(bin_width)
            # print(f"filepath {filepath} filedir {filedir} filename {filename} filefolder {filefolder} datafolder {datafolder} bin_width {bin_width} dbname {dbname}")

            if dbname not in projects:
                projects[dbname] = {
                    "dbname"       : dbname,
                    "datafolder"   : datafolder,
                    "projectfolder": projectfolder,
                    "projectpath"  : projectpath,
                    "bin_widths"    : OrderedDict()
                }

            if bin_width not in projects[dbname]["bin_widths"]:
                projects[dbname]["bin_widths"][bin_width] = {
                    "dbname"       : dbname,
                    "datafolder"   : datafolder,
                    "projectfolder": projectfolder,
                    "projectpath"  : projectpath,
                    "bin_width"     : bin_width,
                    "metrics"      : OrderedDict()
                }

            if metric is None: # root
                projects[dbname]["bin_widths"][bin_width].update({
                    "filepath"  : filepath,
                    "filedir"   : filedir,
                    "filename"  : filename,
                    "filefolder": filefolder,
                })
            
            else: #data
                if metric not in projects[dbname]["bin_widths"][bin_width]["metrics"]:
                    projects[dbname]["bin_widths"][bin_width]["metrics"][metric] = []
                
                projects[dbname]["bin_widths"][bin_width]["metrics"][metric].append({
                    "filepath"       : filepath,
                    "filedir"        : filedir,
                    "filename"       : filename,
                    "filefolder"     : filefolder,
                    "projectfolder"  : projectfolder,
                    "projectpath"    : projectpath,
                    "datafolder"     : datafolder,
                    "bin_width"      : bin_width,
                    "metric"         : metric,
                    "dbname"         : dbname,
                    "chromosome_pos" : chromosome_pos,
                    "chromosome_name": chromosome_name
                })

                projects[dbname]["bin_widths"][bin_width]["metrics"][metric].sort(key=lambda v: v["chromosome_pos"])

        assert len(projects) > 0

        if verbose:
            for dbname, dbdata in projects.items():
                print(f"{'database':23s}: {dbname}")
                
                for dataname, datavalue in dbdata.items():
                    if dataname == 'bin_widths':
                        continue
                    print(f"  {dataname:21s}: {datavalue}")
                
                for bin_width, binvalues in dbdata['bin_widths'].items():
                    print(f"  {'bin_width':21s}: {bin_width}")
                    
                    for binkey, binvalue in binvalues.items():
                        if binkey == 'metrics':
                            continue
                        print(f"    {binkey:19s}: {binvalue}")

                    for metrickey, metricvalues in binvalues['metrics'].items():
                        print(f"    {'metric':19s}: {metrickey}")
                        
                        for metric_pos, metric in enumerate(metricvalues):
                            print(f"      {'file #':17s}: {metric_pos}")
                            for metrickey, metricvalue in metric.items():
                                print(f"        {metrickey:15s}: {metricvalue}")

        return projects


class BGzip():
    """
        http://www.htslib.org/doc/bgzip.html#:~:text=GZI%20FORMAT,in%20the%20uncompressed%20data%20stream.

        The index format is a binary file listing pairs of compressed and
        uncompressed offsets in a BGZF file. Each compressed offset points to
        the start of a BGZF block. The uncompressed offset is the corresponding
        location in the uncompressed data stream.

        All values are stored as little-endian 64-bit unsigned integers.
        
        The file contents are:
            uint64_t number_entries

        followed by number_entries pairs of:
            uint64_t compressed_offset
            uint64_t uncompressed_offset
    """

    """
        rm 360_merged_2.50.vcf.gz.gzi
        rm 360_merged_2.50.vcf.gz.gzj

        import reader
        b = reader.BGzip("../data/360_merged_2.50.vcf.gz")
        b.chromosomes
        for line in b.get_chromosome('SL2.50ch00'):
            print(line)
        for line in b.get_chromosome('SL2.50ch11'):
            print(line)
    """

    def __init__(self, gzip_file: str):
        self.gzip_file = gzip_file
        self.gzi_file  = gzip_file + ".gzi"
        self.gzj_file  = gzip_file + ".gzj"

        assert os.path.exists(self.gzip_file)

        if not os.path.exists(self.gzi_file):
            print(f"index file {self.gzi_file} does not exists. creating")
            print(os.system(f"bgzip -r {self.gzip_file}"))
            print(f"index file created")

        self._data = OrderedDict()

        if os.path.exists(self.gzj_file):
            self._load()
        else:
            print(f"index file {self.gzj_file} does not exists. creating")
            self._parse_gzi()
            self._save()

    @property
    def chromosomes(self) -> ChromosomeNamesType:
        chroms = [k for k in self._data.keys()]
        chroms.sort(key=lambda x: self._data[x]["entry_num"])
        return chroms

    def _save(self):
        print(f"saving gzj to {self.gzj_file}")
        json.dump(self._data, open(self.gzj_file, 'wt'), indent=1)

    def _load(self):
        print(f"reading gzj from {self.gzj_file}")
        self._data = json.load(open(self.gzj_file, 'rt'))

    def _parse_gzi(self):
        print(f"reading gzi from {self.gzi_file}")

        with open(self.gzip_file, 'rb') as gzip_fhd:
            with open(self.gzi_file, 'rb') as fhd:
                number_entries_fmt  = "<Q"
                number_entries_size = struct.calcsize(number_entries_fmt)
                
                offsets_fmt         = "<QQ"
                offsets_size        = struct.calcsize(offsets_fmt)
                
                number_entries      = struct.unpack(number_entries_fmt, fhd.read(number_entries_size))[0]
                print("number_entries", number_entries)
                
                previous_compressed_offset   = 0
                previous_uncompressed_offset = 0
                current_compressed_offset    = 0
                current_uncompressed_offset  = 0
                for entry_num in range(number_entries):
                    next_compressed_offset, next_uncompressed_offset = struct.unpack(offsets_fmt, fhd.read(offsets_size))
                    previous_compressed_size   = current_compressed_offset   - previous_compressed_offset
                    previous_uncompressed_size = current_uncompressed_offset - previous_uncompressed_offset
                    current_compressed_size    = next_compressed_offset      - current_compressed_offset
                    current_uncompressed_size  = next_uncompressed_offset    - current_uncompressed_offset

                    # print("current_compressed_offset", current_compressed_offset, "compressed_offset", compressed_offset, "current_uncompressed_offset", current_uncompressed_offset, "uncompressed_offset", uncompressed_offset, "compressed_size", compressed_size, "uncompressed_size", uncompressed_size)
                    
                    self._get_first_contig(
                        gzip_fhd,
                        entry_num,
                        
                        previous_compressed_offset,
                        previous_uncompressed_offset,
                        previous_compressed_size,
                        previous_uncompressed_size,

                        current_compressed_offset,
                        current_uncompressed_offset,
                        current_compressed_size,
                        current_uncompressed_size,
                    )
                    
                    previous_compressed_offset   = current_compressed_offset
                    previous_uncompressed_offset = current_uncompressed_offset
                    current_compressed_offset    = next_compressed_offset
                    current_uncompressed_offset  = next_uncompressed_offset
    
    def _get_first_contig(self,
            gzip_fhd : typing.IO,
            entry_num: int,
            
            previous_compressed_offset  : int,
            previous_uncompressed_offset: int,
            previous_compressed_size    : int,
            previous_uncompressed_size  : int,

            current_compressed_offset   : int,
            current_uncompressed_offset : int,
            current_compressed_size     : int,
            current_uncompressed_size   : int,
        ):

        gzip_fhd.seek(current_compressed_offset)

        compressed_block        = gzip_fhd.read(current_compressed_size)
        uncompressed_block      = gzip.decompress(compressed_block)
        assert len(uncompressed_block) == current_uncompressed_size
        uncompressed_block_text = uncompressed_block.decode()
        
        firstNewLine = 0
        first_tab     = 0
        chrom_name    = None
        while True:
            if firstNewLine >= len(uncompressed_block_text):
                raise ValueError(f"firstNewLine {firstNewLine} >= {len(uncompressed_block_text)} len(uncompressed_block_text)")

            firstNewLine  = uncompressed_block_text.find("\n", firstNewLine+1)
            if firstNewLine == -1:
                print("no firstNewLine")
                continue
            
            secondNewLine = uncompressed_block_text.find("\n", firstNewLine+1)
            if secondNewLine == -1:
                print("no secondNewLine")
                firstNewLine += 1
                continue
            
            firstLine     = uncompressed_block_text[firstNewLine+1:secondNewLine]
            if len(firstLine) == 0:
                print("no firstLine", firstNewLine, secondNewLine)
                firstNewLine += 1
                continue
            
            if firstLine[0] == "#":
                print("skipping header", firstLine)
                firstNewLine += 1
                continue
            
            first_tab      = firstLine.find("\t")
            if first_tab == -1:
                print("no first_tab", firstLine)
                continue
            
            chrom_name     = firstLine[:first_tab]
            if len(chrom_name) == 0:
                print("no chrom_name", firstNewLine, first_tab, firstLine)
                chrom_name = None
                continue

            break

        if chrom_name is None:
            raise ValueError("No chromosome name found")

        if chrom_name not in self._data:
            print(f"NEW CHROMOSOME")
            print(f"  chrom_name                   {chrom_name:>15s}")
            print(f"  entry_num                    {entry_num:15,d}")

            print(f"  previous_compressed_offset   {previous_compressed_offset:15,d}")
            print(f"  previous_uncompressed_offset {previous_uncompressed_offset:15,d}")
            print(f"  previous_compressed_size     {previous_compressed_size:15,d}")
            print(f"  previous_uncompressed_size   {previous_uncompressed_size:15,d}")

            print(f"  current_compressed_offset    {current_compressed_offset:15,d}")
            print(f"  current_uncompressed_offset  {current_uncompressed_offset:15,d}")
            print(f"  current_compressed_size      {current_compressed_size:15,d}")
            print(f"  current_uncompressed_size    {current_uncompressed_size:15,d}")

            self._data[chrom_name] = {
                "entry_num"                      : entry_num,
                
                "previous_compressed_offset"     : previous_compressed_offset,
                "previous_uncompressed_offset"   : previous_uncompressed_offset,
                "previous_compressed_size"       : previous_compressed_size,
                "previous_uncompressed_size"     : previous_uncompressed_size,

                "current_compressed_offset"      : current_compressed_offset,
                "current_uncompressed_offset"    : current_uncompressed_offset,
                "current_compressed_size"        : current_compressed_size,
                "current_uncompressed_size"      : current_uncompressed_size,
            }

    def get_chromosome(self, chrom_name: str) -> typing.Generator[str, None, None]:
        if chrom_name not in self._data:
            raise ValueError(f"chromosome {chrom_name} does not exists: {list(self._data.keys())}")
        
        chromosome_matrix = self._data[chrom_name]
        file_pos          = chromosome_matrix["previous_compressed_offset"]
        found_chrom       = False

        with open(self.gzip_file, 'rb') as fhd:
            fhd.seek(file_pos)
            ghd = gzip.open(fhd, 'rt')

            for linenum, line in enumerate(ghd):
                if len(line) == 0:
                    continue

                if line[0] == "#":
                    continue

                if file_pos != 0: #first chromosome
                    if linenum == 0: #possible incomplete line
                        # print('INCOMPLETE LINE', line)
                        continue

                first_tab = line.find("\t")
                if first_tab == -1:
                    print('NO TAB', line)
                    continue

                chrom = line[:first_tab]
                if chrom != chrom_name:
                    if found_chrom:
                        # print('WRONG CHROMOSOME TAIL', chrom, line)
                        break
                    else:
                        # print('WRONG CHROMOSOME HEAD', chrom, line)
                        continue

                found_chrom = True
                yield line

                # break



def openFile(filename, mode):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def calculateMatrixSize(sample_count: int) -> int:
    """
        2 - 1   2 - 3
        - -     1 -
        1 -     2 3

        3 - 3   3 - 6
        - - -   1 - -
        1 - -   2 3 -
        2 3 -   4 5 6

        4 - 6    4 - 10
        - - - -  1 - - -
        1 - - -  2 3 - -
        2 3 - -  4 5 6 -
        4 5 6 -  7 8 9 10

        5 - 10      5 - 15
        - - - - -   1  -  -  -  -
        1 - - - -   2  3  -  -  -
        2 3 - - -   4  5  6  -  -
        4 5 6 - -   7  8  9 10  -
        7 8 9 10 - 11 12 13 14 15

        f = lambda sample_count: sum([x for x in range(sample_count)])
        f = lambda sample_count: sum([x for x in range(sample_count)]) + sample_count
    """
    return sum([x for x in range(sample_count)])# + sample_count

def triangleToIndex(size: int) -> typing.Tuple[int, TriangleIndexType]:
    """
        1
        2 3
        4 5 6
        7 8 9 10

        1,1 =  1

        2,1 =  2
        2,2 =  3

        3,1 =  4
        3,2 =  5
        3,3 =  6

        4,1 =  7
        4,2 =  8
        4,3 =  9
        4,4 = 10

        ========

        creating a lower triangle first and a list of counts in the end

        for 4 samples:
        1: count sample 1    -                    -                   -
        2: 1,2               3: count sample 2    -                   -
        4: 1,3               5: 2,3               6: count sample 3   -
        7: 1,4               8: 2,4               9: 3,4              10: count sample 4

        -                    -                    -          -
        1: 1,2               -                    -          -
        2: 1,3               3: 2,3               -          -
        4: 1,4               5: 2,4               6: 3,4     -
        7    8    9    10
        [ cs1, cs2, cs3, cs4 ]

        1   2   3   4   5   6   7   8   9   10
        1,2 1,2 2,3 1,4 2,4 3,4 cs1 cs2 cs3 cs4
    """

    """
        b = np.tril_indices(size)
        # (array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3]), array([0, 0, 1, 0, 1, 2, 0, 1, 2, 3]))
        b = np.tril_indices(4, -1)
        # (array([1, 2, 2, 3, 3, 3]), array([0, 0, 1, 0, 1, 2]))
        a = np.arange(16).reshape(4, 4)
        a[b]
        # array([ 4,  8,  9, 12, 13, 14])
        bp = list(zip(b[0].tolist(), b[1].tolist()))
        bp
        # [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
    """

    indexes = OrderedDict()
    # b       = np.tril_indices(size, -1)
    b       = np.triu_indices(size,  1)
    l       = len(b[0].tolist())
    bp      = list(zip(b[0].tolist(), b[1].tolist()))
    # print("b", b)
    # print("bp", bp)
    
    for bl, bc in enumerate(bp):
        # print("index", bc, bl)
        indexes[bc] = bl

    # for bl in range(size):
    #     indexes[(bl, bl)] = len(bp) + bl
    
    return b, l, indexes

def triangleToMatrix(tri_array) -> np.ndarray:
    """
        ```python
        a = np.array([[1,2,3],[4,5,6],[7,8,9]])
        a
        #
        #array([[1, 2, 3],
        #       [4, 5, 6],
        #       [7, 8, 9]])
        #
        size = a.shape[0]
        i, j = np.tril_indices(size, -1)
        i
        #
        #array([1, 2, 2])
        #
        j
        #
        #array([0, 0, 1])
        #
        a[i,j]
        #
        #array([4, 7, 8])
        #
        b = a[i,j]
        b
        #
        #array([4, 7, 8])
        #
        M = np.zeros([size,size], a.dtype)
        M
        #
        #array([[0, 0, 0],
        #       [0, 0, 0],
        #       [0, 0, 0]])
        #
        M[i,j] = b
        M[j,i] = b
        M
        #
        #array([[0, 4, 7],
        #       [4, 0, 8],
        #       [7, 8, 0]])
        #
        ```
    """

    """
        ```python
        import reader
        import numpy as np
        (i,j), l, od = reader.triangleToIndex(5)
        i
        j
        l
        od
        # OrderedDict([((0, 1), 0), ((0, 2), 1), ((0, 3), 2), ((0, 4), 3), ((1, 2), 4), ((1, 3), 5), ((1, 4), 6), ((2, 3), 7), ((2, 4), 8), ((3, 4), 9)])
        for i in range(5):
            for j in range(i+1, 5):
                od[(i,j)]

        # 0
        # 1
        # 2
        # 3
        # 4
        # 5
        # 6
        # 7
        # 8
        # 9
        ln = np.array([0,1,2,3,4,5,6,7,8,9])
        M = reader.triangleToMatrix(5, ln)
        M
        # array(
        # [[0, 0, 1, 2, 3],
        # [0, 0, 4, 5, 6],
        # [1, 4, 0, 7, 8],
        # [2, 5, 7, 0, 9],
        # [3, 6, 8, 9, 0]])
        ```
    """
    # (i,j), l, _ = triangleToIndex(size)
    # assert len(tri_array) == l

    # M           = np.zeros((size, size), tri_array.dtype)

    # M[i, j] = tri_array
    # M[j, i] = tri_array
    M = squareform(tri_array)
    
    return M

def matrixDistance(matrix: np.ndarray, metric=DEFAULT_METRIC, dtype=np.float64) -> np.ndarray:
    """
        import numpy as np
        from scipy.spatial.distance import pdist, squareform
        ln = np.array([0,1,2,3,4,5,6,7,8,9])
        1/(1 + squareform(pdist(squareform(ln), metric=DEFAULT_METRIC)))
    """
    # print("matrixDistance", matrix.shape, matrix)
    if   len(matrix.shape) == 1:
        # return 1.0/(1.0 + squareform(pdist(squareform(matrix), metric=metric)))
        return (1.0/(1.0 + pdist(squareform(matrix), metric=metric))).astype(dtype)

    elif len(matrix.shape) == 2:
        # return 1.0/(1.0 + squareform(pdist(matrix, metric=metric)))
        return 1.0/(1.0 + pdist(matrix, metric=metric)).astype(dtype)

    else:
        raise ValueError("only 1 and 2 dimension arrays can be used")

def genDiffMatrix(alphabet: typing.List[str]=list(range(4))) -> MatrixType:
    """
        diff_matrixSymetricalHomoExtra = {
            # HomoHomo = 3
            # HeteHete = 2
            # HomoHete = 1
            '00': { '00': 3, '01': 1, '02': 1, '10': 1, '11': 0, '20': 1, '22': 0 },
            '01': {          '01': 2, '02': 1, '10': 2, '11': 1, '20': 1, '22': 0 },
            '02': {                   '02': 2, '10': 2, '11': 0, '20': 2, '22': 1 },
            '10': {                            '10': 2, '11': 1, '20': 1, '22': 0 },
            '11': {                                     '11': 0, '20': 0, '22': 0 },
            '20': {                                              '20': 2, '22': 1 },
            '22': {                                                       '22': 3 }
        }
    """

    diff_matrixSymetricalHomoExtra = OrderedDict()

    for n1v in alphabet:
        for n2v in alphabet:
            # nv = n1v + n2v
            nv = (n1v,n2v)
            for o1v in alphabet:
                for o2v in alphabet:
                    # ov = o1v + o2v
                    ov = (o1v,o2v)
                    # k  = (nv,ov) if nv <= ov else (ov,nv)
                    # k  = f"{n1v}|{n2v}|{o1v}|{o2v}" if nv <= ov else f"{o1v}|{o2v}|{n1v}|{n2v}"
                    k = (nv,ov) if nv <= ov else (ov,nv)
                    if k not in diff_matrixSymetricalHomoExtra:
                        val = 0
                        if n1v == o1v:
                            val += 1
                        if n1v == o2v:
                            val += 1
                        if n2v == o1v:
                            val += 1
                        if n2v == o2v:
                            val += 1
                        diff_matrixSymetricalHomoExtra[k] = val
        # print(diff_matrixSymetricalHomoExtra)

    # diff_matrixSymetricalHomoExtra_simple = {}
    # for k1, v1 in diff_matrixSymetricalHomoExtra.items():
    #     for k2, v2 in v1.items():
    #         diff_matrixSymetricalHomoExtra_simple[(k1,k2)] = v2

    return diff_matrixSymetricalHomoExtra

def genIUPAC() -> IUPACType:
    """
        IUPAC nucleotide code	Base
        A	Adenine
        C	Cytosine
        G	Guanine
        T (or U)	Thymine (or Uracil)
        
        M	A or C
        R	A or G
        W	A or T
        S	C or G
        Y	C or T
        K	G or T

        B	C or G or T
        D	A or G or T
        H	A or C or T
        V	A or C or G
        
        N	any base
        
        . or -	gap
    """

    _IUPAC = [
        [('A', 'A'), 'A'],
        [('C', 'C'), 'C'],
        [('G', 'G'), 'G'],
        [('T', 'T'), 'T'],

        [('A', 'C'), 'M'],
        [('A', 'G'), 'R'],
        [('A', 'T'), 'W'],
        [('C', 'G'), 'S'],
        [('C', 'T'), 'Y'],
        [('G', 'T'), 'K'],
    ]

    IUPAC  = {}
    for pair, val in _IUPAC:
        IUPAC[frozenset(pair)] = val

    return IUPAC



def main():
    genome = Genome(
        sys.argv[1],
        bin_width             = DEFAULT_BIN_SIZE,
        metric                = DEFAULT_METRIC,

        diff_matrix           = genDiffMatrix(alphabet=list(range(4))),
        IUPAC                 = genIUPAC(),

        save_alignment        = DEFAULT_SAVE_ALIGNMENT,
        type_matrix_counter   = DEFAULT_COUNTER_TYPE_MATRIX,
        type_matrix_distance  = DEFAULT_DISTANCE_TYPE_MATRIX,
        type_pairwise_counter = DEFAULT_COUNTER_TYPE_PAIRWISE,
        type_positions        = DEFAULT_POSITIONS_TYPE
    )
    genome.load(threads=6)
    # genome.load(threads=DEFAULT_THREADS if not DEBUG else 1)


if __name__ == "__main__":
    main()
