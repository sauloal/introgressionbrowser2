# Introgression Browser v2.0

## About

Complete rewrite from the group up of Introgression Browser

## Output

```txt
data/
|-- [6.9G]  150_VCFs_2.50.tar.gz
|-- [1022M]  360_merged_2.50.vcf.gz
`-- [594M]  360_merged_2.50.vcf.gz_ib
    `-- [594M]  250000
        |-- [2.6K]  360_merged_2.50.vcf.gz.npz
        |-- [387M]  jaccard
        |   |-- [8.7M]  ib_000000.SL2.50ch00.npz
        |   |-- [ 44M]  ib_000001.SL2.50ch01.npz
        |   |-- [ 26M]  ib_000002.SL2.50ch02.npz
        |   |-- [ 35M]  ib_000003.SL2.50ch03.npz
        |   |-- [ 30M]  ib_000004.SL2.50ch04.npz
        |   |-- [ 33M]  ib_000005.SL2.50ch05.npz
        |   |-- [ 25M]  ib_000006.SL2.50ch06.npz
        |   |-- [ 35M]  ib_000007.SL2.50ch07.npz
        |   |-- [ 33M]  ib_000008.SL2.50ch08.npz
        |   |-- [ 37M]  ib_000009.SL2.50ch09.npz
        |   |-- [ 32M]  ib_000010.SL2.50ch10.npz
        |   |-- [ 21M]  ib_000011.SL2.50ch11.npz
        |   `-- [ 28M]  ib_000012.SL2.50ch12.npz
        `-- [207M]  RAW
            |-- [5.6M]  ib_000000.SL2.50ch00.npz
            |-- [ 24M]  ib_000001.SL2.50ch01.npz
            |-- [ 13M]  ib_000002.SL2.50ch02.npz
            |-- [ 16M]  ib_000003.SL2.50ch03.npz
            |-- [ 19M]  ib_000004.SL2.50ch04.npz
            |-- [ 18M]  ib_000005.SL2.50ch05.npz
            |-- [ 11M]  ib_000006.SL2.50ch06.npz
            |-- [ 16M]  ib_000007.SL2.50ch07.npz
            |-- [ 16M]  ib_000008.SL2.50ch08.npz
            |-- [ 16M]  ib_000009.SL2.50ch09.npz
            |-- [ 17M]  ib_000010.SL2.50ch10.npz
            |-- [ 16M]  ib_000011.SL2.50ch11.npz
            `-- [ 20M]  ib_000012.SL2.50ch12.npz

 8.5G used in 4 directories, 29 files
```

## Requirements

```txt
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
```

## Sizes

```txt
1023M 360_merged_2.50.vcf.gz

w/o alignment
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
```

## Timing

```txt
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
```