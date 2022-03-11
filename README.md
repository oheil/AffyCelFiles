# AffyCelFiles

This is a low level IO package for Affymetrix .CEL files

**References**

https://media.affymetrix.com/support/developer/powertools/changelog/file-formats.html
https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html
https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/generic.html
https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html
https://media.affymetrix.com/support/developer/powertools/changelog/file-format-pgf.html
https://media.affymetrix.com/support/developer/powertools/changelog/file-format-clf.html
https://media.affymetrix.com/support/developer/powertools/changelog/file-format-mps.html

## Currently supported microarrays

This package is in an early stage and the following microarrays .CEL files are tested:
* Clariom_S_Human
* Mouse430_2
* HG-U133_Plus_2
* HuGene-2_0-st

(more to come)

The expression/intensity data was cross checked against R/bioconductor affy and oligo libraries.

## Dependencies

#### Julia versions

* Julia 1.0 or above

#### Third party packages

* none

#### Standard Library packages

* CRC32c

## Usage

In general Affymetrix .CEL can always be read in. To read Affymetrix .CEL files with meaningful information you need one of the following additional file/files:
* .cdf
* .pgf + .clf (+ optional .mps)

Those files provide the mapping from a probe location on the chip to a probe_id and the mapping from probeset_ids to probe_ids. For a specific chip type you
can get those files from the ThermoFisher support page, e.g. for the Clariom S human see https://www.thermofisher.com/order/catalog/product/902927?SID=srch-srp-902927
for download of archive TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1.zip (you need to have a login) which contains the .pgf,.clf and .mps file for this chip.

For biological meaningful analysis you also need annotation data, which maps probeset_ids to, for example, gene names. This is not part of this package, but
affymetrix annotation data is typically just a file of annotations in a table/csv format which can easily be read in (e.g. with CSV.jl) and mapped using probeset_ids or transcript_ids/meta_probeset_ids (in case of .mps used).

In this early stage this package doesn't provide a very convenient API, but this may change in future.

#### Usage examples

```julia
using Pkg
Pkg.add("AffyCelFiles");
#or from github main:
#Pkg.add(url="https://github.com/oheil/AffyCelFiles.jl",rev="main")

using AffyCelFiles
```

Reading .CEL files with corresponding .cdf file:
```julia
cel_file = "example_(Mouse430_2).CEL"
cdf_file = raw"d:\Temp\mouse430_2_libraryfile\CD_Mouse430_2\Full\Mouse430_2\LibFiles\Mouse430_2.cdf"

cel_data = AffyCelFiles.cel_read(cel_file);
cdf_data = AffyCelFiles.cdf_read(cdf_file);

intensities = AffyCelFiles.intensities(cel_data,cdf_data);

#Dict mapping probeset_id to lists of expression values (pm=perfect match, mm=mismatch).
#Depending on the chiptype mm can be empty
intensities.pm
intensities.mm
```
Result:
```julia
julia> intensities.pm
Dict{String, Vector{Float32}} with 45101 entries:
  "1453553_at"   => [93.0, 316.0, 74.0, 214.0, 94.0, 112.0, 52.0, 319.0, 66.0, 74.0, 88.0]
  "1455273_at"   => [1054.0, 57.0, 102.0, 263.0, 180.0, 447.0, 498.0, 683.0, 1233.0, 196.0, 999.0]
  "1425429_s_at" => [59.0, 162.0, 79.0, 286.0, 114.0, 384.0, 305.0, 111.0, 106.0, 138.0, 572.0]
  "1445844_at"   => [129.0, 749.0, 155.0, 156.0, 51.0, 144.0, 53.0, 123.0, 136.0, 102.0, 63.0]
  "1443594_at"   => [82.0, 150.0, 306.0, 447.0, 392.0, 101.0, 79.0, 143.0, 108.0, 531.0, 70.0]
  "1449219_at"   => [123.0, 104.0, 118.0, 74.0, 76.0, 136.0, 116.0, 165.0, 175.0, 116.0, 88.0]
  "1447933_at"   => [68.0, 120.0, 164.0, 42.0, 57.0, 50.0, 228.0, 272.0, 92.0, 68.0, 65.0]
...
julia> intensities.mm
Dict{String, Vector{Float32}} with 45101 entries:
  "1453553_at"   => [124.0, 75.0, 56.0, 317.0, 69.0, 165.0, 43.0, 148.0, 100.0, 65.0, 143.0]
  "1455273_at"   => [171.0, 51.0, 284.0, 93.0, 64.0, 174.0, 284.0, 233.0, 169.0, 48.0, 119.0]
  "1425429_s_at" => [51.0, 200.0, 80.0, 135.0, 94.0, 254.0, 255.0, 98.0, 93.0, 61.0, 161.0]
  "1445844_at"   => [79.0, 255.0, 105.0, 97.0, 65.0, 90.0, 44.0, 152.0, 187.0, 56.0, 40.0]
  "1443594_at"   => [57.0, 171.0, 181.0, 522.0, 244.0, 70.0, 70.0, 87.0, 80.0, 1093.0, 53.0]
...
```

Reading .CEL files with corresponding .pgf and .clf files (you need both):
```julia
cel_file = "example_(Clariom_S_Human).CEL"
pgf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.pgf"
clf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.clf"

cel_data = AffyCelFiles.cel_read(cel_file);
pgf_data = AffyCelFiles.pgf_read(pgf_file);
clf_data = AffyCelFiles.clf_read(clf_file);

intensities = AffyCelFiles.intensities(cel_data, pgf_data, clf_data);

#Dict mapping probeset_id to lists of expression values (pm=perfect match, mm=mismatch).
#Depending on the chiptype mm can be empty
intensities.pm
intensities.mm
```
Result:
```julia
julia> intensities.pm
Dict{String, Vector{Float32}} with 27189 entries:
  "23050629" => [61.0, 43.0, 54.0, 60.0]
  "23054888" => [321.0, 184.0, 169.0, 297.0, 162.0, 217.0, 194.0, 647.0, 124.0, 192.0]
  "23060812" => [56.0, 46.0, 145.0, 58.0, 92.0, 80.0, 204.0, 152.0, 83.0, 118.0]
  "23056335" => [160.0, 227.0, 951.0, 264.0, 286.0, 52.0, 43.0, 53.0, 48.0, 40.0]
  "23059041" => [1043.0, 867.0, 321.0, 739.0, 826.0, 462.0, 330.0, 1190.0, 1485.0, 1537.0]
...
julia> intensities.mm
Dict{String, Vector{Float32}} with 27189 entries:
  "23050629" => []
  "23054888" => []
  "23060812" => []
...
```
Providing a .mps file changes the central ids from probeset_ids to meta_probeset_ids/transcript_ids.
Reading .CEL files with corresponding .pgf, .clf files (you need both) and optional a .mps file:
```julia
cel_file = "example_(Clariom_S_Human).CEL"
pgf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.pgf"
clf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.clf"
mps_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.mps"

cel_data = AffyCelFiles.cel_read(cel_file);
pgf_data = AffyCelFiles.pgf_read(pgf_file);
clf_data = AffyCelFiles.clf_read(clf_file);
mps_data = AffyCelFiles.mps_read(mps_file);

intensities = AffyCelFiles.intensities(cel_data, pgf_data, clf_data, mps_data);

#Dict mapping meta_probeset_id/transcript_id to lists of expression values (pm=perfect match, mm=mismatch).
#Depending on the chiptype mm can be empty
intensities.pm
intensities.mm
```
Result:
```julia
julia> intensities.pm
Dict{String, Vector{Float32}} with 24351 entries:
  "TC0600009248.hg.1" => [42.0, 35.0, 195.0, 43.0, 140.0, 207.0, 54.0, 111.0, 95.0, 70.0]
  "TC0800011018.hg.1" => [47.0, 46.0, 41.0, 51.0, 74.0, 43.0, 62.0, 45.0, 42.0, 34.0]
  "TC0500012822.hg.1" => [299.0, 189.0, 425.0, 337.0, 1465.0, 1790.0, 988.0, 784.0, 471.0, 1094.0]
  "23050629"          => [61.0, 43.0, 54.0, 60.0]
  "TC1200012657.hg.1" => [2004.0, 1886.0, 450.0, 374.0, 175.0, 893.0, 332.0, 1710.0, 1930.0, 1362.0]
  "TC0700012299.hg.1" => [129.0, 77.0, 53.0, 76.0, 96.0, 395.0, 73.0, 130.0, 95.0, 82.0]
...
```









