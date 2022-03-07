# AffyCelFiles
read Affymetrix .CEL files

# Example
using AffyCelFiles

cdf_file = raw"d:\Temp\mouse430_2_libraryfile\CD_Mouse430_2\Full\Mouse430_2\LibFiles\Mouse430_2.cdf"
cdf_data = AffyCelFiles.cdf_read(cdf_file);

cel_file = "example_(Mouse430_2).CEL"
cel_data = AffyCelFiles.cel_read(cel_file);

intensities = AffyCelFiles.intensities(cel_data,cdf_data)

intensities.pm
intensities.mm



