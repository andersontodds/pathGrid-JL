# pathGrid-JL

Julia implementation of the pathGrid MATLAB project.  Until this repository is released, see the pathGrid repository: https://github.com/andersontodds/pathGrid

# TODO
11-22-2021: CompoundDates.canonical() is a convenient way to store date series with high resolution/dynamic range, but extremely slow; profile reveals function canonical() represents the bulk of the getpaths() runtime.  Consider replacing with datenum() or similar.
11-17-2021: All components of getpaths() appear to work, but getpaths(20170906) does not finish within 10 minutes. See benchmark results.
11-16-2021: get time binning and fileoutput working with getpaths.jl (11/17: **DONE**); fill out outline.

# Functions
getpaths(date; fileout="lite")
  getpaths takes a date as input (integer of format YYYYMMDD, e.g. 20170906) and outputs a file or set of files containing the times and endpoints of lightning stroke-to-station paths.
  arguments:
    date    integer


# Outline
main
  pathfile():
    input: date range
    parameters: fileout
    calls: getpaths()
    getpaths()
      input: date
      output: stroke-station pair list
    output: file(s) of stroke-station pairs

or

  getpaths():
    input: date range
    parameter: fileout
    output: file(s) of stroke-station pairs
