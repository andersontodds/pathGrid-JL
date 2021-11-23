# getpaths.jl -- Julia function implementation of getPathsFromAP.m

using MAT
using Printf
using BenchmarkTools
using Profile
using Dates
using CompoundPeriods

function getpaths(startdate; stopdate=startdate)

    # load stations.mat or stations.dat
    statfile = matopen("stations.mat")
    stations = read(statfile,"stations")
    close(statfile)
    # add line numbers to stationlist
    stnum = 1:1:size(stations,1)
    stations = cat(stations[:,1:3], stnum; dims=2)

    APdata = zeros(0,10)        # APfile.data has dimension N x 10
    APpower = zeros(70,0)       # APfile.power has dimension 70 x N

    datestart = DateTime(string(startdate), "yyyymmdd")
    datestop = DateTime(string(stopdate), "yyyymmdd")

    for k = datestart:Day(1):datestop      # possibility of crossing month boundary in range requires indexing through Date rather than Int, then converting back to something parseable by @sprintf
        k_str = Dates.format(k, "yyyymmdd")
        APfilename = @sprintf("C:/Users/tshelbya/Documents/MATLAB/Lightning_matlab/APfiles/AP%s.mat", k_str)
        APfile = matopen(APfilename)
        APfiledata = read(APfile,"data")
        APfilepower = read(APfile,"power")

        # concatenate total data, power with current file data, power
        APdata = cat(APdata, APfiledata; dims=1)
        APpower = cat(APpower, APfilepower; dims=2)
        close(APfile)
    end

    # return APdata

    # convert station ID from 0-start (in stations.mat, stations.dat) to 1-start; except station 0 (Dunedin)
    # since stations are listed in order of ascending stID in each column of APpower[], can add 1 to all station IDs, then find all station IDs not on the first row and set them back to 0.  The fastest alternative that I know of would be to add 1 to all elements of the first row of strows[], then add 1 to all nonzero rows of @view strows[2:end, :].
    strows = APpower[1:2:end,:]
    strows .+= 1.0
    notfirstrow = @view strows[2:end,:]
    falseones = findall(x -> x == 1.0, notfirstrow)
    notfirstrow[falseones] .= 0.0

    nullrows = NaN.*ones(size(strows))
    stmat = ones(size(APpower))
    stmat[1:2:end, :] = strows;
    stmat[2:2:end, :] = nullrows;   #QUESTION is it necessary to pad stmat with NaNs in order to make its size equal to APpower? → useful if station energy needs to be used later.

    # find every stationID; as of 4 Oct 2018 there are 122 entries in stations.dat
    # NOTE: stations.dat and .mat index from stID == 1 (Dunedin), not 0
    pairlist = zeros(0,13)    # "strokelist", in MATLAB version
    for stID = 1:122
        stindex = Tuple.(findall(x -> x==stID, stmat)) # defaults to CartesianIndex without Tuple.()
        foundrows = first.(stindex)
        foundcols = last.(stindex)

        APdata_stID = APdata[foundcols,:]

        st_lat = stations[stID, 1]
        st_lon = stations[stID, 2]

        strokecount_stID = size(foundcols, 1)
        strokelist_stID = cat(APdata_stID, stID*ones(strokecount_stID,1), st_lat*ones(strokecount_stID,1), st_lon*ones(strokecount_stID,1); dims=2)

        pairlist = cat(pairlist, strokelist_stID; dims=1)

    end

    #return pairlist

    # TODO 11/16/2021: get pairlist_lite_10m and file output working.

    # collapse Nx6 time to Nx1 DateTime
    stroke_time = (pairlist[:,1:6])
    tcp = Year.(stroke_time[:,1]) + Month.(stroke_time[:,2]) + Day.(stroke_time[:,3]) +
        + Hour.(stroke_time[:,4]) + Minute.(stroke_time[:,5]) +
        + Nanosecond.(round.(Int,1e9*stroke_time[:,6]))
    t = canonical.(tcp)
    stroke_dt = DateTime.(t)
    # could preserve microsecond accuracy by either keeping only Nx1 CompoundPeriods, or Nx2 Date,Time
    #t_date = Date.(t)
    #t_time = Time.(t)

    stroke_lat = pairlist[:,7]
    stroke_lon = pairlist[:,8]
    stat_lat = pairlist[:,12]
    stat_lon = pairlist[:,13]
    pairlist_lite = cat(stroke_dt, stroke_lat, stroke_lon, stat_lat, stat_lon, dims=2)


    bin_edges_10min = datestart:Minute(10):datestop

    for bin in bin_edges_10min[1:end-1]

        index_10m = findall(x-> (x >= bin) & (x <= bin+Minute(10)), timecol)
        pairlist_lite_10m = pairlist_lite[index_10m, :]
        save_df = Dates.format(bin, "yyyymmddHHMM")
        savename = @sprintf("pairlist_lite_10m_%s.jld2", save_df)
        jldsave(savename, pairlist=pairlist_lite_10m)
        mv(savename, @sprintf("test/%s",savename)) # move data file to data/ directory

    end

    return pairlist_lite

end

## scratch

# test function
# NOTE 11/17: need to profile smaller parts of function.  In current form, getpaths(20170906) does not finish in ~10 minutes.
# NOTE 11/18: JLD2 files are much larger than MAT files; e.g.:
#   pairlist_lite_10m_201709060000.jld2       9,679 KB
#   strokelist_lite_10m_201709060000.mat        376 KB
# JLD2 is in HDF5 format and supports "pretty much arbitrary Julia objects".  Perhaps try csv? Also try profiling saving different file types
# Also: when searching for times (e.g. building 10 minute pairlists), instead of doing findall(times in 10-minute window), first sort by time, then iterate down array until time window is exceeded.  This should reduce number of times array is traversed from N to 2.


Juno.@profiler getpaths(20170906)

## profiling
# building pairlist
@benchmark begin
    pairlist = zeros(0,13);
    for stID = 1:122
        stindex = Tuple.(findall(x -> x==stID, stmat)) # defaults to CartesianIndex without Tuple.()
        foundrows = first.(stindex)
        foundcols = last.(stindex)

        APdata_stID = APdata[foundcols,:]

        st_lat = stations[stID, 1]
        st_lon = stations[stID, 2]

        strokecount_stID = size(foundcols, 1)
        strokelist_stID = cat(APdata_stID, stID*ones(strokecount_stID,1), st_lat*ones(strokecount_stID,1), st_lon*ones(strokecount_stID,1); dims=2)

        pairlist = cat(pairlist, strokelist_stID; dims=1)

    end
end
# @benchmark result:
  # memory estimate:  29.02 GiB
  # allocs estimate:  16253
  # --------------
  # minimum time:     33.025 s (10.85% GC)
  # median time:      33.025 s (10.85% GC)
  # mean time:        33.025 s (10.85% GC)
  # maximum time:     33.025 s (10.85% GC)
  # --------------
  # samples:          1
  # evals/sample:     1

# saving various file formats
# JLD2:
  # memory estimate:  33.93 MiB
  # allocs estimate:  991500
  # --------------
  # minimum time:     266.250 ms (0.00% GC)
  # median time:      324.783 ms (6.04% GC)
  # mean time:        325.528 ms (4.98% GC)
  # maximum time:     363.220 ms (4.88% GC)
  # --------------
  # samples:          16
  # evals/sample:     1
testarray = rand(Float64, (30000,5))
savename = "testfloats.jld2"
@benchmark jldsave(savename, ta=testarray)
#  memory estimate:  9.05 KiB
#  allocs estimate:  103
#  --------------
#  minimum time:     6.583 ms (0.00% GC)
#  median time:      8.965 ms (0.00% GC)
#  mean time:        34.106 ms (0.00% GC)
#  maximum time:     882.296 ms (0.00% GC)
#  --------------
#  samples:          160
#  evals/sample:     1
savename = "testfloats.csv"
@benchmark @save savename testarray

##
using Plots

APdata = getpaths(20170906)
APpower = power


strows = power[1:2:end,:]
strows = strows .+ 1

notfirstrow = @view strows[2:end,:]
falseones = findall(x -> x == 1.0, notfirstrow)
notfirstrow[falseones] .= 0.0


nullrows = NaN.*ones(size(strows))

stmat = ones(size(power))
stmat[1:2:end, :] = strows;
stmat[2:2:end, :] = nullrows;

stmat
#strows
#trueones = findall(x -> x == 1.0, strows)

statfile = matopen("stations.mat")
stations = read(statfile,"stations")
close(statfile)

#stnum = 1:1:size(stations,1)
#stations = cat(stnum, stationlist[:,1:2]; dims=2)

stID = 11
stindex = Tuple.(findall(x -> x==stID, stmat)) # how to separate rows and columns in vector of CartesianIndex{2}?
foundrows = first.(stindex)
foundcols = last.(stindex)

APdata_stID = APdata[foundcols,:]

st_lat = stations[stID, 1]
st_lon = stations[stID, 2]

strokecount_stID = size(foundcols, 1)
strokelist_stID = cat(APdata, st_lat*ones(strokecount_stID,1), st_lon*ones(strokecount_stID,1); dims=2)

strokelist = cat(strokelist, strokelist_stID; dims=1)

# benchmark for loop
pairlist = zeros(0,13)    # "strokelist", in MATLAB version
for stID = 1:122
    stindex = Tuple.(findall(x -> x==stID, stmat)) # defaults to CartesianIndex without Tuple.()
    foundrows = first.(stindex)
    foundcols = last.(stindex)

    APdata_stID = APdata[foundcols,:]

    st_lat = stations[stID, 1]
    st_lon = stations[stID, 2]

    strokecount_stID = size(foundcols, 1)
    strokelist_stID = cat(APdata_stID, stID*ones(strokecount_stID,1), st_lat*ones(strokecount_stID,1), st_lon*ones(strokecount_stID,1); dims=2)

    pairlist = cat(pairlist, strokelist_stID; dims=1)
end

pl_20170906 = getpaths(20170906)

k = 20170906
APfilename = @sprintf("C:/Users/tshelbya/Documents/MATLAB/Lightning_matlab/APfiles/AP%d.mat", k)
APfile = matopen(APfilename)
APfiledata = read(APfile,"data")

in_len = size(APfiledata,1)
out_len = size(pl_20170906,1)

out_len/in_len

# DateTime only supports millisecond precision.  Can either
#   a) round stroke_time to nearest millisecond,
#   b) split stroke_time into Date (YYYYMMDD) and Time (hh:mm:ss.)
stroke_time = (pairlist[:,1:6])
tcp = Year.(stroke_time[:,1]) + Month.(stroke_time[:,2]) + Day.(stroke_time[:,3]) +
    + Hour.(stroke_time[:,4]) + Minute.(stroke_time[:,5]) +
    + Nanosecond.(round.(Int,1e9*stroke_time[:,6]))

t = canonical.(tcp)

stroke_dt = DateTime.(t)
t_date = Date.(t)
t_time = Time.(t)

stroke_lat = pairlist[:,7]
stroke_lon = pairlist[:,8]
stat_lat = pairlist[:,12]
stat_lon = pairlist[:,13]
pairlist_lite = cat(stroke_dt, stroke_lat, stroke_lon, stat_lat, stat_lon, dims=2)

times = pairlist

size(times, 1)

times = [   2017 09 06 00 00 04.911360;
            2017 09 06 04 15 55.193727;
            2017 09 06 22 55 12.256655]

sdate = Date(times[1,1]):Day(1):Date(times[3,1])
stime = Time(times[1,4]):Second(1):Time(times[3,4])

sdate = Vector{Date}

for i = 1:size(times, 1)
    #stroke_date = times[1,1:3]
    years   = times[i,1]
    months  = times[i,2]
    days    = times[i,3]
    hours   = times[i,4]
    minutes = times[i,5]
    seconds = floor(times[i,6])
    millis  = floor((times[i,6] - seconds)*1000)
    micros  = floor(((times[i,6] - seconds)*1000 - millis)*1000)

    sdate[i] = Date(years, months, days)
    stime[i] = Time(hours, minutes, seconds, millis, micros)

end

i = 1
years   = times[i,1]
months  = times[i,2]
days    = times[i,3]
hours   = times[i,4]
minutes = times[i,5]
seconds = floor(times[i,6])
millis  = floor((times[i,6] - seconds)*1000)
micros  = floor(((times[i,6] - seconds)*1000 - millis)*1000)

sdate = Date(years, months, days)
stime = Time(hours, minutes, seconds, millis, micros)

date_vector = [[2017 09 06], [2017 09 07], [2017 09 08]]

sdate = Date.(date_vector)

fracseconds = times[1, 6]

datestr = @sprintf("%f %f %f %f %f %.3f", years, months, days, hours, minutes, fracseconds)

## @bernhard solution

using Dates
dfmt = DateFormat("yyyy mm dd HH MM SS.s")
DateTime("2017 09 06 00 00 04.911",dfmt) #works for me
DateTime("2017 09 06 00 00 04.911360",dfmt) #fails

#How about this:
using Dates
times=["2017 09 06 00 00 04.911360","2017 09 06 04 15 55.193727"]
dfmt = DateFormat("yyyy mm dd HH MM SS")
times_wo_mus = map(x->x[1:findfirst(".",x)[1]-1],times)
times_mus = map(x->x[findfirst(".",x)[1]+1:end],times)

#Parse DateTime (second precision only)
dt = DateTime.(times_wo_mus,dfmt)
sdate = Date.(dt)
stime = Time.(dt)

#add Microseconds
mus = Microsecond.(parse.(Int,times_mus))
stime .= stime .+ mus

## @rafael.guerra solution

using Dates, CompoundPeriods

T = [   2017 09 06 00 00 04.911360;
        2017 09 06 04 15 55.193727;
        2017 09 06 22 55 12.256655]

tcp = Year.(T[:,1]) + Month.(T[:,2]) + Day.(T[:,3]) + Hour.(T[:,4]) +
      + Minute.(T[:,5]) + Nanosecond.(round.(Int,1e9*T[:,6]))

t = canonical.(tcp)

t_dt = DateTime.(t)
t_date = Date.(t)
t_time = Time.(t)

## combine/solve

times = [   2017 09 06 00 00 04.911360;
            2017 09 06 04 15 55.193727;
            2017 09 06 22 55 12.256655]

years   = times[:,1]
months  = times[:,2]
days    = times[:,3]
hours   = times[:,4]
minutes = times[:,5]
seconds = times[:,6]

dfmt = DateFormat("yyyy mm dd HH MM SS.s")

timestring = fill("",size(times, 1),1)

for i = 1:size(times, 1)

    timestring[i] = @sprintf("%.0f %02.0f %02.0f %02.0f %02.0f %06.3f", years[i], months[i], days[i], hours[i], minutes[i], seconds[i])

end

dt = DateTime.(timestring,dfmt)

## bin into 10-minute files

dt_start = DateTime(2017,09,06,00,00,00)
dt_stop = DateTime(2017,09,07,00,00,00)
bin_edges_10min = dt_start:Minute(10):dt_stop

bin_edges_10min[2]

bin = bin_edges_10min[2]-Minute(10)
timecol = @view pairlist_lite[:,1]
size(pairlist_lite)
index_10m = findall(x-> (x >= bin) & (x < bin+Minute(10)), timecol)
timecol[index_10m[end]]
pairlist_10m_test = pairlist_lite[index_10m, :]

using Plots
plotly()
plot(1:1:length(index_10m),index_10m)

A = [1 2 3 4;
     5 6 7 8;
     9 10 11 12;]

firstcol = @view A[:,1]
fiveind = findall(x->x==5, firstcol)
A_rowwithfive = A[fiveind,:]

for bin in bin_edges_10min[1:end-1]

    index_10m = findall(x-> (x >= bin) & (x <= bin+Minute(10)), timecol)
    pairlist_lite_10m = pairlist_lite[index_10m, :]
    save_df = Dates.format(bin, "yyyymmddHHMM")
    savename = @sprintf("pairlist_lite_10m_%s.jld2", save_df)
    jldsave(savename, pairlist=pairlist_lite_10m)
    mv(savename, @sprintf("data/%s",savename)) # move data file to data/ directory

end

using FileIO
using JLD2



index_10m = findall(x-> (x >= bin) & (x <= bin+Minute(10)), timecol)
pairlist_lite_10m = pairlist_lite[index_10m, :]
save_df = Dates.format(bin,"yyyymmddHHMM")
savename = @sprintf("pairlist_lite_10m_%s.jld2",save_df)
jldsave(savename,pairlist=pairlist_lite_10m)
mv(savename, @sprintf("data/%s",savename)) # move data file to data/ directory


pairlist_lite_10m = pairlist_lite[index_10m, :]
save_df = Dates.format(bin,"yyyymmddHHMM")
savename = @sprintf("pairlist_lite_10m_%s.jld2",save_df)
jldsave(savename,pairlist=pairlist_lite_10m)
mv(savename, @sprintf("data/%s",savename))

save("pairlist_lite_10m_test.jld2", "pairlist_lite_10m_test", pairlist_lite[index_10m, :])
jldsave("example_pairlist_lite_10m.jld2", pl_10m=pairlist_lite[index_10m, :])
pl = jldopen("example_pairlist_lite_10m.jld2")
read(pl, "pl_10m")
