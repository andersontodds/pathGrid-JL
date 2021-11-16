# getpaths.jl -- Julia function implementation of getPathsFromAP.m

using MAT
using Printf
using BenchmarkTools

function getpaths(startdate; stopdate=startdate, frames=144)

    # load stations.mat or stations.dat
    statfile = matopen("stations.mat")
    stations = read(statfile,"stations")
    close(statfile)
    # add line numbers to stationlist
    stnum = 1:1:size(stations,1)
    stations = cat(stationlist[:,1:3], stnum; dims=2)

    APdata = zeros(0,10)        # APfile.data has dimension N x 10
    APpower = zeros(70,0)       # APfile.power has dimension 70 x N

    for k = startdate:stopdate
        APfilename = @sprintf("C:/Users/tshelbya/Documents/MATLAB/Lightning_matlab/APfiles/AP%d.mat", k)
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

    return pairlist

end

## scratch

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
