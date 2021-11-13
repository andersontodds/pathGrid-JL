# getpaths.jl -- Julia function implementation of getPathsFromAP.m

using MAT
using Printf

function getpaths(startdate; stopdate=startdate, frames=144)

    # load stations.mat or stations.dat
    statfile = matopen("stations.mat")
    stations = read(statfile,"stations")
    close(statfile)

    APdata = zeros(0,10)        # APfile.data has dimension N x 10
    APpower = zeros(70,0)       # APfile.power has dimension 70 x N

    for k = startdate:stopdate
        APfilename = @sprintf("C:/Users/tshelbya/Documents/MATLAB/Lightning_matlab/APfiles/AP%d.mat", k)
        APfile = matopen(APfilename)
        APfiledata = read(APfile,"data")
        APfilepower = read(APfile,"power")

        # concatenate total data, power with current file data, power
        APdata = [APdata; APfiledata]
        APpower = [APpower APfilepower]
        close(APfile)
    end

    return APpower

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
    stmat[2:2:end, :] = nullrows;   # is it necessary to pad stmat with NaNs in order to make its size equal to APpower?

    # find every stationID; as of 4 Oct 2018 there are 122 entries in stations.dat
    # NOTE: stations.dat and .mat index from stID == 1 (Dunedin), not 0
    pairlist = zeros(0,13)    # "strokelist", in MATLAB version
    for stID = 1:122
        #NOTE: MATLAB code below.  Need to find the columns in APpower in which stID exists, then find the APdata corresponding to those detections, and thereby compile subset of APdata that includes only strokes detected by current stID.

        stindex = findall(x -> x==stID, stmat)
        #lidx = sub2ind(size(APpower),row + 1, col)
        #APdata_stID = flrAPdata(col,:);
        #APpower_stID = flrAPpower(lidx);

        # return list of stroke time, location for each station
        # need dimension stID, or vectorize and concatenate

        #strokecount_stID = size(APdata_stID,1);
        #c = ones(strokecount_stID,1);

        #st_lat = stations{stID,1};
        #st_lon = stations{stID,2};

        #strokelist_stID = cat(2,APdata_stID,stID*c,st_lat*c,st_lon*c);

        #strokelist = cat(1,strokelist,strokelist_stID);
    end


end

## scratch

power = getpaths(20170906)

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



stID = 11
stindex = findall(x -> x==stID, stmat) # how to separate rows and columns in vector of CartesianIndex{2}?
