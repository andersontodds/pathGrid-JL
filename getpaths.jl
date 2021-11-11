# getpaths.jl -- Julia function implementation of getPathsFromAP.m

using MAT
using Printf

function getpaths(startdate; stopdate=startdate, frames=144)

    # load stations.mat or stations.dat
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

    end

    return APpower
    #strows = APpower[1:2:end,:];

end

power = getpaths(20170906)
