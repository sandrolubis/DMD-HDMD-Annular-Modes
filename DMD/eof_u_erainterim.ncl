begin

print("Reading the input file...")

f = addfile("data/u.anom.erainterim.nc", "r")

u = f->u(time|:, {lev|100:1000}, {lat|-20:-80}, lon|0)

u@_FillValue = -999

print("Weighting latitude...")

rad = 4.0 * atan(1.0) / 180.0
wgt = sqrt(cos(u&lat * rad))

u = (/u * conform(u, wgt, 2)/)

print("Calculating EOFs...")

neof = 2

eof = eofunc_n_Wrap(u, neof, False, 0)
pc = eofunc_ts_n_Wrap(u, eof, False, 0)

printVarSummary(eof)
printVarSummary(pc)

print("Opening the output file...")

system("rm data/eof.u.erainterim.nc")
g1 = addfile("data/eof.u.erainterim.nc", "c")

system("rm data/pc.u.erainterim.nc")
g2 = addfile("data/pc.u.erainterim.nc", "c")

print("Writing the output file...")

g1->eof = eof
g2->pc = pc

print("Closing the output file...")

delete(g1)
delete(g2)

end
