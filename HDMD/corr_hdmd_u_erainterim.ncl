begin

d = 2
p = 11
r = 100

print("Reading the input file...")

f1 = addfile("data/eof.u.erainterim.nc", "r")
f2 = addfile("data/hdmd_u_erainterim_d_" + d + "_p_" + p + "_r_" + r + ".nc", "r")

eof = f1->eof

dmd_real = dble2flt(f2->dmd_real)
dmd_imag = dble2flt(f2->dmd_imag)

print("Calculating correlation...")

eof_1d = new(2 * dimsizes(eof&lat) * dimsizes(eof&lev), float, -999)
eof_1d(0:(dimsizes(eof&lat) * dimsizes(eof&lev) - 1)) = ndtooned(eof(0, :, :))
eof_1d((dimsizes(eof&lat) * dimsizes(eof&lev)):(2 * dimsizes(eof&lat) * dimsizes(eof&lev) - 1)) = ndtooned(eof(1, :, :))

eof_1d = abs(eof_1d)

dmd_1d = new(2 * dimsizes(eof&lat) * dimsizes(eof&lev), float, -999)

corr = new(dimsizes(dmd_real&ndmd), float)
corr!0 = "ndmd"
corr&ndmd = dmd_real&ndmd

do i = 0, (dimsizes(corr) - 1)
	dmd_real_i = dmd_real(i, :, :)
	dmd_imag_i = dmd_imag(i, :, :)
	
	dmd_1d(0:(dimsizes(eof&lat) * dimsizes(eof&lev) - 1)) = ndtooned(dmd_real_i)
	dmd_1d((dimsizes(eof&lat) * dimsizes(eof&lev)):(2 * dimsizes(eof&lat) * dimsizes(eof&lev) - 1)) = ndtooned(dmd_imag_i)
	
	dmd_1d = abs(dmd_1d)
	
	corr(i) = (/escorc(eof_1d, dmd_1d)/)
end do

corr_sort = abs(corr)
copy_VarCoords(corr, corr_sort)

qsort(corr_sort)

corr_sort_unique = get_unique_values(corr_sort)

n = 3

corr_sort_unique_n = corr_sort_unique((dimsizes(corr_sort_unique) - 1):(dimsizes(corr_sort_unique) - 1 - (n - 1)))

index_corr_sort_unique = new(n, integer)

do i = 0, (n - 1)
	index = ind(abs(corr).eq.corr_sort_unique_n(i))
	if (dimsizes(index).gt.1) then
		index_corr_sort_unique(i) = index(0)
	else
		index_corr_sort_unique(i) = index
	end if
	delete(index)
end do

corr_sort_unique_n!0 = "ndmd"
corr_sort_unique_n&ndmd = index_corr_sort_unique + 1

print("DMD " + corr_sort_unique_n&ndmd + ": r = " + corr_sort_unique_n)

print("Opening the output file...")

system("rm data/corr_hdmd_u_erainterim_d_" + d + "_p_" + p + "_r_" + r + ".nc")
g = addfile("data/corr_hdmd_u_erainterim_d_" + d + "_p_" + p + "_r_" + r + ".nc", "c")

print("Writing the output file...")

g->corr = corr_sort_unique_n

print("Closing the output file...")

delete(g)

end
