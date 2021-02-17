begin

d = 2
p = 11
r = 100

f1 = addfile("data/hdmd_u_erainterim_d_" + d + "_p_" + p + "_r_" + r + ".nc", "r")
f2 = addfile("data/corr_hdmd_u_erainterim_d_" + d + "_p_" + p + "_r_" + r + ".nc", "r")

dmd_real = f1->dmd_real

dmd_real&lev@units = "hPa"

dmd_imag = f1->dmd_imag

dmd_imag&lev@units = "hPa"

lambda_real = f1->lambda_real
lambda_imag = f1->lambda_imag

corr = f2->corr

dmd_real_n = dmd_real({corr&ndmd}, :, :)
dmd_imag_n = dmd_imag({corr&ndmd}, :, :)
lambda_real_n = lambda_real({corr&ndmd})
lambda_imag_n = lambda_imag({corr&ndmd})

wks = gsn_open_wks("x11", "hdmd_u_erainterim_d_" + d + "_p_" + p + "_r_" + r)
plot = new(dimsizes(corr) * 2, graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnFillOn = True
res@cnFillPalette = "MPL_bwr"
res@cnLinesOn = True
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@lbLabelBarOn = False
res@gsnStringFontHeightF = 0.035

panel_string = new(dimsizes(plot), string)

do i = 0, (dimsizes(corr) - 1)
	res@gsnRightString = "~F33~l~F25~~B~R~N~ = " + sprintf("%7.4f", lambda_real_n(i))
	plot(i * 2) = gsn_csm_pres_hgt(wks, dmd_real_n(i, :, :), res)
	
	panel_string(i * 2) = "DMD " + corr&ndmd(i) + " Real"
	
	res@gsnRightString = "~F33~l~F25~~B~I~N~ = " + sprintf("%7.4f", lambda_imag_n(i))
	plot(i * 2 + 1) = gsn_csm_pres_hgt(wks, dmd_imag_n(i, :, :), res)
	
	panel_string(i * 2 + 1) = "DMD " + corr&ndmd(i) + " Imag"
end do

resP = True
resP@gsnMaximize = True
resP@gsnPanelMainString = "d = " + d + ", p = " + p + ", r = " + r
resP@gsnPanelFigureStrings = panel_string
resP@gsnPanelFigureStringsFontHeightF = 0.0125
resP@amJust = "TopLeft"
resP@gsnPanelLabelBar = True
resP@lbLabelFontHeightF = 0.015
resP@pmLabelBarWidthF = 0.6
resP@pmLabelBarHeightF = 0.05
resP@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/dimsizes(corr), 2/), resP)

end
