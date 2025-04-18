from snplot import *
from snplot.data import eulerdata
from snplot.plot import inverse_pole_figure_contour, pole_figure_contour

data_ = read_csv("../random_tex.txt", skiprows=4, sep="\s+")
data_.columns = [["phi1","PHI","phi2","weight"]]
E = eulerdata(data_["phi1"].values,data_["PHI"].values,data_["phi2"].values,"tex")
E.to_inverse_pole_figure("z")
p = inverse_pole_figure_contour([E],"origin")
p.rc_params['snplot.contour.lim']=(0,2.5)
p.have_colorbar = True
p.save()

strain_rates = ["7E-5","1E-1","1E3","5E7"]
for rate in strain_rates:
    data_ = read_csv(f"./Tex_{rate}.tex", skiprows=4, sep="\s+")
    data_.columns = [["phi1","PHI","phi2","weight"]]
    E = eulerdata(data_["phi1"].values,data_["PHI"].values,data_["phi2"].values,"tex")
    E.to_inverse_pole_figure("z")
    p = inverse_pole_figure_contour([E],rate+"_ipf")
    p.rc_params['snplot.contour.lim']=(0,4)
    p.have_colorbar = True
    p.save()
