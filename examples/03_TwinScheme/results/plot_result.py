from snplot import *

strain_rates = ["7E-5", "1E-1", "1E3"]
d = []
for rate in strain_rates:
    filename = f"str_{rate}.csv"
    file_data = read_csv(filename).loc[:, ["E33", "S33"]]
    ld = linedata(file_data["E33"],file_data["S33"],rate)
    d.append(ld)
p = xyplot(d, fig_name=f"stress", xlabel="Strain", ylabel="Stress (MPa)")
p.set_xlim((0,0.2))
p.set_ylim((0,150))
p.save()

d = []
for rate in strain_rates:
    filename = f"defect_{rate}.csv"
    file_data = read_csv(filename).loc[:, ["EVM", "Mode 1"]]
    ld = linedata(file_data["EVM"],file_data["Mode 1"],rate)
    d.append(ld)
p = xyplot(d, fig_name=f"dislocation", xlabel="Strain", ylabel="Dislocation Density (1/m2)")
p.set_yscale("log")
p.set_xlim((0,0.2))
p.set_ylim((1E10,1E16))
p.save()

d = []
for rate in strain_rates:
    filename = f"defect_{rate}.csv"
    file_data = read_csv(filename).loc[:, ["EVM", "Mode 2"]]
    ld = linedata(file_data["EVM"],file_data["Mode 2"],rate)
    d.append(ld)
p = xyplot(d, fig_name=f"twin", xlabel="Strain", ylabel="Twin Volume Fraction")
p.set_xlim((0,0.2))
p.set_ylim((0,0.5))
p.save()
