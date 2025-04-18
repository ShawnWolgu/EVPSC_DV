from snplot import *

strain_rates = ["7E-5", "1E-1", "1E3", "5E7"]
d = []
for rate in strain_rates:
    filename = f"str_{rate}.csv"
    file_data = read_csv(filename).loc[:, ["E33", "S33"]]
    ld = linedata(file_data["E33"],file_data["S33"],rate)
    d.append(ld)
p = xyplot(d, fig_name=f"stress", x_label="Strain", y_label="Stress (MPa)")
p.set_yscale("log")
p.set_xlim((0,0.2))
p.set_ylim((0.1,10000))
p.save()

