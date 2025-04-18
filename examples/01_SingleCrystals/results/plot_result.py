from snplot import *

axes = ["001", "111", "112", "213"]
strain_rates = ["7E-5", "1E-1", "1E3", "5E7"]
for rate in strain_rates:
    d = []
    for axis in axes:
        filename = f"str_{axis}_{rate}.csv"
        file_data = read_csv(filename).loc[:, ["E33", "S33"]]
        ld = linedata(file_data["E33"],file_data["S33"],axis)
        d.append(ld)
    p = xyplot(d, fig_name=f"strain rate_{rate}", x_label="Strain", y_label="Stress (MPa)")
    p.save()

