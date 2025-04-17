from snplot.data import linedata
from snplot.plot import xyplot
from pandas import read_csv

case_folders = ["./1E-3", "./1E3", "./1E5", "./1E6", "./1E7", "./1E8"]
d = []
for case_folder in case_folders:
    data = read_csv(case_folder + "/str_str.csv")[["E33","S33"]]
    iline = linedata(data["E33"], data["S33"], case_folder)
    d.append(iline)
p = xyplot(d)
p.set_xlabel("Strain")
p.set_ylabel("Stress (MPa)")
p.set_xlim((0, 0.2))
p.show_tk()
