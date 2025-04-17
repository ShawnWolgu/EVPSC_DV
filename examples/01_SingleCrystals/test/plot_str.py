from snplot.data import linedata
from snplot.plot import xyplot
from pandas import read_csv

ds = []
for icase in ["001", "111", "112", "213"]:
    df = read_csv(f"./str_{icase}.csv")
    d = linedata(df["E33"], df["S33"], label=icase)
    ds.append(d)
p = xyplot(ds, xlabel="strain", ylabel="stress (MPa)", title="Stress-strain curve")
p.set_fig_name("str_5E7")
p.expand_args()
p.show_tk()
