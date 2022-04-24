import pandas as pd
import numpy
from matplotlib import pyplot as plt
from matplotlib import style

fig = plt.figure()
style.use("ggplot")
plot_axis = fig.add_subplot(111)
plot_axis.set(title="experimental data", ylabel="indicator", xlabel="temperature")

experimental_data = pd.read_csv("../src_data/test_data.csv")

plot_axis.scatter(experimental_data.iloc[:, [1]], experimental_data.iloc[:, [0]], color="blue",
                  label="experimental_data",
                  marker='*')

plot_axis.scatter(-95.3644 * numpy.log10(experimental_data.iloc[:, [0]].to_numpy()) + 207.492,
               experimental_data.iloc[:, [0]], color="red", label="processed")

plot_axis.legend()
plt.show()
