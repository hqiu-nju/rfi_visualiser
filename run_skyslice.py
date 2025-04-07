#### this is a script to run oskar on the sky models and generate slices for the every interval
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
import oskar


settings = oskar.SettingsTree("oskar_sim_interferometer")
settings.from_dict(params)