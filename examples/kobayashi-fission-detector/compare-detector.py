"""Compare detector capture rates with and without fission (bands: ±1σ)."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

directory = Path(__file__).resolve().parent
fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")

for filename, label in [
    ("without_fission.h5", "Without fission"),
    ("with_fission.h5", "With fission"),
]:
    with h5py.File(directory / filename, "r") as f:
        tally = f["tallies/detector"]
        time = tally["grid/time"][:]
        # Normalize each file using its own detector time-bin widths.
        rate = tally["capture/mean"][:] / np.diff(time)
        sd = tally["capture/sdev"][:] / np.diff(time)

    mid = 0.5 * (time[:-1] + time[1:])
    (line,) = ax.plot(mid, rate, label=label)
    ax.fill_between(mid, rate - sd, rate + sd, color=line.get_color(), alpha=0.2)

ax.set_xlabel("Time [s]")
ax.set_ylabel("Detector capture rate per source neutron [s$^{-1}$]")
ax.set_yscale("log")
ax.set_title("Detector response (shading: ±1σ)")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.savefig(directory / "detector-comparison.png", dpi=200)
plt.show()
