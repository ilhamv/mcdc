# Kobayashi dog-leg examples

- `steady_state`: original steady-state shielding benchmark from [Kobayashi et al. (2001)](https://doi.org/10.1016/S0149-1970%2801%2900007-5).
- `pulsed`: adaptation of the [time-dependent benchmark on Zenodo](https://doi.org/10.5281/zenodo.15069882), with a 0–50 s source pulse, an outlet detector, and tallies through 200 s.
- `pulsed_with_fission`: pulsed source with a fuel sphere at the second turn, an outlet detector, and tallies through 500 s, emphasizing the later fission contributions.

The fuel sphere is purely fissioning, with a fission cross section of 0.1 cm⁻¹, zero scattering and capture, and a prompt fission multiplicity of 2.5.

Each folder contains `input.py` and `process-output.py`.
Run from that folder:

```sh
python input.py
python process-output.py
```

Both pulsed examples use a cylindrical detector parallel to the outlet channel, centered at x=z=35 cm and spanning y=90–100 cm, with diameter 8 cm (four-fifths of the channel side).
The detector is filled with shield material (capture and scattering each 0.05 cm⁻¹), and the surrounding channel retains its low-density material.
Their processors show a flux map above density and detector response curves, with uncertainty bands and moving time markers.
Detector capture uses 1 s bins; flux and density use 10 s bins.
