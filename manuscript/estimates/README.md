# Order-of-magnitude estimates quoted in Sec. 4.1

Scaling estimates from stated case inputs, not CFD results. Re-run with
`python3 manuscript/estimates/<script>.py`; the `.txt` files are the outputs
quoted in the text.

- `dimensionless.py`: Re, radial diffusion time, wall Damkoehler number (Eq. damkohler).
- `buoyancy_soret.py`: wall-to-bulk temperature lag, axial-gradient Grashof and
  Richardson numbers, compositional buoyancy, Soret drift bounds (alpha_T = 1).
