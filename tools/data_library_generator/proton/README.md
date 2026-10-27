# Proton data preparation

The proton generator converts TENDL ACE and optional NIST PSTAR data into the
MC/DC `mcdc-data` H1 layout.

For requirements, environment variables, command-line options, supported
reactions, and output details, see the
[Proton data preparation guide](../../../docs/source/user_guide/data_preparation/proton.rst).

Run the generator from the repository root:

```bash
python tools/data_library_generator/proton/generate.py
```
