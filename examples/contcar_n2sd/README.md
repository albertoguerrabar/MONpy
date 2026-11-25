# CONTCAR_n2SD dimer

This case has 174 atoms split into two fragments of 87 atoms each. The ATOM
numbers of the two molecules are interleaved. The TOML records the complete
partition. Fragment A contains atom 1.

Run the bundled signed-CI NDOL output from the repository root:

```bash
python main.py --config examples/contcar_n2sd/monpy.toml
```

The NDOL file is about 51 MB. The example uses `color_gamma = 0.5` to make
small fragment contributions easier to see.

| State | Energy (eV) | CT fraction | Character |
| --- | ---: | ---: | --- |
| S1 | 2.6224 | 0.000832 | Local on B |
| S9 | 3.9278 | 0.766304 | B → A |
| T15 | 3.9415 | 0.646234 | B → A |

[Reference heatmap for S9](../../docs/source/_static/monpy_dimer_s9.png)
