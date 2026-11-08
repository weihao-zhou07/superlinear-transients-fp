# superlinear-transients-fp

Code and data for the paper

> **Superlinear transients of a 3x+1-type polynomial map over F_p**

The manuscript is submitted separately. All mathematical proofs are by hand.
This repository contains the programs and outputs used for

- the computational consistency checks of hand proofs,
- the observations for other primes in Section 8,
- the evidence for the conjectures.

No proof depends on these computations.

## For referees: quick start

```bash
pip install -r requirements.txt
cd data
python ../scripts/gen_census_data_v6_4.py
```

A successful run ends with

```
triples: 1652
all checks passed
```

To run the verification suites (optional, about ten minutes):

```bash
cd ../logs
python ../scripts/verify_v5.py
python ../scripts/verify_v5_supp.py
```

## Where the data in the paper are

| Paper item | File |
|---|---|
| Table 4: smallest rigid levels, rho, C_rho | `data/table_v7_observations.csv` |
| p=11, W=822: complete Lambda and Z_c data | `data/finite_data_p11_W822.csv` |
| p=23, W=158: complete Lambda and Z_c data | `data/finite_data_p23_W158.csv` |
| Generator output: all checks passed | `logs/gen_census_data.log` |
| Standard verification run (A-J) | `logs/verify_v5_log_std.txt` |
| Extended verification run (E,G) | `logs/verify_v5_log.txt` |
| Supplementary verification run (S0-S3) | `logs/verify_v5_supp_log.txt` |

The generator `scripts/gen_census_data_v6_4.py` reproduces the first two CSVs
and checks Tables 3-5, the p=5 hand proof data, the `m0(99573)` decomposition
and the finite-data appendix.

## Repository layout

```
scripts/
  gen_census_data_v6_4.py   generates the observation tables and CSV data;
                            checks Tables 3-5 and the finite data
  verify_v5.py              A-J verification suite (exact integer/numpy)
  verify_v5_supp.py         S0-S3 supplementary checks
  m0_regularity.py          experiments on p-regularity of m0
  m0_probe*.py, m0_traj.py, m0_dichotomy.py
                            exploratory drivers for m0_regularity
  hankel2.py                Hankel-rank test for p-regularity
data/
  table_v7_observations.csv the ten rows of Table 4 in the v7 manuscript
  finite_data_p11_W822.csv  complete Lambda/Z_c data, p=11, W=822
  finite_data_p23_W158.csv  complete Lambda/Z_c data, p=23, W=158
logs/
  gen_census_data.log       output of the generator (all checks passed)
  verify_v5_log_std.txt     standard verification run (A-J)
  verify_v5_log.txt         extended verification run (E,G)
  verify_v5_supp_log.txt    supplementary verification run (S0-S3)
```

Additional generated files in `data/` (`tab_superlinear.tex`,
`decomp_99573.tex`, `count_triples.tex`, `appendix_rigid.tex`) and
`logs/static_check.txt` come from the earlier v6.4 LaTeX layout. They are kept
for provenance and are **not needed** to read the v7 paper or its data.

The verification scripts keep their original names (`verify_v5*.py`) and logs;
only the v7 manuscript and its data map above are current.

## Requirements

- Python 3.8 or newer
- `numpy` (used by `scripts/verify_v5.py`)

The checked-in logs were produced with Python 3.13.14 and numpy 2.5.4. The
standard run took about 550 s, the extended E,G run about 29 s, and the
supplementary run about 44 s.

## Exploratory scripts

`m0_regularity.py`, `m0_probe*.py`, `m0_traj.py`, `m0_dichotomy.py` and
`hankel2.py` are exploratory and are not needed for the paper's claims. Their
usage is documented in their file headers, for example

```bash
python scripts/hankel2.py 5 9 30
```

## License

MIT; see the `LICENSE` file.
