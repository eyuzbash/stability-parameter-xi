# Stability Parameter ξ Calculator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15550102.svg)](https://doi.org/10.5281/zenodo.15550102)

This repository provides Mathematica, Python, and MATLAB codes for computing the stability parameter **ξ** from the Eliashberg function α²F(ω). The stability condition ξ < ξ_c ≤ 1 sets fundamental upper limits on the electron-phonon interaction strength and superconducting transition temperature in metals.

According to current experimental data, ξ_c ≈ 1/2 in hydrides and is significantly lower in conventional metals. These results are based on the theoretical analysis in:

> D. Semenok, B. L. Altshuler, and E. A. Yuzbashyan,  
> [arXiv:2407.12922](https://arxiv.org/abs/2407.12922)

---

## 📁 Repository Structure

```bash
.
├── README.md               # This file
├── LICENSE                 # MIT License
├── .gitignore              # Common excludes for MATLAB, Python, Mathematica
├── data/                   # Sample α²F(ω) input files
├── python/                 # Python script implementation
├── matlab/                 # MATLAB script implementation
├── mathematica/            # Mathematica notebook
└── docs/                   # Optional documentation or figures
```

---

## ⚙️ Requirements

### Python

- Python ≥ 3.8
- Libraries: `numpy`, `scipy`, `mpmath`, `matplotlib`

Install via:

```bash
pip install numpy scipy mpmath matplotlib
```

Then run:

```bash
python python/compute_xi.py
```

### MATLAB

Run the script `matlab/compute_xi.m` from within MATLAB. No special toolboxes are required.

### Mathematica

Open and evaluate `mathematica/compute_xi.nb` cell-by-cell in Wolfram Mathematica.

---

## 📝 Citation

If you use this code, please cite:

> D. Semenok, B. L. Altshuler, and E. A. Yuzbashyan.  
> *Fundamental Limits on the Electron-Phonon Coupling and Superconducting Tc*.  
> [arXiv:2407.12922](https://arxiv.org/abs/2407.12922) (to appear in *Advanced Materials*, 2025)

You may also use the following BibTeX entry:

```bibtex
@article{semenok2025xi,
  author    = {Dmitrii Semenok and Boris L. Altshuler and Emil A. Yuzbashyan},
  title     = {Fundamental Limits on the Electron-Phonon Coupling and Superconducting Tc},
  journal   = {Advanced Materials},
  year      = {2025},
  note      = {arXiv:2407.12922}
}
```

---

## 📬 Contact

For questions or feedback, please contact:
- [dmitrii.semenok@hpstar.ac.cn](mailto:dmitrii.semenok@hpstar.ac.cn)
- [eyuzbash@physics.rutgers.edu](mailto:eyuzbash@physics.rutgers.edu)
