# Circular-Segment Quadrature Rules (QR, QS, and QR*)

Reference implementations and benchmark datasets accompanying the manuscript:  
**"Circular-segment quadrature rules based on the three-point circumcircle"**  
*Authors: Md. Abdur Rahman, Md. Sindid Safwanur Kouly, Dr. Mostak Ahmed*

## Repository Contents
- `qr_quadrature.py`: Zero-dependency IEEE 754 double-precision Python implementation of single-panel and composite $Q_R$, simplified algebraic $Q_S$, Richardson-extrapolated $Q_R^*$, and equal-cost classical Newton–Cotes rules.
- `circular_quadrature_30digit.wl` & `circular_quadrature_30digit.nb`: 30-digit arbitrary-precision Wolfram Mathematica implementation including Lemma 1 graph-validity checks, Maclaurin series evaluation for small curvature ($t < 1/4$), and convergence order analysis.
- `Supplementary Material - 100-Function Quadrature Benchmark.pdf`: Complete benchmark dataset and family-wise performance summary across 100 test integrands.

## Usage
- **Python:** Run `python qr_quadrature.py`
- **Mathematica:** Open `circular_quadrature_30digit.nb` and press `Shift + Enter` to evaluate all cells.
