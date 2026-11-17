"""
Circular-Segment Quadrature Rules (QR, QS, and Extrapolated QR*)
Zero-Dependency Reference Implementation in Python (IEEE 754 Double Precision)

Paper: Circular-segment quadrature rules based on the three-point circumcircle
Authors: Md. Abdur Rahman, Md. Sindid Safwanur Kouly, Dr. Mostak Ahmed
"""

import math

PHI_COEFFS = (
    2.0 / 3.0,
    1.0 / 5.0,
    3.0 / 28.0,
    5.0 / 72.0,
    35.0 / 704.0,
    63.0 / 1664.0,
    231.0 / 7680.0,
    429.0 / 17408.0,
    6435.0 / 311296.0,
    12155.0 / 688128.0,
    46189.0 / 3014656.0,
    88179.0 / 6553600.0,
    676039.0 / 56623104.0,
    1300075.0 / 121634816.0,
)


def check_graph_validity(h: float, y0: float, y1: float, y2: float, eps: float = 1e-12) -> bool:
    """
    Lemma 1 (Graph Criterion): Checks whether the circumcircle arc through
    (a, y0), (m, y1), (b, y2) is the graph of a single-valued function over [a, b].
    Condition: D^2 <= 2 * (h^2 + d_min^2) * (1 + eps).
    """
    d = y0 - 2.0 * y1 + y2
    if d == 0.0:
        return True
    sigma_d = 1.0 if d > 0.0 else -1.0
    d_min = min(sigma_d * (y0 - y1), sigma_d * (y2 - y1))
    return (d * d) <= 2.0 * (h * h + d_min * d_min) * (1.0 + eps)


def qr_exact_panel(f, a: float, b: float, warn_invalid: bool = False) -> float:
    """
    Algorithm 1 (Exact Rule QR): Evaluates the circular-segment quadrature rule
    on a single panel [a, b] using 3 equispaced samples.
    """
    h = 0.5 * (b - a)
    m = 0.5 * (a + b)
    y0, y1, y2 = float(f(a)), float(f(m)), float(f(b))

    d = y0 - 2.0 * y1 + y2
    trap = h * (y0 + y2)

    if d == 0.0:
        return trap

    if warn_invalid and not check_graph_validity(h, y0, y1, y2):
        raise ValueError(f"Panel [{a}, {b}] violates Lemma 1 graph-validity criterion.")

    l1 = math.hypot(h, y1 - y0)
    l2 = math.hypot(h, y2 - y1)
    c2 = 4.0 * h * h + (y2 - y0) ** 2
    t = min(1.0, (h * abs(d)) / (l1 * l2))

    # Stable evaluation of segment area S (Section 2.5)
    if t < 0.25:
        t2 = t * t
        poly = 0.0
        for coeff in reversed(PHI_COEFFS):
            poly = poly * t2 + coeff
        seg = 0.25 * c2 * t * poly
    else:
        theta = 2.0 * math.asin(t)
        r = math.sqrt(c2) / (2.0 * t)
        seg = 0.5 * r * r * (theta - math.sin(theta))

    sigma = -1.0 if d > 0.0 else 1.0
    return trap + sigma * seg


def qr_simplified_panel(f, a: float, b: float) -> float:
    """
    Algorithm 1 (Simplified Algebraic Rule QS): Evaluates the algebraic
    circular-arc rule on a single panel [a, b] without trigonometric functions.
    """
    h = 0.5 * (b - a)
    m = 0.5 * (a + b)
    y0, y1, y2 = float(f(a)), float(f(m)), float(f(b))

    d = y0 - 2.0 * y1 + y2
    trap = h * (y0 + y2)
    if d == 0.0:
        return trap

    l1 = math.hypot(h, y1 - y0)
    l2 = math.hypot(h, y2 - y1)
    c2 = 4.0 * h * h + (y2 - y0) ** 2
    mu = c2 / (l1 * l2)

    return trap - (h * mu * d) / 6.0


# =========================================================================
# Composite QR, QS & Richardson Extrapolation
# =========================================================================

def composite_qr(f, A: float, B: float, n_panels: int) -> float:
    """Composite QR rule over N panels (2N + 1 function evaluations)."""
    h_panel = (B - A) / n_panels
    return sum(
        qr_exact_panel(f, A + i * h_panel, A + (i + 1) * h_panel)
        for i in range(n_panels)
    )


def composite_qs(f, A: float, B: float, n_panels: int) -> float:
    """Composite QS rule over N panels (2N + 1 function evaluations)."""
    h_panel = (B - A) / n_panels
    return sum(
        qr_simplified_panel(f, A + i * h_panel, A + (i + 1) * h_panel)
        for i in range(n_panels)
    )


def extrapolated_qr(f, A: float, B: float, n_panels: int) -> float:
    """
    Order-6 Richardson-extrapolated rule QR* = (16 * QR_{2N} - QR_N) / 15.
    Uses 4N + 1 function evaluations (Section 5.4).
    """
    q_n = composite_qr(f, A, B, n_panels)
    q_2n = composite_qr(f, A, B, 2 * n_panels)
    return (16.0 * q_2n - q_n) / 15.0


# =========================================================================
# Classical Newton-Cotes Rules (Equal Nodes: 2 * n_panels + 1)
# =========================================================================

def composite_trapezoid(f, A: float, B: float, n_panels: int) -> float:
    m_sub = 2 * n_panels
    hh = (B - A) / m_sub
    y = [float(f(A + i * hh)) for i in range(m_sub + 1)]
    return hh * (0.5 * y[0] + sum(y[1:-1]) + 0.5 * y[-1])


def composite_simpson13(f, A: float, B: float, n_panels: int) -> float:
    m_sub = 2 * n_panels
    hh = (B - A) / m_sub
    y = [float(f(A + i * hh)) for i in range(m_sub + 1)]
    return (hh / 3.0) * (y[0] + 4.0 * sum(y[1:-1:2]) + 2.0 * sum(y[2:-2:2]) + y[-1])


def composite_simpson38(f, A: float, B: float, n_panels: int) -> float:
    m_sub = 2 * n_panels
    if m_sub % 3 != 0:
        raise ValueError("2 * n_panels must be divisible by 3 for Simpson 3/8 rule.")
    hh = (B - A) / m_sub
    y = [float(f(A + i * hh)) for i in range(m_sub + 1)]
    total = sum(
        y[i] + 3.0 * y[i + 1] + 3.0 * y[i + 2] + y[i + 3]
        for i in range(0, m_sub, 3)
    )
    return (3.0 * hh / 8.0) * total


def composite_boole(f, A: float, B: float, n_panels: int) -> float:
    m_sub = 2 * n_panels
    if m_sub % 4 != 0:
        raise ValueError("2 * n_panels must be divisible by 4 for Boole's rule.")
    hh = (B - A) / m_sub
    y = [float(f(A + i * hh)) for i in range(m_sub + 1)]
    total = sum(
        7.0 * y[i] + 32.0 * y[i + 1] + 12.0 * y[i + 2] + 32.0 * y[i + 3] + 7.0 * y[i + 4]
        for i in range(0, m_sub, 4)
    )
    return (2.0 * hh / 45.0) * total


def composite_weddle(f, A: float, B: float, n_panels: int) -> float:
    m_sub = 2 * n_panels
    if m_sub % 6 != 0:
        raise ValueError("2 * n_panels must be divisible by 6 for Weddle's rule.")
    hh = (B - A) / m_sub
    y = [float(f(A + i * hh)) for i in range(m_sub + 1)]
    total = sum(
        y[i] + 5.0 * y[i + 1] + y[i + 2] + 6.0 * y[i + 3] + y[i + 4] + 5.0 * y[i + 5] + y[i + 6]
        for i in range(0, m_sub, 6)
    )
    return (3.0 * hh / 10.0) * total


# =========================================================================
# All-in-One Comparison Runner
# =========================================================================

def compare_all_methods(title: str, f, A: float, B: float, exact_val: float, n_panels: int = 12):
    methods = [
        ("Composite Trapezoid", composite_trapezoid(f, A, B, n_panels)),
        ("Composite Simpson 1/3", composite_simpson13(f, A, B, n_panels)),
        ("Composite Simpson 3/8", composite_simpson38(f, A, B, n_panels)),
        ("Composite Boole", composite_boole(f, A, B, n_panels)),
        ("Composite Weddle", composite_weddle(f, A, B, n_panels)),
        ("Composite QR (Exact)", composite_qr(f, A, B, n_panels)),
        ("Composite QS (Simplified)", composite_qs(f, A, B, n_panels)),
        ("Extrapolated QR* (4N+1)", extrapolated_qr(f, A, B, n_panels)),
    ]

    print(f"\n=== {title} on [{A}, {B}] at N = {n_panels} ({2*n_panels + 1} nodes) ===")
    print(f"{'Method':<27} | {'Computed Value':<20} | {'Absolute Error'}")
    print("-" * 68)
    for name, val in methods:
        err = abs(val - exact_val)
        print(f"{name:<27} | {val:<20.15f} | {err:.3e}")


if __name__ == "__main__":
    N = 12  # 25 equispaced nodes (multiple of 6 to run all rules simultaneously)

    # Test 1: Smooth analytic integrand (Function 20, Section 7.1)
    compare_all_methods("Function 20: exp(x)", math.exp, 0.0, 1.0, math.e - 1.0, N)

    # Test 2: Circular arc with endpoint singularity (Function 1, Table 5)
    f_arc = lambda x: math.sqrt(max(0.0, 1.0 - x * x))
    compare_all_methods("Function 1: sqrt(1 - x^2)", f_arc, 0.0, 1.0, math.pi / 4.0, N)
