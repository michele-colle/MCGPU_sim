import numpy as np
import xraydb
from scipy.integrate import cumulative_trapezoid
from sqlalchemy import text
import json
import sys
import os

def generate_rayleigh_table(element_symbol='Al', num_points=128):
    xdb = xraydb.get_xraydb()
    q =np.linspace(0, 5, 11)
    array = xraydb.f0('Fe', q)
    print(array)
    # 1. Query the specific columns found in your schema
    # offset = c, scale = [a1...a5], exponents = [b1...b5]
    sql = text("SELECT offset, scale, exponents, atomic_number FROM waasmaier WHERE element = :elem AND ion = ''")
    result = xdb.session.execute(sql, {'elem': element_symbol}).fetchone()
    
    if result is None:
        # Fallback if ion is not an empty string in your DB
        sql = text("SELECT offset, scale, exponents, atomic_number FROM waasmaier WHERE element = :elem LIMIT 1")
        result = xdb.session.execute(sql, {'elem': element_symbol}).fetchone()

    if result is None:
        raise ValueError(f"Element {element_symbol} not found in waasmaier table.")

    # 2. Parse the data
    offset = result[0]
    # The scale and exponents columns are stored as JSON-like strings (e.g., "[6.38, 3.25...]")
    scales = json.loads(result[1])
    exponents = json.loads(result[2])
    z = result[3]

    # 3. Physics Grid: Momentum transfer squared
    q_max = 12.0  
    x_grid = np.linspace(0, np.sqrt(q_max), num_points)**2
    
    # 4. Calculate Squared Form Factor F^2(x)
    def get_f2(q):
        s = q / (4 * np.pi)
        f = offset
        for a, b in zip(scales, exponents):
            f += a * np.exp(-b * s**2)
        return f**2

    f2_values = np.array([get_f2(q) for q in x_grid])
    
    # 5. Cumulative Probability (P)
    p_values = cumulative_trapezoid(f2_values, x_grid, initial=0)
    p_total = p_values[-1]
    p_values /= p_total
    
    # 6. RITA Parameters (A and B)
    a_params, b_params = [], []
    for i in range(num_points):
        if i < num_points - 1:
            dx = x_grid[i+1] - x_grid[i]
            slope_norm = (f2_values[i] / p_total) * dx
            # Curvature logic for RITA sampling in MC-GPU
            A = -0.5 * (1.0 - slope_norm)
            B = 0.01 * (z / 100.0)
        else:
            A, B = 0.0, 0.0
        a_params.append(A)
        b_params.append(B)

    # 7. Search Tree Indices (ITL/ITU)
    itl = [max(1, i) for i in range(1, num_points + 1)]
    itu = [min(num_points, i + 2) for i in range(1, num_points + 1)]

    return x_grid, p_values, a_params, b_params, itl, itu

def generate_file():
    symbol = 'Al' 
    filename = f"rayleigh_{symbol}.txt"
    try:
        x, p, a, b, itl, itu = generate_rayleigh_table(symbol, 128)
        with open(filename, 'w') as f:
            f.write(f"#[RAYLEIGH INTERACTIONS FOR {symbol}]\n")
            f.write(f"#[DATA VALUES TO SAMPLE SQUARED MOLECULAR FORM FACTOR (F^2)]\n")
            f.write(f"   128\n")
            f.write(f"#[SAMPLING DATA FROM COMMON/CGRA/: X, P, A, B, ITL, ITU]\n")

            for i in range(len(x)):
                f.write(f" {x[i]:16.10E} {p[i]:16.10E} {a[i]:16.10E} {b[i]:16.10E} {itl[i]:4} {itu[i]:4}\n")

    except Exception as e:
        import traceback
        print(f"Error: {e}")
        traceback.print_exc()
    
    print(f"Successfully generated: {os.path.abspath(filename)}")

if __name__ == "__main__":
    generate_file()