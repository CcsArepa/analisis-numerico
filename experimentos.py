"""
experimentos.py
===============
Experimentos numéricos del proyecto (Parte 3).

Ejecutar con:
    python experimentos.py

Genera:
    figuras/    -> todas las gráficas (.png)
    resultados/ -> tablas en formato .csv

Los experimentos son deterministas (no hay números aleatorios),
por lo que al ejecutarlos de nuevo se obtienen exactamente los mismos resultados.

scipy.integrate.quad se usa SOLO como referencia para comparar.
"""

import os
import csv
import math
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg") # guarda las figuras sin abrir ventanas
import matplotlib.pyplot as plt
from scipy.integrate import quad # solo para comparar

from simpson_adaptativo import simpson_adaptativo, simpson_compuesto

os.makedirs("figuras", exist_ok=True)
os.makedirs("resultados", exist_ok=True)

# 1. Funciones a integrar
def f0(x):
    """Ejemplo de clase (diapositiva 15, clase 25): oscila cada vez más rápido cerca de x = 4."""
    return (x + 1) ** 2 * np.cos((2 * x + 1) / (x - 4.3))

def f1(x):
    """Función suave."""
    return np.exp(x)

def f2(x):
    """Función oscilatoria."""
    return np.sin(20 * x)

def f3(x):
    """Función con variación localizada."""
    return np.exp(-100 * (x - 0.5) ** 2)

def f4(x):
    """Función con derivadas no acotadas en x = 0."""
    return np.sqrt(x)


# Valor de referencia para I0: no tiene primitiva elemental, así que se usa
# scipy.integrate.quad con una tolerancia muy estricta (solo como referencia).
# (quad avisa que 1e-14 está al límite del redondeo; se ignora ese aviso.
#  El valor se confirma con Simpson adaptativo a tol = 1e-10: difieren en ~2e-13.)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    I0_ref = quad(f0, 0, 4, epsabs=1e-14, epsrel=1e-14, limit=500)[0]

# Cada problema: (nombre corto, título, función, a, b, valor exacto, fórmula LaTeX)
problemas = [
    ("I0", "Ejemplo de clase", f0, 0.0, 4.0, I0_ref,
     r"$(x+1)^2\cos\left(\frac{2x+1}{x-4.3}\right)$"),
    ("I1", "Función suave", f1, 0.0, 1.0, math.e - 1,
     r"$e^{x}$"),
    ("I2", "Función oscilatoria", f2, 0.0, 2.0, (1 - math.cos(40)) / 20,
     r"$\sin(20x)$"),
    ("I3", "Variación localizada", f3, 0.0, 1.0, math.sqrt(math.pi) / 10 * math.erf(5),
     r"$e^{-100(x-0.5)^2}$"),
    ("I4", "Derivadas no acotadas", f4, 0.0, 1.0, 2.0 / 3.0,
     r"$\sqrt{x}$"),
]

# Tolerancias que se van a probar (se usa la misma para la parte absoluta y la relativa)
tolerancias = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
TOL_GRAFICA = 1e-6     # tolerancia usada para dibujar los subintervalos
MAX_PROF = 50


def evaluaciones_compuesto_necesarias(f, a, b, exacto, error_objetivo):
    """Busca el menor n (par) tal que Simpson compuesto alcance error <= error_objetivo.
    Primero duplica n y luego hace bisección. Devuelve el número de evaluaciones (n + 1)."""
    n = 2
    while abs(simpson_compuesto(f, a, b, n)[0] - exacto) > error_objetivo:
        n *= 2
        if n > 2 ** 24:
            return None
    bajo, alto = n // 2, n # el n buscado está entre bajo y alto
    while alto - bajo > 2:
        medio = (bajo + alto) // 2
        if medio % 2 == 1:
            medio += 1
        if abs(simpson_compuesto(f, a, b, medio)[0] - exacto) <= error_objetivo:
            alto = medio
        else:
            bajo = medio
    return alto + 1

# 2. Experimento A: distintas tolerancias, comparación con compuesto y quad
print("=" * 100)
print("Valor de referencia de I0 (quad, tol = 1e-14):", repr(I0_ref))
print("=" * 100)

filas = []
for nombre, titulo, f, a, b, exacto, _ in problemas:
    print(f"\n{nombre}: {titulo}   (valor exacto = {exacto:.15f})")
    print(f"{'tol':>8} {'adapt. error':>13} {'err. estim.':>12} {'eval':>6} {'subint':>6} {'prof':>5}"
          f" | {'comp. mismo n':>13} | {'eval comp.':>10} | {'quad error':>11} {'eval':>5}")
    for tol in tolerancias:
        # Simpson adaptativo
        r = simpson_adaptativo(f, a, b, tol_abs=tol, tol_rel=tol, max_prof=MAX_PROF)
        err_adapt = abs(r["integral"] - exacto)

        # Simpson compuesto con (aprox.) el mismo número de evaluaciones
        n_mismo = r["evaluaciones"] - 1
        if n_mismo % 2 == 1:
            n_mismo += 1
        err_comp_mismo = abs(simpson_compuesto(f, a, b, n_mismo)[0] - exacto)

        # Evaluaciones que necesita el compuesto para lograr el mismo error que el adaptativo
        objetivo = max(err_adapt, 1e-15)
        eval_comp = evaluaciones_compuesto_necesarias(f, a, b, exacto, objetivo)

        # scipy.integrate.quad (solo referencia)
        q, _, info_q = quad(f, a, b, epsabs=tol, epsrel=tol, limit=200, full_output=1)
        err_quad = abs(q - exacto)

        filas.append({
            "problema": nombre, "tol": tol,
            "adapt_integral": r["integral"], "adapt_error": err_adapt,
            "adapt_error_estimado": r["error_estimado"],
            "adapt_evaluaciones": r["evaluaciones"],
            "adapt_subintervalos": len(r["intervalos"]),
            "adapt_prof_max": r["prof_max"], "adapt_sin_converger": r["sin_converger"],
            "comp_mismas_eval_error": err_comp_mismo,
            "comp_eval_para_mismo_error": eval_comp,
            "quad_error": err_quad, "quad_evaluaciones": info_q["neval"],
        })
        print(f"{tol:8.0e} {err_adapt:13.3e} {r['error_estimado']:12.3e} {r['evaluaciones']:6d}"
              f" {len(r['intervalos']):6d} {r['prof_max']:5d} | {err_comp_mismo:13.3e} |"
              f" {str(eval_comp):>10} | {err_quad:11.3e} {info_q['neval']:5d}")

with open("resultados/tabla_tolerancias.csv", "w", newline="", encoding="utf-8") as archivo:
    escritor = csv.DictWriter(archivo, fieldnames=list(filas[0].keys()))
    escritor.writeheader()
    escritor.writerows(filas)

# 3. Experimento B: error en función del número de evaluaciones
fig, ejes = plt.subplots(3, 2, figsize=(11, 12))
ejes = ejes.flatten()
for k, (nombre, titulo, f, a, b, exacto, formula) in enumerate(problemas):
    # Simpson compuesto: n = 2, 4, 8, ..., 2^17
    ev_c, er_c = [], []
    for p in range(1, 18):
        valor, nev = simpson_compuesto(f, a, b, 2 ** p)
        ev_c.append(nev)
        er_c.append(max(abs(valor - exacto), 1e-17))

    # Simpson adaptativo y quad para varias tolerancias
    ev_a, er_a, ev_q, er_q = [], [], [], []
    for tol in [10.0 ** (-p) for p in range(1, 13)]:
        r = simpson_adaptativo(f, a, b, tol_abs=tol, tol_rel=tol, max_prof=MAX_PROF)
        ev_a.append(r["evaluaciones"])
        er_a.append(max(abs(r["integral"] - exacto), 1e-17))
        q, _, info_q = quad(f, a, b, epsabs=tol, epsrel=tol, limit=200, full_output=1)
        ev_q.append(info_q["neval"])
        er_q.append(max(abs(q - exacto), 1e-17))

    ax = ejes[k]
    ax.loglog(ev_c, er_c, "s-", label="Simpson compuesto")
    ax.loglog(ev_a, er_a, "o-", label="Simpson adaptativo")
    ax.loglog(ev_q, er_q, "^--", label="scipy quad (referencia)")
    ax.set_title(f"{nombre}: {formula}")
    ax.set_xlabel("Número de evaluaciones de f")
    ax.set_ylabel("Error absoluto")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
ejes[5].axis("off")
fig.tight_layout()
fig.savefig("figuras/error_vs_evaluaciones.png", dpi=150)
plt.close(fig)

# 4. Experimento C: tolerancia pedida vs error real y vs evaluaciones
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
for nombre, titulo, f, a, b, exacto, formula in problemas:
    datos = [fila for fila in filas if fila["problema"] == nombre]
    tols = [d["tol"] for d in datos]
    ax1.loglog(tols, [max(d["adapt_error"], 1e-17) for d in datos], "o-", label=nombre)
    ax2.loglog(tols, [d["adapt_evaluaciones"] for d in datos], "o-", label=nombre)
ax1.loglog(tolerancias, tolerancias, "k--", label="error = tol")
ax1.set_xlabel("Tolerancia solicitada")
ax1.set_ylabel("Error real")
ax1.set_title("Error real vs tolerancia (Simpson adaptativo)")
ax2.set_xlabel("Tolerancia solicitada")
ax2.set_ylabel("Evaluaciones de f")
ax2.set_title("Costo vs tolerancia (Simpson adaptativo)")
for ax in (ax1, ax2):
    ax.invert_xaxis()
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("figuras/tolerancia_error_costo.png", dpi=150)
plt.close(fig)

# 5. Experimento D: gráfica de cada función con sus subintervalos
for nombre, titulo, f, a, b, exacto, formula in problemas:
    r = simpson_adaptativo(f, a, b, tol_abs=TOL_GRAFICA, tol_rel=TOL_GRAFICA, max_prof=MAX_PROF)
    intervalos = r["intervalos"]
    nodos = [intervalos[0][0]] + [fin for (_, fin) in intervalos]

    x = np.linspace(a, b, 2000)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1.3]})
    # Arriba: función + bordes de los subintervalos
    for xn in nodos:
        ax1.axvline(xn, color="tab:red", lw=0.5, alpha=0.5)
    ax1.plot(x, f(x), color="tab:blue", lw=1.8, label=f"f(x) = {formula}")
    ax1.plot(nodos, f(np.array(nodos)), "o", color="tab:red", ms=2.5,
             label=f"bordes de subintervalos ({len(intervalos)} subintervalos)")
    ax1.set_ylabel("f(x)")
    ax1.set_title(f"{nombre} - {titulo}  (tol = {TOL_GRAFICA:.0e}, "
                  f"{r['evaluaciones']} evaluaciones)")
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(alpha=0.3)
    # Abajo: ancho de cada subintervalo
    centros = [(i + j) / 2 for (i, j) in intervalos]
    anchos = [j - i for (i, j) in intervalos]
    ax2.semilogy(centros, anchos, "o-", color="tab:green", ms=3)
    ax2.set_xlabel("x")
    ax2.set_ylabel("ancho h")
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"figuras/subintervalos_{nombre}.png", dpi=150)
    plt.close(fig)

# Las 5 funciones juntas (vista general)
fig, ejes = plt.subplots(3, 2, figsize=(11, 10))
ejes = ejes.flatten()
for k, (nombre, titulo, f, a, b, exacto, formula) in enumerate(problemas):
    x = np.linspace(a, b, 2000)
    ejes[k].plot(x, f(x))
    ejes[k].set_title(f"{nombre}: {formula} en [{a:g}, {b:g}]")
    ejes[k].grid(alpha=0.3)
ejes[5].axis("off")
fig.tight_layout()
fig.savefig("figuras/funciones.png", dpi=150)
plt.close(fig)

# 6. Resumen con tol = 1e-6
print("\n" + "=" * 100)
print(f"RESUMEN con tol = {TOL_GRAFICA:.0e}")
print("=" * 100)
for fila in filas:
    if fila["tol"] == TOL_GRAFICA:
        print(f"{fila['problema']}: integral = {fila['adapt_integral']:.12f}, "
              f"error = {fila['adapt_error']:.2e}, eval = {fila['adapt_evaluaciones']}, "
              f"subintervalos = {fila['adapt_subintervalos']}, "
              f"compuesto necesita {fila['comp_eval_para_mismo_error']} eval, "
              f"quad usa {fila['quad_evaluaciones']} eval")

print("\nFiguras guardadas en 'figuras/' y tablas en 'resultados/'.")
