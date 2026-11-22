"""
simpson_adaptativo.py
=====================
Implementación DESDE CERO de:

  * simpson_simple      -> regla de Simpson en un solo intervalo [a, b]
  * simpson_compuesto   -> regla de Simpson compuesta con partición uniforme
  * simpson_adaptativo  -> integración adaptativa basada en la regla de Simpson

No se usa scipy.integrate en este archivo. Solo Python puro (y math).

Idea del método adaptativo
--------------------------
En cada intervalo [a, b] se calculan dos aproximaciones:
    S1 = Simpson con un solo paso   h = (b - a) / 2
    S2 = Simpson con dos pasos      h = (b - a) / 4   (mitad izquierda + mitad derecha)
El error de S2 se estima como  |S2 - S1| / 15.
  - Si esa estimación cumple la tolerancia, el intervalo se ACEPTA.
  - Si no, el intervalo se divide en dos mitades y se repite el proceso en cada una.
"""

import math


# ---------------------------------------------------------------------------
# 1. Regla de Simpson simple
# ---------------------------------------------------------------------------
def simpson_simple(a, b, fa, fm, fb):
    """Regla de Simpson en [a, b] usando los valores ya calculados
    fa = f(a), fm = f((a+b)/2), fb = f(b).

    S = (b - a)/6 * [ f(a) + 4 f(m) + f(b) ]
    """
    return (b - a) / 6.0 * (fa + 4.0 * fm + fb)


# ---------------------------------------------------------------------------
# 2. Regla de Simpson compuesta (partición uniforme)
# ---------------------------------------------------------------------------
def simpson_compuesto(f, a, b, n):
    """Simpson compuesto con n subintervalos iguales (n debe ser par).

    Devuelve (aproximación, número de evaluaciones de f).
    Usa n + 1 evaluaciones de la función.
    """
    if n % 2 != 0:
        n = n + 1  # Simpson necesita un número par de subintervalos
    h = (b - a) / n

    suma = f(a) + f(b)
    for i in range(1, n):
        x = a + i * h
        if i % 2 == 1:
            suma += 4.0 * f(x)   # puntos impares
        else:
            suma += 2.0 * f(x)   # puntos pares
    return h / 3.0 * suma, n + 1


# ---------------------------------------------------------------------------
# 3. Simpson adaptativo
# ---------------------------------------------------------------------------
def simpson_adaptativo(f, a, b, tol_abs=1e-6, tol_rel=1e-6, max_prof=50):
    """Integración adaptativa de Simpson.

    Parámetros
    ----------
    f        : función a integrar, f(x) -> float
    a, b     : extremos del intervalo (a < b)
    tol_abs  : tolerancia absoluta global
    tol_rel  : tolerancia relativa
    max_prof : profundidad máxima de recursión (límite de subdivisiones)

    Criterio de aceptación de un subintervalo [a_i, b_i] de profundidad k:

        |S2 - S1| / 15  <=  max( tol_abs / 2^k ,  tol_rel * |S2| )

    - La parte absoluta se reparte entre las dos mitades (se divide entre 2
      en cada nivel), así la suma de los errores locales no supera tol_abs.
    - La parte relativa compara el error con el tamaño de la integral local.

    Devuelve un diccionario con:
        'integral'        : aproximación de la integral
        'error_estimado'  : suma de los errores locales estimados
        'evaluaciones'    : número total de evaluaciones de f
        'intervalos'      : lista de subintervalos finales aceptados (a_i, b_i)
        'prof_max'        : profundidad máxima alcanzada
        'sin_converger'   : subintervalos aceptados solo por llegar a max_prof
    """
    # Diccionario donde se van guardando los registros del algoritmo
    info = {
        "evaluaciones": 0,
        "intervalos": [],
        "error_estimado": 0.0,
        "prof_max": 0,
        "sin_converger": 0,
    }

    def evaluar(x):
        """Evalúa f(x) y cuenta la evaluación."""
        info["evaluaciones"] += 1
        return f(x)

    def paso_recursivo(a, b, fa, fm, fb, S1, tol_a, prof):
        """Procesa el intervalo [a, b].

        fa, fm, fb y S1 ya vienen calculados del nivel anterior,
        por eso solo se necesitan 2 evaluaciones nuevas por intervalo.
        """
        m = (a + b) / 2.0
        ml = (a + m) / 2.0          # punto medio de la mitad izquierda
        mr = (m + b) / 2.0          # punto medio de la mitad derecha
        fml = evaluar(ml)           # <- únicas evaluaciones nuevas
        fmr = evaluar(mr)

        S_izq = simpson_simple(a, m, fa, fml, fm)
        S_der = simpson_simple(m, b, fm, fmr, fb)
        S2 = S_izq + S_der

        error = abs(S2 - S1) / 15.0                 # estimación del error de S2
        tolerancia = max(tol_a, tol_rel * abs(S2))  # tolerancia absoluta/relativa

        info["prof_max"] = max(info["prof_max"], prof)

        if error <= tolerancia or prof >= max_prof:
            # Intervalo ACEPTADO
            if error > tolerancia:
                info["sin_converger"] += 1
            info["intervalos"].append((a, b))
            info["error_estimado"] += error
            # Extrapolación de Richardson: mejora S2 usando la estimación del error
            return S2 + (S2 - S1) / 15.0

        # Intervalo RECHAZADO: se divide en dos y se reutilizan fa, fml, fm, fmr, fb
        izquierda = paso_recursivo(a, m, fa, fml, fm, S_izq, tol_a / 2.0, prof + 1)
        derecha = paso_recursivo(m, b, fm, fmr, fb, S_der, tol_a / 2.0, prof + 1)
        return izquierda + derecha

    # Evaluaciones iniciales: extremos y punto medio
    fa = evaluar(a)
    fb = evaluar(b)
    fm = evaluar((a + b) / 2.0)
    S1 = simpson_simple(a, b, fa, fm, fb)

    integral = paso_recursivo(a, b, fa, fm, fb, S1, tol_abs, 0)

    info["integral"] = integral
    info["intervalos"].sort()
    return info


# ---------------------------------------------------------------------------
# Pequeña prueba al ejecutar este archivo directamente
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    r = simpson_adaptativo(math.exp, 0.0, 1.0, tol_abs=1e-8, tol_rel=1e-8)
    exacto = math.e - 1
    print("Integral de e^x en [0, 1]")
    print("  Simpson adaptativo :", r["integral"])
    print("  Valor exacto       :", exacto)
    print("  Error real         :", abs(r["integral"] - exacto))
    print("  Evaluaciones de f  :", r["evaluaciones"])
    print("  Subintervalos      :", len(r["intervalos"]))
