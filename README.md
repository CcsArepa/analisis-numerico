# Proyecto Final – Integración adaptativa mediante la regla de Simpson

## Contenido

| Archivo / carpeta | Descripción |
|---|---|
| `simpson_adaptativo.py` | Implementación **desde cero** de Simpson simple, Simpson compuesto y Simpson adaptativo. |
| `experimentos.py` | Ejecuta todos los experimentos (5 integrales, varias tolerancias) y genera tablas y gráficas. |
| `figuras/` | Gráficas generadas (`.png`). |
| `resultados/tabla_tolerancias.csv` | Tabla completa de resultados. |
| `reporte.tex` | Reporte técnico en LaTeX (4 capítulos), listo para Overleaf. |
| `reporte_overleaf.zip` | `reporte.tex` + `figuras/` comprimidos para subir a Overleaf. |
| `requirements.txt` | Dependencias. |

## Requisitos

- Python 3.9 o superior
- NumPy, Matplotlib, SciPy (SciPy solo se usa como **referencia** de comparación)

```
pip install -r requirements.txt
```

## Cómo ejecutar

Abrir la carpeta en VS Code y, en la terminal:

```
python simpson_adaptativo.py
```
Prueba rápida del algoritmo con ∫₀¹ eˣ dx.

```
python experimentos.py
```
Corre todos los experimentos (tarda ~2 minutos, sobre todo por la búsqueda del número de
evaluaciones que necesita Simpson compuesto para √x). Imprime las tablas en la terminal y guarda
las figuras en `figuras/` y la tabla en `resultados/`. Los resultados son deterministas:
cada ejecución produce exactamente los mismos números.

## Uso de la función

```python
from simpson_adaptativo import simpson_adaptativo
import math

r = simpson_adaptativo(math.sin, 0, math.pi, tol_abs=1e-8, tol_rel=1e-8, max_prof=50)
r["integral"]        # aproximación
r["error_estimado"]  # error estimado por el algoritmo
r["evaluaciones"]    # número de evaluaciones de f
r["intervalos"]      # lista de subintervalos finales [(a1, b1), (a2, b2), ...]
r["prof_max"]        # profundidad máxima alcanzada
r["sin_converger"]   # intervalos aceptados solo por llegar a max_prof
```

## Reporte en Overleaf

1. En Overleaf: **New Project → Upload Project** y subir `reporte_overleaf.zip`.
2. Compilar con pdfLaTeX (opción por defecto).
3. Cambiar `Integrante 1` e `Integrante 2` por los nombres del equipo (línea `\author{...}`).
4. Descargar el PDF.
