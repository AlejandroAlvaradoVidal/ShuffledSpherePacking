# Shuffled Sphere Packing en Python/PyTorch

Este proyecto es una traducción del script de MATLAB `GenCSSPv4_NF31` a Python, utilizando PyTorch para las operaciones de tensores y la optimización.

El objetivo es generar patrones de apertura codificada 3D optimizados mediante un algoritmo de empaquetamiento de esferas barajadas (Shuffled Sphere Packing).

## Requisitos

- Python 3.8+
- PyTorch
- NumPy
- Matplotlib
- SciPy

## Instalación

1.  Crea un entorno virtual (recomendado).
2.  Instala las dependencias listadas en `requirements.txt`:
    ```bash
    pip install -r requirements.txt
    ```

## Uso

Ejecuta el script principal desde la terminal. El script iterará de `M=1` a `M=16` por defecto, guardando los patrones y gráficos resultantes en las carpetas `Patrones_pt/` y `graphs_pt/`.

```bash
python main.py
```