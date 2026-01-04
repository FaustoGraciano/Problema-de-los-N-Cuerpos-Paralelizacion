# 🌌 Paralelización del Problema de los N Cuerpos

## 📌 Descripción
Este proyecto implementa la simulación del **Problema de los N Cuerpos**, el cual modela el movimiento de un conjunto de cuerpos que interactúan entre sí mediante la **ley de gravitación universal de Newton**.

Dado que el cálculo de fuerzas entre todos los cuerpos tiene una complejidad **O(N²)**, el objetivo principal del trabajo es **reducir el tiempo de ejecución y mejorar la escalabilidad** mediante técnicas de **programación paralela**.

Se desarrollaron y compararon dos versiones paralelas del algoritmo:
- Una implementación basada en **Pthreads** (memoria compartida).
- Una implementación **híbrida MPI + Pthreads** (memoria distribuida y compartida).

---

## ▶️ Cómo ejecutar el programa
1. Compilar el código según la versión deseada:
   - **Pthreads:** usando un compilador compatible con `pthread` (por ejemplo `gcc`).
   - **MPI + Pthreads:** usando un compilador MPI (`mpicc`) con soporte para `pthread`.
2. Ejecutar el programa indicando los parámetros de simulación (cantidad de cuerpos, pasos de tiempo, etc.).
3. Para la versión híbrida, ejecutar utilizando `mpirun` o `mpiexec`, indicando la cantidad de procesos MPI.
4. El programa imprime los tiempos de ejecución y permite comparar resultados con la versión secuencial.

---

## 🧪 Validación
Para validar la correcta paralelización, se comparan las **posiciones finales** de los cuerpos obtenidas con las versiones paralelas frente a la versión secuencial.  
Se calcula la distancia entre posiciones finales y se verifica que el error promedio sea bajo, garantizando resultados consistentes.

---

## ⚙️ Implementación (resumen técnico)
- **Lenguaje:** C
- **Modelos de paralelismo:**
  - **Pthreads:** paralelismo sobre memoria compartida
  - **MPI + Pthreads:** paralelismo híbrido (distribuido + compartido)
- **Estrategias utilizadas:**
  - Paradigma **SPMD**
  - Distribución de trabajo por **stripes**
  - Matrices de fuerzas privadas por hilo
  - **Barreras de sincronización** para evitar condiciones de carrera
  - Comunicación **sistólica** entre procesos MPI
- **Arquitecturas soportadas:**
  - Sistemas de un solo nodo (Pthreads)
  - Sistemas distribuidos multinodo (MPI + Pthreads)

---

## 📊 Resultados
Ambas versiones paralelas muestran una **mejora significativa del tiempo de ejecución** respecto a la versión secuencial.  
El análisis de eficiencia y speedup indica que las soluciones son **débilmente escalables**, ajustándose al modelo de **Gustafson**, especialmente para tamaños de problema grandes.

---

## 📚 Contexto académico
Trabajo Final  
**Sistemas Distribuidos y Paralelos – UNLP**  
Año 2025

---

## 📄 Licencia
Este proyecto se distribuye bajo la **MIT License**.

