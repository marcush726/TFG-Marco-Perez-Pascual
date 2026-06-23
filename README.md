# TFG: Distribuciones cuasiestacionarias en cadenas de Markov en tiempo continuo absorbentes

**Autor:** Marco Pérez Pascual  
**Grado:** Matemáticas y Ciencia de Datos (Universidad Complutense de Madrid)  
**Curso:** 2025-2026  

## Descripción del Proyecto

Este repositorio contiene el código y la documentación correspondientes a mi Trabajo Fin de Grado (TFG). El trabajo estudia la **distribución cuasiestacionaria (QS)** en cadenas de Markov en tiempo continuo (CMTC) absorbentes, una herramienta matemática que describe el comportamiento de un sistema antes de su absorción (por ejemplo, antes de la extinción de una epidemia).

Dado que el análisis límite clásico resulta insuficiente en estos casos (ya que siempre conduce al estado absorbente), se explora la existencia, unicidad y velocidad de convergencia al régimen cuasiestacionario en espacios de estados finitos e infinitos. Además, se introduce la **distribución del cociente de medias (RE)** como alternativa descriptiva.

Toda la base teórica se aplica a tres modelos epidemiológicos clásicos, ilustrados mediante simulaciones numéricas:
* **Modelo SIS** (Susceptible-Infectado-Susceptible)
* **Modelo SIR** (Susceptible-Infectado-Recuperado)
* **Modelo SIRS** (Susceptible-Infectado-Recuperado-Susceptible)

## Estructura del Repositorio

El proyecto está organizado de la siguiente manera:

* `codigo_R/`: Contiene los scripts utilizados para las simulaciones numéricas y la generación de gráficos.
  * `ModeloSIS.R`: Simulación de trayectorias, cálculo de distribuciones y gráficos del modelo SIS.
  * `ModeloSIR.R`: Subgenerador y distribución del cociente de medias (RE) para el modelo SIR.
  * `ModeloSIRS.R`: Análisis del régimen cuasiestacionario y el efecto de la pérdida de inmunidad en el modelo SIRS.
* `Memoria_TFG.pdf`: Documento completo con el desarrollo teórico y los resultados de los experimentos matemáticos. *(Nota: Sube tu PDF final aquí)*.

## Tecnologías Utilizadas

* **R**: Para la simulación estocástica de procesos de Markov, cálculo de autovalores/autovectores de las matrices de transición y visualización de datos.
* **LaTeX**: Para la redacción y maquetación formal de la memoria matemática.