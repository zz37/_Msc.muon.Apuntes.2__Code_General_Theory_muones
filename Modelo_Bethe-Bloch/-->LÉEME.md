# Archivos en el Folder
## Teoría de Bethe-Bloch

* Teóría de Bethe-Bloch.

* **NOTA**: Hay variantes de la forma de Bethe-Bloch donde algunos escritos donde se grafica en función de la E (energía) , p (momento) o en este caso β․ɣ (unidades de momento MeV/c). 
* **NOTA**: Para más información ver la notas de ONEnote.


### Cuadernos

* En los cuadernos de *jupyer.ipynb* y *mathematica.nb* se escribe la forma de Bethe-Bloch como 𝑓 -->  f(β․ɣ). Para un muon y material de elemento Silicio.
    * `Jupyter.ipynb`, `v. Python3.10`.
        * Calcula la formula de Bethe-Bloch $\frac{dE}{dx}$, para el término simple, correcciones de densidad (parametrización Sternheimer’s) y sus gráficas en función $𝑓 (\beta \cdot \gamma)$. 

    * `Mathematica.nb`, `v. 13.0`.
        * Calcula la formula de Bethe-Bloch $\frac{dE}{dx}$, para el término simple, que no incluye las correcciones ni efecto de densidad. Es para observar la forma de como escribir $\frac{dE}{dx}$ --> $𝑓 (\beta \cdot \gamma)$.

    * Ambos códigos incluyen la forma de Bethe-Bloch con el denominador del argumento del $Ln () --> I²ᐧ(1-β²)$.
    * El archivo *C++* calcula la forma de Bethe-Bloch con correcciones pero en función de la energía. Y el término del denominador es  I².