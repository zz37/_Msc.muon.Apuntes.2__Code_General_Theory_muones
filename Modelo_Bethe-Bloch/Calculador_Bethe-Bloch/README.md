# BetheBloch
Implementación simple de solo encabezado de Bethe-Bloch para muones

## Uso
Cree una instancia de `BetheBloch_Calculator`, configurando el material (por ejemplo, `Material::kIron`), y llame a `Calc_dEdx(energy)` donde `energy` es la energía del muón en MeV. El retorno será el `<dE/dx>` de Bethe-Bloch calculado utilizando las descripciones y tablas del Grupo de datos de partículas (https://pdg.lbl.gov/2020/).

## Requisitos
* ROOT >= `v.6.24`.
* Cmake `v.3.14`.

## Ejemplo
Para correlo es:
* `cd <folder_src>`.  
* `mkdir build`.
* `BetheBloch_Example.cpp` se compile en `BetheBloch_example.exe`, que simplemente se ejecuta mediante `./appBetheBloch.exe`.

El ejemplo requiere ROOT (https://root.cern.ch/) para hacer gráficos, aunque las tablas se pueden crear fácilmente evitando las dependencias de ROOT.

* Some of the code does not belong to me, so I made a few changes to improve results. Feel free to use it, at your own discretion.