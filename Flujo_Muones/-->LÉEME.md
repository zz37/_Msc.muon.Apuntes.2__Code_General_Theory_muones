

## Código de flujo de muones

* Calculos del flujo de muones con `python v3.0`, `c++`, `ROOT v6.24/06` y `Cmake v3.0`.


* `FlujoM.cpp` calcula el flujo de muones según la formula de Gaisser.
	* Para correrlo, es con `CMAKE`, usando el `CMakeLists.txt`:
        1. `mkdir build`
	    2. `cd build`
	    3. `cmake ..`
	    4. `make`
	    5. `./RootMain`
	

* `Flujo_Muones.ipynb`. Es el calculo del flujo de muones con `jupyter-notebook` y `python v3.0`.
	* Incluye el cómputo de la forma de Gaisser y la versión modificada de Gaisser/Tang.
	* Falta terminar de implementar el flujo integral que marca errores.

* `IntegratedMuonFlux.cpp`, calcula el flujo integral de muones.
    * Para correrlo con `ROOT` es:   
        1.`root -l IntegratedMuonFlux.cpp`
        2. Ya en ROOT, `muonFlux3()`
        
        
* `ModGaisser.cpp`, calcula el flujo de muones, con la forma de Gaisser modificada.
    * Para correrlo con `ROOT` es:   
        1.`root -l ModGaisser.cpp`
        2. Ya en ROOT, `muonFlux2()`
        

