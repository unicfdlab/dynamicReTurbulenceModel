# dynamicRekEpsilom Turbulence Model
An improved k-epsilon turbulence model applied to incompressible round jets.
The main advantage of the proposed model is its dependence on a single input parameter - the Reynolds number which is based on the nozzle exit diameter. Other coefficients are determined as a function of the Reynolds number.

# Usage
* Set turbulence model: RASModel dynamicRekEpsilon;
* Set Reynolds number: Re 15000;

# Building

This branch is laid out as an OpenFOAM (openfoam.com) plugin.

- Standalone build into user locations (`$FOAM_USER_LIBBIN`):
  ```bash
  ./Allwmake
  ```
- As a plugin: place (or add as a git submodule) this repository in
  `$WM_PROJECT_DIR/plugins/dynamicReTurbulenceModel` and run `./Allwmake-plugins` from
  `$WM_PROJECT_DIR`. It is then installed into `$FOAM_LIBBIN`.
- Custom location:
  `./Allwmake -prefix=user|group|openfoam|/some/path`

# For citation
If you have found the library useful for your research, please cite next sources:

* Main paper: A. Epikhin, D. Romanova, S. Elistratov; An improved k-e turbulence model applied to incompressible round jets. Physics of Fluids 1 September 2025; 37 (9): 095189. https://doi.org/10.1063/5.0283896
* Validation research: Epikhin, A., Romanova, D. (2026). Validation of Improved k–e Turbulence Model Based on Reynolds Number. Lecture Notes in Mechanical Engineering. Springer, Cham. https://doi.org/10.1007/978-3-032-16138-3_41
