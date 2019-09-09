# sprayFlameletFoam
This solver implements the **Flamelet Progress-Variable (FPV)** approach for modeling spray combustion in [OpenFOAM-6](https://openfoam.org/version/6/). Detailed chemical processes are mapped onto two trajectory variables, namely the mixture fraction and the progress variable.

## Installation
```bash
mkdir ~/OpenFOAM/sprayFlameletFoam/
git clone https://github.com/ZX114/sprayFlameletFoam.git ~/OpenFOAM/sprayFlameletFoam/
cd ~/OpenFOAM/sprayFlameletFoam/
./Allwmake
```

## Example -- aachenBomb
- Flamelet tables can be generated from solutions of the counterflow diffusion flame. And the laminar results are stored in `tables/Zeta_0/`.
- To account for the effect of turbulence, a presumed beta-PDF integration is then performed. This is done through `tableTest` provided in the `utilities` directory, generating flamelet data at different *Zeta*.
- With the complete set of tablse, run **sprayFlameletFoam** with
```bash
sprayFlameletFoam > log &
```
