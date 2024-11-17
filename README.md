# SCIAN-Force (English)

**SCIAN-Force** is software designed to analyze the deformation of Biocompatible Microdroplets (BMD) caused by cellular forces, enabling the quantification of anisotropic stresses on their surfaces. This non-invasive analysis is essential for studying cell migration processes, particularly in tumor progression and metastasis research.

## References

[![SCIAN-Drop & SCIAN-Force: Quantification of Membrane Deformation via Image Analysis of Biocompatible Fluorocarbon MicroDroplets](https://img.shields.io/badge/Poster-PDF-blue)](https://github.com/L3N73J4M4N/SCIAN-Force/releases/download/SCIAN-Force/2024_Poster_MCarvajal_SCIANForce_SCIANDrop.pdf)

**[1]** M. Carvajal et al., "SCIAN-Drop & SCIAN-Force: Quantification of Membrane Deformation via Image Analysis of Biocompatible Fluorocarbon MicroDroplets," XXXVI REUNIÓN ANUAL SOCIEDAD DE BIOLOGÍA CELULAR DE CHILE, 2024. Available in: [Poster SCIAN-Drop & SCIAN-Force](https://github.com/L3N73J4M4N/SCIAN-Force/releases/download/SCIAN-Force/2024_Poster_MCarvajal_SCIANForce_SCIANDrop.pdf).

**[2]** R. Lucio, A. Arashiro, and N. R. Demarquette, "Pendant drop method: A generalization of the Laplace equation to measure interfacial tension," in *Materials Research*, vol. 18, no. 1, pp. 23–32, 2015. DOI: [10.1590/1516-1439.267314](https://doi.org/10.1590/1516-1439.267314).

**[3]** O. Càmpas et al., "Quantifying cell-generated forces within living embryonic tissues," *Nature Methods*, vol. 11, no. 2, pp. 183–189, 2014. DOI: [10.1038/nmeth.2761](https://doi.org/10.1038/nmeth.2761).

**[4]** S. Rusinkiewicz, "Estimating Curvatures and Their Derivatives on Triangle Meshes," in *Proceedings of the IEEE Symposium on 3D Data Processing, Visualization, and Transmission*, Thessaloniki, Greece, 2004, pp. 486–493. DOI: [10.1109/TDPVT.2004.1335277](https://doi.org/10.1109/TDPVT.2004.1335277).

**[5]** M. Mayer et al., "Visualization of fluid dynamics using triangular meshes and curvature-driven smoothing," in *Mathematics and Visualization*, Springer, 2003, pp. 191–205. DOI: [10.1007/978-3-642-55566-6](https://doi.org/10.1007/978-3-642-55566-6).
   
## Description

SCIAN-Force performs advanced analyses, including:
- **3D Reconstruction** of microdroplet surfaces using the `Mesh3D` class.
- **Calculation of average curvatures and anisotropic stresses** using the `Forces3D` class.
- **Interactive visualization of results** using the `Inter3D` class.
- Various utility functions to work with triangular meshes.

## System Requirements

To run the software, ensure you have Python 3.x installed along with the following libraries:

- NumPy
- Meshfix
- Rtree
- Plotly
- SciPy
- Pymeshfix
- Pandas
- Tkinter
- Openpyxl
- Pillow (PIL)
- Skimage
- Tifffile

Install the dependencies using:
```bash
pip install -r requirements.txt
```

## Key Features

- **3D Reconstruction**: Generates 3D models from segmented images.
- **Stress Calculation**: Determines anisotropic pressure at each surface point.
- **Interactive Interface**: Allows intuitive parameter adjustment and result exportation.
- **Experimental Validation**: Tested with synthetic images and experimental data.

## Applications

- Study of tumor progression and metastasis.
- Quantification of forces exerted by cells in 3D environments.
- Analysis of interfacial stresses in biocompatible fluids.

## Usage

1. **Input:** Load a set of segmented images (stacks) using `Mesh3D`, specifying the voxel size. This step is optional if the 3D reconstruction is obtained by other means.
2. **Analysis:** SCIAN-Force calculates anisotropic stresses using `Forces3D`.
3. **Visualization:** Explore the interactive results, including stress and curvature maps.
4. **Exportation:** Save the results in compatible formats such as spreadsheets.

## Download SCIAN-Force

You can download the latest version of **SCIAN-Force** for Windows as an executable file. Click the link below to start the download:

[Download SCIAN-Force.exe (updated 2024/11/16)](https://github.com/L3N73J4M4N/SCIAN-Force/releases/download/SCIAN-Force/SCIAN-Force.v2024.11.08.exe)

## License

This project is licensed under the MIT License.

---

# SCIAN-Force (Español)

**SCIAN-Force** es un software diseñado para analizar la deformación de Microgotas Biocompatibles (BMD) causadas por fuerzas celulares, permitiendo cuantificar las tensiones anisotrópicas en su superficie. Este análisis no invasivo es clave en la investigación del proceso de migración celular, particularmente en estudios sobre progresión tumoral y metástasis.

## Referencias

[![Poster](https://img.shields.io/badge/Poster-PDF-blue)](https://github.com/L3N73J4M4N/SCIAN-Force/releases/download/SCIAN-Force/2024_Poster_MCarvajal_SCIANForce_SCIANDrop.pdf)

**[1]** R. Lucio, A. Arashiro, and N. R. Demarquette, "Pendant drop method: A generalization of the Laplace equation to measure interfacial tension," in *Materials Research*, vol. 18, no. 1, pp. 23–32, 2015. DOI: [10.1590/1516-1439.267314](https://doi.org/10.1590/1516-1439.267314).

**[2]** O. Càmpas et al., "Quantifying cell-generated forces within living embryonic tissues," *Nature Methods*, vol. 11, no. 2, pp. 183–189, 2014. DOI: [10.1038/nmeth.2761](https://doi.org/10.1038/nmeth.2761).

**[3]** S. Rusinkiewicz, "Estimating Curvatures and Their Derivatives on Triangle Meshes," in *Proceedings of the IEEE Symposium on 3D Data Processing, Visualization, and Transmission*, Thessaloniki, Greece, 2004, pp. 486–493. DOI: [10.1109/TDPVT.2004.1335277](https://doi.org/10.1109/TDPVT.2004.1335277).

**[4]** M. Mayer et al., "Visualization of fluid dynamics using triangular meshes and curvature-driven smoothing," in *Mathematics and Visualization*, Springer, 2003, pp. 191–205. DOI: [10.1007/978-3-642-55566-6](https://doi.org/10.1007/978-3-642-55566-6).
   
## Descripción

SCIAN-Force realiza análisis avanzados, que incluyen:
- Reconstrucción 3D de la superficie de las microgotas a través de la clase `Mesh3D`.
- Cálculo de curvaturas promedio y tensiones anisotrópicas a través de la clase `Forces3D`.
- Visualización interactiva de los resultados a través de la clase `Inter3D`
- Diversas funciones útiles para trabajar con mallas triángulares.

## Características principales

- **Reconstrucción 3D**: Genera modelos tridimensionales a partir de imágenes segmentadas.
- **Cálculo de tensiones**: Determina la presión anisotrópica en cada punto de la superficie.
- **Interfaz interactiva**: Permite ajustar parámetros y exportar resultados de forma intuitiva.
- **Validación experimental**: Probado con imágenes sintéticas y datos experimentales.

## Aplicaciones

- Estudio de la progresión tumoral y metástasis.
- Cuantificación de fuerzas ejercidas por células en entornos tridimensionales.
- Análisis de tensiones interfaciales en fluidos biocompatibles.

## Requisitos del sistema

Para ejecutar el software, asegúrate de tener Python 3.x instalado, junto con las siguientes bibliotecas:

- NumPy
- Meshfix
- Rtree
- Plotly
- SciPy
- Pymeshfix
- Pandas
- Tkinter
- Openpyxl
- Pillow (PIL)

Instala las dependencias usando:
```bash
pip install -r requirements.txt
```

## Uso
1. **Entrada:** Carga un conjunto de imágenes segmentadas (stacks) con `Mesh3D`, indicando el tamaño del voxel. Este paso es optativo, puede obtener la reconstrucción 3D de otra forma.
2. **Análisis:** SCIAN-Force calcula las tensiones anisotrópicas con `Forces3D`.
3. **Visualización:** Explora los resultados interactivos generados, incluyendo mapas de tensiones y curvaturas.
4. **Exportación:** Guarda los resultados en formatos compatibles, como hojas de cálculo.

## Descargar SCIAN-Force

Puedes descargar la última versión de **SCIAN-Force** para Windows como un archivo ejecutable. Haz clic en el siguiente enlace para iniciar la descarga:

[Descargar SCIAN-Force.exe (updated 2024/11/16)](https://github.com/L3N73J4M4N/SCIAN-Force/releases/download/SCIAN-Force/SCIAN-Force.v2024.11.08.exe)

## Licencia
Este proyecto está licenciado bajo la Licencia MIT.
