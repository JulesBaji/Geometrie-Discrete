# Géométrie discrète : analyse de grains de riz

Projet de géométrie discrète du master Informatique ISICG (Université de Limoges, 2024-2025), en C++ avec la bibliothèque [DGtal](https://dgtal.org/).

À partir d'images binaires de grains de riz de trois variétés (camargue, japonais, basmati), le programme mesure la forme de chaque grain et compare les variétés.

![Grains de riz basmati, camargue et japonais : contours puis polygones, colorés selon la circularité](docs/geometrie.gif)

*Basmati, camargue et japonais (de gauche à droite). Les grains qui touchent le bord sont écartés (en rouge), puis le contour de chaque grain est suivi et polygonisé ; la couleur finale indique la circularité, de l'orange (allongé) au bleu (rond). Contours et polygones sont ceux calculés par le programme ; l'animation est générée à partir de ces résultats.*

## Étapes

1. Lecture des images et extraction des **composantes connexes** (adjacence 8-4).
2. Suppression des grains qui touchent le bord de l'image.
3. **Suivi de contour** dans l'espace de Khalimsky.
4. **Polygonisation** du contour par segmentation en segments de droite discrets (code de Freeman, segmentation gloutonne).
5. Mesure de l'**aire**, du **périmètre** (discrets et polygonaux) et de la **circularité**.
6. Moyennes et variances par variété.

Chaque grain est exporté en PDF (composante, contour, polygone) dans des dossiers par image. Le rapport est dans `Geometrie_Discrete.pdf`.

## Compilation

Nécessite DGtal et CMake.

```bash
mkdir build-local && cd build-local
cmake ..
make
./main
```

## Technologies

C++, DGtal, CMake.
