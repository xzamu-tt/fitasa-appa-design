# Figura 4: Superposición WT 1DKL vs predicciones AF3 de las 3 variantes finales
# Uso:  pymol -cq fig4_wt_vs_variants.pml
# Produce: fig4_wt_vs_variants.png
#
# Requiere los .cif de AF3 para cada variante (rutas locales del lab):
#   Fitasa-01: t03_2025_07_03_11_59_3/fold_2025_07_03_11_59_model_0.cif
#   Fitasa-02: t02_2025_07_03_11_48/fold_2025_07_03_11_48_model_0.cif
#   Fitasa-03: t02_2025_07_03_11_50_3/fold_2025_07_03_11_50_model_0.cif

fetch 1dkl, async=0, name=WT

# 1DKL es homodímero; solo se muestra la cadena A (los modelos AF3 ya son monómeros)
remove chain B, WT

# Comentar/descomentar según la disponibilidad de los .cif en el equipo local:
# load path/to/t03_2025_07_03_11_59_3/fold_2025_07_03_11_59_model_0.cif, Fitasa01
# load path/to/t02_2025_07_03_11_48/fold_2025_07_03_11_48_model_0.cif,   Fitasa02
# load path/to/t02_2025_07_03_11_50_3/fold_2025_07_03_11_50_model_0.cif, Fitasa03

align Fitasa01, WT
align Fitasa02, WT
align Fitasa03, WT

hide everything
show cartoon

color grey70, WT
color salmon, Fitasa01
color skyblue, Fitasa02
color forest, Fitasa03

set cartoon_transparency, 0.3, WT

bg_color white
set ray_opaque_background, 1
orient WT
zoom WT, 3
set ray_trace_mode, 1
set ray_shadows, 0
ray 1600, 1200
png fig4_wt_vs_variants.png, dpi=300
