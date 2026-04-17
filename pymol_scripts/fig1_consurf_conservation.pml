# Figura 1: Mapa de conservación ConSurf sobre 1DKL
# Uso:  pymol -cq fig1_consurf_conservation.pml
# Produce: fig1_consurf_conservation.png

# PDB con B-factor = grado de conservación (1-9) producido por ConSurf
load ../data/1DKL_ATOMS_section_With_ConSurf.pdb, appa

# 1DKL es homodímero. Solo se muestra la cadena A
remove appa and chain B

hide everything
show cartoon, appa
set cartoon_transparency, 0

# Espectro azul (variable) -> blanco -> rojo (conservado) sobre B-factor
spectrum b, blue_white_red, appa and name CA

# Resaltar ligando si existe
show sticks, resn 6PL+IHP+FIT+INS and appa
color yellow, resn 6PL+IHP+FIT+INS and appa

bg_color white
set ray_opaque_background, 1
orient appa
zoom appa, 2
set ray_trace_mode, 1
set ray_shadows, 0
ray 1600, 1200
png fig1_consurf_conservation.png, dpi=300
