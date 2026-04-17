# Figura 2: Residuos del sitio activo de AppA (7 Å del ligando)
# Uso:  pymol -cq fig2_active_site.pml
# Produce: fig2_active_site.png

fetch 1dkl, async=0

# 1DKL es homodímero; solo se muestra la cadena A
remove chain B, 1dkl

hide everything
show cartoon, 1dkl
color grey70, 1dkl

# 39 residuos del sitio activo
select active_site, chain A and resi 16+17+18+20+21+22+23+24+45+88+90+92+93+125+126+129+204+212+213+215+216+217+219+220+223+250+253+254+258+267+302+303+304+305+306+325+327+328+329

show sticks, active_site
color salmon, active_site

# Motivos canónicos HAP: RHGXRXP (R16-H17-G18-R20) y HD (H303-D304)
select hap_motif_N, chain A and resi 16+17+18+20
select hap_motif_C, chain A and resi 303+304
color red, hap_motif_N
color firebrick, hap_motif_C
show sticks, hap_motif_N hap_motif_C

label hap_motif_N and name CA, "%s%s" % (resn, resi)
label hap_motif_C and name CA, "%s%s" % (resn, resi)
set label_size, 18

bg_color white
set ray_opaque_background, 1
orient active_site
zoom active_site, 5
set ray_trace_mode, 1
set ray_shadows, 0
ray 1600, 1200
png fig2_active_site.png, dpi=300
