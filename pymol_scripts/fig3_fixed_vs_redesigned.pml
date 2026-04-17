# Figura 3: Residuos fijados (rojo) vs rediseñados (azul)
# Uso:  pymol -cq fig3_fixed_vs_redesigned.pml
# Produce: fig3_fixed_vs_redesigned.png

fetch 1dkl, async=0

# 1DKL es homodímero; solo se muestra la cadena A
remove chain B, 1dkl

hide everything
show cartoon, 1dkl
color grey80, 1dkl

# 201 posiciones REDISEÑADAS (ver data/FIXED_residues_format.txt)
select redesigned, chain A and resi 1+2+4+8+9+11+12+13+14+15+16+17+18+19+20+21+22+23+24+33+34+37+38+40+41+42+43+45+47+48+49+50+51+54+58+59+61+63+66+70+71+72+73+76+77+78+81+83+87+88+90+91+92+93+94+96+97+99+101+103+104+105+106+107+108+118+122+123+124+125+126+127+128+129+130+131+133+136+140+142+143+149+150+155+159+163+166+169+170+175+178+179+188+190+195+198+204+205+207+209+211+212+213+214+215+216+217+219+220+221+222+223+224+225+226+227+228+230+232+233+234+235+236+243+246+247+249+250+251+252+253+254+255+256+257+258+259+260+261+262+263+264+265+267+270+271+272+275+278+279+290+297+300+302+303+304+305+306+308+309+310+312+313+314+315+316+318+320+321+322+323+324+325+327+328+329+330+331+332+333+334+335+336+337+338+344+348+349+350+354+355+356+357+359+360+361+365+366+368+371+372+377+380+381+382+389+391+393+396+400+403+406+407+408+410

# Fijadas = complemento
select fixed, chain A and not redesigned and polymer

color skyblue, redesigned
color firebrick, fixed

show spheres, fixed and name CA
show spheres, redesigned and name CA
set sphere_scale, 0.6, name CA

bg_color white
set ray_opaque_background, 1
orient 1dkl
zoom 1dkl, 3
set ray_trace_mode, 1
set ray_shadows, 0
ray 1600, 1200
png fig3_fixed_vs_redesigned.png, dpi=300
