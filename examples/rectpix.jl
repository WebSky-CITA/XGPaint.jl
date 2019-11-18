using XGPaint
using PyCall
using PyPlot
using WCS
enmap = pyimport("pixell.enmap")

wcs2 = WCSTransform(2;
           cdelt = [1., 1.], #cdelt = [5/60., 5/60.],
           ctype = ["RA---CAR", "DEC--CAR"],
           crpix = [181.0, 91.0],
           crval = [0., 0.])
##
l = collect(1:1000)
cl = l .^ -2.5
shape,wcs = enmap.geometry(shape=(180, 360),
                           res=deg2rad(1),pos=(0,0))
modl = enmap.modlmap(shape, wcs);
cmb = enmap.rand_map(shape, wcs,
   reshape(XGPaint.ellpad(cl),(1,1,length(cl)+1)));
##
imshow(cmb)
gcf()

##

coords = zeros(2,3)
coords = [
   0.0  90.0  180.0;
   0.0  00.0  00.0
]
(world_to_pix(wcs2, coords ))

##
coords = [
   1.0  5.0  180.0;
   1.0  1.0  1.0
]
(pix_to_world(wcs2, coords ))
