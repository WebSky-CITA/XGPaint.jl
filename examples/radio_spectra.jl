using XGPaint
using Healpix

nside = 4096
m = Map{Float64,RingOrder}(nside)
XGPaint.catalog2map!(m, [1.0], [0.2], [0.3])
