var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = XGPaint","category":"page"},{"location":"#XGPaint","page":"Home","title":"XGPaint","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [XGPaint]","category":"page"},{"location":"#XGPaint.AbstractForegroundModel","page":"Home","title":"XGPaint.AbstractForegroundModel","text":"All foreground models inherit from this type.\n\n\n\n\n\n","category":"type"},{"location":"#XGPaint.CIB_Planck2013","page":"Home","title":"XGPaint.CIB_Planck2013","text":"CIB_Planck2013{T}(model parameters...)\n\nDefine CIB model parameters. Defaults are from Viero et al. 2013.\n\nmodel = CIB_Planck2013{Float32}(shang_Mpeak=10^12.4)\n\n\n\n\n\n","category":"type"},{"location":"#XGPaint.Enmap","page":"Home","title":"XGPaint.Enmap","text":"Map type, contains an AbstractArray and a WCS object, but behaves like the AbstractArray it contains for array operations.\n\nIt only implements the subset of Base.Array operations which are common on maps. You should work with the data directly using enmap_instance.data if you need additional Array functions.\n\n\n\n\n\n","category":"type"},{"location":"#XGPaint.Radio_Sehgal2009","page":"Home","title":"XGPaint.Radio_Sehgal2009","text":"Radio_Sehgal2009{T}(model parameters...)\n\nDefine CIB model parameters. Defaults are from Viero et al. 2013.\n\nmodel = CIBModel{Float32}(a_0=0.4)\n\n\n\n\n\n","category":"type"},{"location":"#XGPaint.build_c_lnm2r_interpolator-Union{Tuple{}, Tuple{T}} where T","page":"Home","title":"XGPaint.build_c_lnm2r_interpolator","text":"Generates an interpolator r(c, lnm)\n\nGenerate a LinearInterpolation object that turns concentration and ln(M_halo) into satellite radius.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.build_r2z_interpolator-Union{Tuple{T}, Tuple{T, T, Cosmology.AbstractCosmology}} where T","page":"Home","title":"XGPaint.build_r2z_interpolator","text":"Construct a fast r2z linear interpolator.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.build_shang_interpolator-Union{Tuple{T}, Tuple{T, T, XGPaint.AbstractCIBModel}} where T","page":"Home","title":"XGPaint.build_shang_interpolator","text":"Build a linear interpolation function which maps log(Mh) to Nsat.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.build_sigma_sat_ln_interpolator-Union{Tuple{T}, Tuple{T, XGPaint.AbstractCIBModel}} where T","page":"Home","title":"XGPaint.build_sigma_sat_ln_interpolator","text":"Build a linear interpolator that takes in ln(M_halo) and returns sigma.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.chunk-Tuple{Any, Integer}","page":"Home","title":"XGPaint.chunk","text":"Generates a list of tuples which describe starting and ending chunk indices. Useful for parallelizing an array operation.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.ellpad-Union{Tuple{Array{T, N}}, Tuple{N}, Tuple{T}} where {T, N}","page":"Home","title":"XGPaint.ellpad","text":"Utility function which prepends some zeros to an array. It makes a copy instead of modifying the input.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.enfft!-Union{Tuple{XGPaint.Enmap{T, N, AA}}, Tuple{AA}, Tuple{N}, Tuple{T}} where {T, N, AA}","page":"Home","title":"XGPaint.enfft!","text":"Physically normalized enmap FFT\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.enifft!-Union{Tuple{XGPaint.Enmap{T, N, AA}}, Tuple{AA}, Tuple{N}, Tuple{T}} where {T, N, AA}","page":"Home","title":"XGPaint.enifft!","text":"Physically normalized enmap inverse FFT\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.generate_sources-Union{Tuple{TH}, Tuple{T}, Tuple{XGPaint.AbstractCIBModel{T}, Cosmology.FlatLCDM{T}, AbstractMatrix{TH}, AbstractVector{TH}}} where {T, TH}","page":"Home","title":"XGPaint.generate_sources","text":"Produce a source catalog from a model and halo catalog.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.generate_sources-Union{Tuple{TH}, Tuple{T}, Tuple{XGPaint.AbstractRadioModel{T}, Cosmology.FlatLCDM{T}, AbstractMatrix{TH}, AbstractVector{TH}}} where {T, TH}","page":"Home","title":"XGPaint.generate_sources","text":"Produce a source catalog from a model and halo catalog.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.generate_subhalo_offsets-Tuple{Any}","page":"Home","title":"XGPaint.generate_subhalo_offsets","text":"Generate an array where the value at index i corresponds to the index of the first source of halo i. Takes an array where the value at index i corresponds to the number of subhalos that halo i has.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.get_angles-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"Home","title":"XGPaint.get_angles","text":"Compute angles of halos\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.get_basic_halo_properties-Union{Tuple{T}, Tuple{Matrix{T}, XGPaint.AbstractForegroundModel, Cosmology.FlatLCDM{T}, Healpix.Resolution}} where T","page":"Home","title":"XGPaint.get_basic_halo_properties","text":"Fill in basic halo properties.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.get_cosmology-Union{Tuple{Type{T}}, Tuple{T}} where T","page":"Home","title":"XGPaint.get_cosmology","text":"Construct a background cosmology.\n\nThis function duplicates the cosmology() function in Cosmology.jl, but with typing. The type of the cosmology will the type of h and OmegaM. This is primarily for keeping the code entirely in Float32 or Float64.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.get_interpolators-Union{Tuple{T}, Tuple{XGPaint.AbstractCIBModel, Cosmology.FlatLCDM{T}, T, T}} where T","page":"Home","title":"XGPaint.get_interpolators","text":"Construct the necessary interpolator set.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.hod_sehgal-Union{Tuple{T}, Tuple{Any, Any, Radio_Sehgal2009{T}}} where T","page":"Home","title":"XGPaint.hod_sehgal","text":"Populate halos with radio sources according to the HOD in Sehgal et al. 2009.\n\nThe optional rng parameter provides an array of random number generators, one for each thread.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.integrand_L-Tuple{Any, Any, XGPaint.AbstractCIBModel}","page":"Home","title":"XGPaint.integrand_L","text":"<L_sat> interpolation values\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.l2f-Union{Tuple{T}, Tuple{T, T, T}} where T","page":"Home","title":"XGPaint.l2f","text":"Inverse square law with redshift dependence.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.m2r-Union{Tuple{T}, Tuple{T, Cosmology.FlatLCDM{T}}} where T","page":"Home","title":"XGPaint.m2r","text":"Convert virial mass to virial radius.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.mz2c-Union{Tuple{T}, Tuple{T, T, Cosmology.FlatLCDM{T}}} where T","page":"Home","title":"XGPaint.mz2c","text":"Compute concentration factor from Duffy et al. 2008.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.paint!-Union{Tuple{T}, Tuple{Healpix.Map{T, Healpix.RingOrder, AA} where AA<:AbstractVector{T}, T, XGPaint.AbstractCIBModel, Any}} where T","page":"Home","title":"XGPaint.paint!","text":"Paint a source catalog onto a map.\n\nThis function creates the arrays for you.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.paint!-Union{Tuple{T}, Tuple{T_map}, Tuple{Healpix.Map{T_map, Healpix.RingOrder, AA} where AA<:AbstractVector{T_map}, T, XGPaint.AbstractCIBModel, Any, AbstractArray, AbstractArray}} where {T_map, T}","page":"Home","title":"XGPaint.paint!","text":"Paint a source catalog onto a map, recording the fluxes.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.process_centrals!-Union{Tuple{T}, Tuple{XGPaint.AbstractCIBModel{T}, Cosmology.FlatLCDM{T}, Healpix.Resolution}} where T","page":"Home","title":"XGPaint.process_centrals!","text":"Fill up arrays with information related to CIB central sources.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.process_sats!-Union{Tuple{T}, Tuple{XGPaint.AbstractCIBModel{T}, Cosmology.FlatLCDM{T}, Healpix.Resolution}} where T","page":"Home","title":"XGPaint.process_sats!","text":"Fill up arrays with information related to CIB satellites.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.read_halo_catalog_hdf5-Tuple{Any}","page":"Home","title":"XGPaint.read_halo_catalog_hdf5","text":"Utility function to read an HDF5 table with x, y, z, M_h as the four rows. The hdf5 record is \"halos\".\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.sehgal_LF!-Union{Tuple{T}, Tuple{Vector{T}, Any, Any, T, Any, Vector{T}}} where T","page":"Home","title":"XGPaint.sehgal_LF!","text":"Fills the result array with draws from the luminosity function.\n\n\n\n\n\n","category":"method"},{"location":"#XGPaint.shang_z_evo-Union{Tuple{T}, Tuple{T, XGPaint.AbstractCIBModel}} where T","page":"Home","title":"XGPaint.shang_z_evo","text":"Compute redshift evolution factor for LF.\n\n\n\n\n\n","category":"method"}]
}
