module Series

import DICOM
using Images
using ImageView
using ImageCore

export load_dicom_volume, apply_window, window_preset, view_slice

struct Volume{T, A<:AbstractArray{T,3}}
    data::A
    
end

"""
    apply_window(img, WW, WL) -> Array{Float32,2}

Apply CT window width (WW) and level (WL) to a 2D slice and map to [0, 1].
"""
function apply_window(img::AbstractArray{<:Real,2}, WW::Real, WL::Real)
    WW == 0 && error("WW (window width) cannot be 0.")
    lower = WL - WW/2
    upper = WL + WW/2
    imgc  = clamp.(float.(img), lower, upper)
    return Float32.((imgc .- lower) ./ WW)
end

"""
    window_preset(name::Symbol) -> (WW, WL)

Common presets: :lung, :brain, :bone.
"""
function window_preset(name::Symbol)
    name === :lung  && return (1500.0, -600.0)
    name === :brain && return (  80.0,    40.0)
    name === :bone  && return (2000.0,   500.0)
    error("Unknown preset: $name. Try :lung, :brain, or :bone.")
end

# -- Internal: convert one DICOM slice to HU -----------------------------------

function _slice_to_hu(ds)
    # Try the built-in decoder first (handles compressed PixelData if codecs available)
    img = try
        DICOM.pixeldata(ds)
    catch
        # Manual fallback for uncompressed data
        rows   = Int(ds[DICOM.tag"Rows"])
        cols   = Int(ds[DICOM.tag"Columns"])
        bits   = Int(get(ds, DICOM.tag"BitsAllocated", 16))
        signed = get(ds, DICOM.tag"PixelRepresentation", 0) == 1  # 0=unsigned, 1=signed
        T = bits == 16 ? (signed ? Int16 : UInt16) : UInt8

        raw = ds[DICOM.tag"PixelData"]
        reshape(reinterpret(T, raw), (cols, rows))'   # (rows, cols)
    end

    slope     = float(get(ds, DICOM.tag"RescaleSlope",     1.0))
    intercept = float(get(ds, DICOM.tag"RescaleIntercept", 0.0))
    return Float32.(slope .* float.(img) .+ intercept)     # Hounsfield units
end

"""
    load_dicom_volume(folder::AbstractString) -> Array{Float32,3}

Read a DICOM series from `folder`, sort slices, convert to HU, and stack into a 3D volume.
"""
function load_dicom_volume(folder::AbstractString)
    slices = DICOM.dcmdir_parse(folder)
    isempty(slices) && error("No DICOM files found in: $folder")

    # Prefer sorting by InstanceNumber; else fall back to ImagePositionPatient (z)
    has_inst = all(!isnothing(get(s, DICOM.tag"InstanceNumber", nothing)) for s in slices)
    if has_inst
        slices = sort(slices, by = s -> get(s, DICOM.tag"InstanceNumber", 0))
    else
        zpos(s) = try
            pos = get(s, DICOM.tag"ImagePositionPatient", nothing)
            pos === nothing ? 0.0 : float(pos[3])
        catch
            0.0
        end
        slices = sort(slices, by = zpos)
    end

    vols = Vector{Array{Float32,2}}(undef, length(slices))
    @inbounds for i in eachindex(slices)
        vols[i] = _slice_to_hu(slices[i])
    end

    # Sanity check: all slices same size
    sz = size(vols[1])
    @assert all(size(v) == sz for v in vols) "Slice sizes differ across the series."

    return cat(vols...; dims=3)   # Float32 volume (rows, cols, slices)
end

"""
    view_slice(img2d)

Render a 2D slice. Safely maps to Gray{N0f8} before calling `imshow`.
"""
function view_slice(img2d::AbstractArray{<:Real,2})
    imshow(Gray.(N0f8.(clamp01.(float.(img2d)))))
    return nothing
end

end