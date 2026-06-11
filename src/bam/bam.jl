# BAM File Format
# ===============

module BAM

using BioGenerics
using GenomicFeatures
using XAM.SAM
import ..XAM: flags, XAMRecord, XAMReader, XAMWriter,
	ismapped, isprimaryalignment, ispositivestrand, isnextmapped #TODO: Deprecate import of flag queries. These were imported to preseve existing API.

import BGZFLib
import BioAlignments
import Indexes
import BioSequences
import BioGenerics: isfilled, header

import GenomicFeatures: eachoverlap

# BGZFLib uses the BufferIO interface (not Base.IO), so Julia's generic
# read(io, T) / write(io, T) for primitives don't dispatch. Add the minimal
# overloads needed by BAM header reading/writing.

Base.read(io::BGZFLib.BGZFReader, nb::Integer) =
    (v = Vector{UInt8}(undef, nb); read!(io, v); v)

function Base.read(io::BGZFLib.BGZFReader, ::Type{T}) where T<:Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}
    ref = Ref{T}()
    GC.@preserve ref unsafe_read(io, Ptr{UInt8}(Base.unsafe_convert(Ptr{T}, ref)), UInt(sizeof(T)))
    return ltoh(ref[])
end

function Base.write(io::BGZFLib.BGZFWriter, x::T) where T<:Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}
    v = htol(x)
    ref = Ref(v)
    GC.@preserve ref return unsafe_write(io, Ptr{UInt8}(Base.unsafe_convert(Ptr{T}, ref)), UInt(sizeof(T)))
end

Base.write(io::BGZFLib.BGZFWriter, v::AbstractVector{UInt8}) =
    unsafe_write(io, pointer(v), UInt(length(v)))

Base.write(io::BGZFLib.BGZFWriter, s::AbstractString) =
    unsafe_write(io, pointer(s), UInt(ncodeunits(s)))

Base.write(io::BGZFLib.BGZFWriter, c::Char) = write(io, UInt8(c))


include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
