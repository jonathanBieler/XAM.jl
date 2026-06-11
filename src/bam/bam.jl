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

# BGZFLib exposes BufferIO's `AbstractBufReader`/`AbstractBufWriter` interface
# rather than `Base.IO`, and that interface intentionally omits `read(io, T)`
# for multi-byte primitives. We provide the little-endian primitive/byte IO that
# BAM header and record parsing needs as *local* helpers: defining methods on
# `Base.read`/`Base.write` for the BGZFLib-owned reader/writer types would be
# type piracy (neither the function nor the argument types are owned here).

# Read exactly `nb` bytes. BufferIO's `read!` throws on a truncated stream,
# unlike `read(io, nb)`, which would silently return a short vector.
function bam_read(io::BGZFLib.BGZFReader, nb::Integer)
    data = Vector{UInt8}(undef, nb)
    read!(io, data)
    return data
end

bam_read(io::BGZFLib.BGZFReader, ::Type{UInt8}) = read(io, UInt8)

function bam_read(io::BGZFLib.BGZFReader, ::Type{T}) where {T<:Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}}
    ref = Ref{T}()
    # BufferIO's `unsafe_read` returns a short count at EOF instead of throwing.
    n = GC.@preserve ref unsafe_read(io, Ptr{UInt8}(Base.unsafe_convert(Ptr{T}, ref)), UInt(sizeof(T)))
    n < sizeof(T) && throw(EOFError())
    return ltoh(ref[])
end

function bam_write(io::BGZFLib.BGZFWriter, x::T) where {T<:Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}}
    v = htol(x)
    ref = Ref(v)
    return GC.@preserve ref unsafe_write(io, Ptr{UInt8}(Base.unsafe_convert(Ptr{T}, ref)), UInt(sizeof(T)))
end

bam_write(io::BGZFLib.BGZFWriter, data::AbstractVector{UInt8}) =
    unsafe_write(io, pointer(data), UInt(length(data)))

bam_write(io::BGZFLib.BGZFWriter, s::AbstractString) =
    unsafe_write(io, pointer(s), UInt(ncodeunits(s)))

bam_write(io::BGZFLib.BGZFWriter, x::UInt8) = write(io, x)

bam_write(io::BGZFLib.BGZFWriter, c::Char) = bam_write(io, UInt8(c))


include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
