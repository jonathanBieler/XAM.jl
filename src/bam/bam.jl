# BAM File Format
# ===============

module BAM

using BioGenerics
using GenomicFeatures
using XAM.SAM
import ..XAM: flags, XAMRecord, XAMReader, XAMWriter,
	ismapped, isprimaryalignment, ispositivestrand, isnextmapped #TODO: Deprecate import of flag queries. These were imported to preseve existing API.

import BGZFLib
import MemoryViews
import BioAlignments
import Indexes
import BioSequences
import BioGenerics: isfilled, header

import GenomicFeatures: eachoverlap

# Indexes.Chunk.start/stop are BGZFStreams.VirtualOffset (a 64-bit primitive type).
# Convert to BGZFLib.VirtualOffset without importing BGZFStreams.
@inline function _to_virtual_offset(vo)
    u = reinterpret(UInt64, vo)
    BGZFLib.VirtualOffset(u >> 16, u & 0xffff)
end

# BGZFLib uses the BufferIO interface (not Base.IO), so Base.read(io, ::Type{T})
# for non-UInt8 primitives is not provided. These private helpers own both the
# function name and the dispatch, avoiding type piracy.

@inline function _bam_read(io::BGZFLib.BGZFReader, ::Type{T}) where T <: Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}
    ref = Ref{T}()
    GC.@preserve ref BGZFLib.BufferIO.read_all!(io, MemoryViews.MemoryView(unsafe_wrap(Array, Ptr{UInt8}(Base.unsafe_convert(Ptr{T}, ref)), sizeof(T))))
    return ltoh(ref[])
end

@inline function _bam_write(io, x::T) where T <: Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}
    v = htol(x)
    ref = Ref(v)
    GC.@preserve ref return unsafe_write(io, Ptr{UInt8}(Base.unsafe_convert(Ptr{T}, ref)), UInt(sizeof(T)))
end

include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
