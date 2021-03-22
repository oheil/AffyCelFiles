module AffyCelFiles

#http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html


struct DataHeader
    isDefined::Bool
    dataTypeIdentifier::String
    uniqueID::String
    creationDateTime::String
    locale::String
    nParameters::Int
    names::Array{String,1}
    rawvalues::Array{String,1}
    values::Array{String,1}
    types::Array{String,1}
    nParentFileHeaders::Int
    parentFileHeaders::Array{DataHeader,1}
    function DataHeader(isDefined=false)
        new(
            isDefined
        )
    end
    function DataHeader(args...)
        new(
            args...
        )
    end
end

struct DataGroup
    isDefined::Bool
    test::Int
    function DataGroup(isDefined=false)
        new(
            isDefined
        )
    end
    function DataGroup(args...)
        new(
            args...
        )
    end
end

struct Cel
    isDefined::Bool
    nGroups::Int
    dataHeader::DataHeader
    dataGroups::Array{DataGroup,1}
    function Cel(isDefined=false)
		new(
			isDefined
		)
	end
    function Cel(args...)
		new(
			args...
		)
	end
end

function cel_read(file::AbstractString)
	io=open(file,"r")
	cel=cel_read(io)
	close(io)
	cel
end

function cel_read(io::IO)::Cel
	seekstart(io)
    magic=read(io,Int8)
    if magic == 59
        #Command Console generic data file format
        #  http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/generic.html
        return cel_read_generic(io)
    end
    
	seekstart(io)
    magic=read(io,Int32)
    if magic == 64
        #Version 4 Format
        #  http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html
        return cel_read_v4(io)
    end

	@warn "Unknown CEL file version"
	return Cel()
end

function cel_read_v4(io::IO)
    #not implemented, need example file
	return Cel()
end

function cel_read_generic(io::IO)
    version=read(io,Int8)
    if version != 1
        @warn "Version number of Command Console generic data file is "*string(version)*". Expected 1."
    end
    nGroups=ntoh(read(io,Int32))
    firstGroupPos=ntoh(read(io,UInt32))

    dataHeader=cel_read_generic_header(io)
    #println(firstGroupPos)
    #println(position(io))

    if firstGroupPos != position(io)
        @warn "current io position is "*string(position(io))*", should be "*string(firstGroupPos)*", performing seek to "*string(firstGroupPos)
        seek(io,firstGroupPos)
    end

    dataGroups=Array{DataGroup}(undef,nGroups)
    nextGroupPos = firstGroupPos
    while nextGroupPos > 0
        nextGroupPos=ntoh(read(io,UInt32))

    end

    return Cel(true,nGroups,dataHeader,dataGroups)
end

function cel_read_generic_header(io::IO)
    #read generic data header
    data_type_identifier=read_string(io)
    guid=read_string(io)
    creation=read_string(io)
    locale=read_wstring(io)
    nParameters=ntoh(read(io,Int32))
    names=Array{String}(undef,nParameters)
    rawvalues=Array{String}(undef,nParameters)
    values=Array{String}(undef,nParameters)
    types=Array{String}(undef,nParameters)
    for i in 1:nParameters
        #println(position(io))
        names[i]=read_wstring(io)
        rawvalues[i]=read_value(io)
        types[i]=read_wstring(io)
        values[i]=calvin_convert(rawvalues[i], Val(Symbol(types[i])) )
    end
    nParentFileHeaders=ntoh(read(io,Int32))
    parentFileHeaders=Array{DataHeader}(undef,nParentFileHeaders)
    if nParentFileHeaders > 0
        for i in 1:nParentFileHeaders
            parentFileHeaders[i]=cel_read_generic_header(io)
        end
    end

    return DataHeader(
        true,
        data_type_identifier,
        guid,
        creation,
        locale,
        nParameters,
        names,
        rawvalues,
        values,
        types,
        nParentFileHeaders,
        parentFileHeaders
    )
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-integer-8")})
    return string(reinterpret(Int8,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-unsigned-integer-8")})
    return string(reinterpret(UInt8,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-integer-16")})
    return string(reinterpret(Int16,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-unsigned-integer-16")})
    return string(reinterpret(UInt16,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-integer-32")})
    return string(reinterpret(Int32,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-unsigned-integer-32")})
    return string(reinterpret(UInt32,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/x-calvin-float")})
    return string(reinterpret(Float32,UInt16.(collect(v)[2:-1:1]))[1])
end

function calvin_convert(v::String,::Val{Symbol("text/plain")})
    chars=collect(v)
    eos=findfirst(c->c=='\0',chars)
    return isnothing(eos) ? v : v[1:(eos-1)] 
end

function calvin_convert(v::String,::Val{Symbol("text/ascii")})
    chars=Char.(reinterpret(UInt8,UInt16.(collect(v))))
    eos=findfirst(c->c=='\0',chars)
    return isnothing(eos) ? join(chars) : join(chars[1:(eos-1)])
end

function read_string(io::IO)
    length=ntoh(read(io,Int32))
    c=read(io,length)
    return join(Array{Char}(c))
end

function read_wstring(io::IO)
    length=2*ntoh(read(io,Int32))
    b=zeros(UInt8,length)
    readbytes!(io,b,length)
    b=ntoh.(reinterpret(UInt16,b))
    return transcode(String,b)
end

function read_value(io::IO)
    length=ntoh(read(io,Int32))
    bufferLength=length
    if !iseven(length)
        bufferLength+=1
    end
    b=zeros(UInt8,bufferLength)
    readbytes!(io,b,length)
    b=ntoh.(reinterpret(UInt16,b))
    return transcode(String,b)
end

end # module
