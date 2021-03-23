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

struct DataSet
    isDefined::Bool
    dataSetName::String
    nParameters::Int
    parameterNames::Array{String,1}
    parameterRawvalues::Array{String,1}
    parameterValues::Array{String,1}
    parameterTypes::Array{String,1}
    dataSetColNames::Array{String,1}
    dataSetColTypeCodes::Array{Int8,1}
    dataSetColTypeSizes::Array{Int32,1}
    dataSetColTypes::Array{Type,1}
    dataSetColumns::Dict{Int,AbstractArray}
    function DataSet(isDefined=false)
        new(
            isDefined
        )
    end
    function DataSet(args...)
        new(
            args...
        )
    end
end

struct DataGroup
    isDefined::Bool
    nDataSets::Int
    dataSets::Array{DataSet,1}
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
        @warn "current io position is "*string(position(io))*", should be first data group at "*string(firstGroupPos)*", performing seek to "*string(firstGroupPos)
        seek(io,firstGroupPos)
    end

    dataGroups=Array{DataGroup}(undef,nGroups)
    dataGroupsIndex=1
    nextGroupPos=firstGroupPos
    while nextGroupPos > 0 && dataGroupsIndex <= nGroups
        nextGroupPos=ntoh(read(io,UInt32))
        nextDataSet=ntoh(read(io,UInt32))
        nDataSets=ntoh(read(io,Int32))
        groupName=read_wstring(io)
        dataSets=Array{DataSet}(undef,nDataSets)
        dataSetsIndex=1
        while dataSetsIndex <= nDataSets
            if nextDataSet != position(io)
                #@warn "current io position is "*string(position(io))*", should be next dataset at "*string(nextDataSet)*", performing seek to "*string(nextDataSet)
                seek(io,nextDataSet)
            end
            dataElementPos=ntoh(read(io,UInt32))
            nextDataSet=ntoh(read(io,UInt32))
            dataSetName=read_wstring(io)
            nParameters=ntoh(read(io,Int32))
            names=Array{String}(undef,nParameters)
            rawvalues=Array{String}(undef,nParameters)
            values=Array{String}(undef,nParameters)
            types=Array{String}(undef,nParameters)
            for i in 1:nParameters
                names[i]=read_wstring(io)
                rawvalues[i]=read_value(io)
                types[i]=read_wstring(io)
                values[i]=calvin_convert(rawvalues[i], Val(Symbol(types[i])) )
            end
            nColumns=ntoh(read(io,UInt32))
            colNames=Array{String}(undef,nColumns)
            colTypeCodes=Array{Int8}(undef,nColumns)
            colTypeSizes=Array{Int32}(undef,nColumns)
            colTypes=Array{Type}(undef,nColumns)
            rowSize=0
            for i in 1:nColumns
                colNames[i]=read_wstring(io)
                colTypeCodes[i]=read(io,Int8)
                colSize=ntoh(read(io,Int32))
                colTypeSizes[i]=colSize
                rowSize+=colSize
                colTypes[i]=decode_type(colTypeCodes[i],colSize)
            end
            nRows=ntoh(read(io,UInt32))
            if dataElementPos != position(io)
                @warn "current io position is "*string(position(io))*", should be next data element at "*string(dataElementPos)*", performing seek to "*string(dataElementPos)
                seek(io,dataElementPos)
            end
            dataSetColumns=Dict{Int,AbstractArray}()
            buffer=zeros(UInt8,nRows*rowSize)
            readbytes!(io,buffer,nRows*rowSize)
            #permBuffer=PermutedDimsArray(reshape(buffer,(rowSize,nRows)),(2,1))
            #startColIndex=0
            #endColIndex=0
            #for i in 1:nColumns
            #    startColIndex=endColIndex+1
            #    endColIndex=startColIndex+colTypeSizes[i]-1
            #    if colTypeCodes[i] <= 6
            #        column=map( row -> reinterpret(colTypes[i],reverse(row[startColIndex:endColIndex]))[1], eachrow(permBuffer))
            #    elseif colTypeCodes[i] == 7 || colTypeCodes[i] == 8
            #        @warn "type string not implemented"
            #        column=Array{String}(undef,0)
            #    else
            #        @error "unknown type code, code="*string(colTypeCodes[i])
            #        column=Array{String}(undef,0)
            #    end
            #    dataSetColumns[i]=column
            #end
            for i in 1:nColumns
                if colTypeCodes[i] <= 6
                    dataSetColumns[i]=Array{colTypes[i]}(undef,nRows)
                elseif colTypeCodes[i] == 7 || colTypeCodes[i] == 8
                    @warn "type string not implemented"
                    column=Array{String}(undef,0)
                else
                    @error "unknown type code, code="*string(colTypeCodes[i])
                    column=Array{String}(undef,0)
                end
            end
            for j in 1:nRows
                startColIndex=(j-1)*rowSize
                endColIndex=startColIndex
                for i in 1:nColumns
                    startColIndex=endColIndex+1
                    endColIndex=startColIndex+colTypeSizes[i]-1
                    if colTypeCodes[i] <= 6
                        dataSetColumns[i][j]=reinterpret(colTypes[i],reverse(buffer[startColIndex:endColIndex]))[1]
                    end
                end   
            end
            dataSets[dataSetsIndex]=DataSet(true,dataSetName,nParameters,names,rawvalues,values,types,colNames,colTypeCodes,colTypeSizes,colTypes,dataSetColumns)
            dataSetsIndex+=1
        end
        dataGroups[dataGroupsIndex]=DataGroup(true,nDataSets,dataSets)
        dataGroupsIndex+=1
    end
    return Cel(true,nGroups,dataHeader,dataGroups)
end

function decode_type(c::Int8,s::Int32)
    if c==0 && s==1
        return Int8
    elseif  c==1 && s==1
        return UInt8
    elseif  c==2 && s==2
        return Int16
    elseif  c==3 && s==2
        return UInt16
    elseif  c==4 && s==4
        return Int32
    elseif  c==5 && s==4
        return UInt32
    elseif  c==6 && s==2
        return Float16
    elseif  c==6 && s==4
        return Float32
    elseif  c==6 && s==8
        return Float64
    elseif  c==7 || c==8
        return String
    else
        @error "unknown type code and size combination, code="*string(c)*" and size="*string(s)*", falling back to String"
    end
    return String
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
