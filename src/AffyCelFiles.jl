module AffyCelFiles

#http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html

mutable struct Cel

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
    firstPos=ntoh(read(io,UInt32))

    #read generic data header
    data_type_identifier=read_string(io)

    


end

function read_string(io::IO)
    length=ntoh(read(io,Int32))
    c=read(io,length)
    return join(Array{Char}(c))
end

end # module
