module AffyCelFiles

#https://media.affymetrix.com/support/developer/powertools/changelog/file-formats.html

using CRC32c

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

#https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html

struct Cel
    isDefined::Bool
    nGroups::Int
    dataHeader::DataHeader
    dataGroups::Array{DataGroup,1}
    function Cel(isDefined=false)
		new(
			isDefined,
            -1
		)
	end
    function Cel(args...)
		new(
			args...
		)
	end
end

struct Qc
    qc::Dict{String,String}
    cells::Vector{Vector{Int}}  #only INDEX
    function Qc()
        new(
            Dict{String,String}(),
            Vector{Vector{Int}}(undef,1)
        )
    end    
end

struct Block
    block::Dict{String,String}
    cells::Vector{Vector{Int}}  #only INDEX
    pm::Vector{Vector{Bool}}    #PBASE!=TBASE => perfect match
    function Block()
        new(
            Dict{String,String}(),
            Vector{Vector{Int}}(undef,1),
            Vector{Vector{Bool}}(undef,1)
        )
    end    
end

struct Unit
    unit::Dict{String,String}
    blocks::Vector{Vector{Block}}
    function Unit()
        new(
            Dict{String,String}(),
            Vector{Vector{Block}}(undef,1)
        )
    end
end

#https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html

struct Cdf
    cdf::Dict{String,String}
    chip::Dict{String,String}
    qc::Vector{Qc}
    units::Vector{Unit}
    definedUnitIndices::Dict{Int,Bool}
    function Cdf()
        new(
            Dict{String,String}(),
            Dict{String,String}(),
            Vector{Qc}(),
            Vector{Unit}(),
            Dict{Int,Bool}()
            )
    end
    function Cdf(args...)
		new(
			args...
		)
	end
end

cdfs=Dict{String,AffyCelFiles.Cdf}()
cdfs_files=Dict{String,String}()
cdfs_checksums=Dict{String,UInt32}()

#https://media.affymetrix.com/support/developer/powertools/changelog/file-format-pgf.html

struct Pgf
    chip_name::String
    header_tags::Vector{String}
    header_values::Vector{String}
    probe_ids::Dict{String,Vector{String}}       #probeset_id -> [probe_id1, probe_id2, ...]
    pm::Dict{String,Bool}               #probe_id -> is perfect match(true)
    function Pgf()
        new(
            "",
            Vector{String}(),
            Vector{String}(),
            Dict{String,Vector{String}}(),
            Dict{String,Bool}()
            )
    end
    function Pgf(args...)
		new(
			args...
		)
	end
end

pgfs=Dict{String,AffyCelFiles.Pgf}()
pgfs_files=Dict{String,String}()
pgfs_checksums=Dict{String,UInt32}()

#https://media.affymetrix.com/support/developer/powertools/changelog/file-format-clf.html

struct Clf
    chip_name::String
    header_tags::Vector{String}
    header_values::Vector{String}
    coordinates::Dict{String,Tuple{Int,Int}}    #probe_id -> (x,y)
    max_x::Int
    max_y::Int
    function Clf()
        new(
            "",
            Vector{String}(),
            Vector{String}(),
            Dict{String,Tuple{Int,Int}}(),
            -1,
            -1
            )
    end
    function Clf(args...)
		new(
			args...
		)
	end
end

clfs=Dict{String,AffyCelFiles.Clf}()
clfs_files=Dict{String,String}()
clfs_checksums=Dict{String,UInt32}()

#https://media.affymetrix.com/support/developer/powertools/changelog/file-format-mps.html

struct Mps
    chip_name::String
    header_tags::Vector{String}
    header_values::Vector{String}
    other_ids::Dict{String,Vector{NamedTuple}} # probeset_id (from probeset_list) => (meta_probeset_id1,transcript_cluster_id1),(meta_probeset_id2,transcript_cluster_id2),...
    meta::Dict{String,Vector{String}}          # meta_probeset_id => list of probeset_id's
    function Mps()
        new(
            "",
            Vector{String}(),
            Vector{String}(),
            Dict{String,Vector{NamedTuple}}(),
            Dict{String,Vector{String}}()
            )
    end
    function Mps(args...)
		new(
			args...
		)
	end
end

mpss=Dict{String,AffyCelFiles.Mps}()
mpss_files=Dict{String,String}()
mpss_checksums=Dict{String,UInt32}()

function get_chiptype(cel::Cel)
    chiptype=""
    chiptype_index=findfirst(x->x=="affymetrix-array-type",cel.dataHeader.names)
    if chiptype_index == Nothing
        @warn "chiptype is not part of headers in cel data"
        chiptype="unknown"
    else
        chiptype=cel.dataHeader.values[chiptype_index]
    end
    return chiptype
end

function get_intensities(cel::Cel)
    intensities=Vector{Float32}()
    dataGroupIndex=0
    dataSetIndex=0
    dataSetColIndex=0
    for dataGroupsIndex in eachindex(cel.dataGroups)
        for dataSetsIndex in eachindex(cel.dataGroups[dataGroupsIndex].dataSets)
            if cel.dataGroups[dataGroupsIndex].dataSets[dataSetsIndex].dataSetName == "Intensity"
                dataGroupIndex=dataGroupsIndex
                dataSetIndex=dataSetsIndex
                for dataSetColNamesIndex in eachindex(cel.dataGroups[dataGroupIndex].dataSets[dataSetIndex].dataSetColNames)
                    if cel.dataGroups[dataGroupIndex].dataSets[dataSetIndex].dataSetColNames[dataSetColNamesIndex] == "Intensity"
                        dataSetColIndex=dataSetColNamesIndex
                        break
                    end
                end
                break
            end
        end
        if dataGroupIndex>0
            break
        end
    end
    if dataGroupIndex > 0 && dataSetIndex > 0 && dataSetColIndex > 0
        intensities=cel.dataGroups[dataGroupIndex].dataSets[dataSetIndex].dataSetColumns[dataSetColIndex]
    end
    return intensities
end

function intensities(cel::Cel, pgf::Pgf, clf::Clf, mps::Mps=Mps())
    pm=Dict{String,Vector{Float32}}()
    mm=Dict{String,Vector{Float32}}()

    cel_chiptype=AffyCelFiles.get_chiptype(cel)
    clf_chiptype=clf.chip_name
    pgf_chiptype=pgf.chip_name
    if cel_chiptype != clf_chiptype || cel_chiptype != pgf_chiptype || clf_chiptype != pgf_chiptype
        @error "chiptype of cel data ("*cel_chiptype*"), clf data ("*clf_chiptype*") and pgf data ("*pgf_chiptype*") don't match"
        return (pm=pm,mm=mm)
    end

    use_meta_probeset_id=false
    if ! isempty(mps.chip_name)
        use_meta_probeset_id=true
        if mps.chip_name != cel_chiptype
            @error "chiptype of mps data ("*mps.chip_name*") and cel data ("*cel_chiptype*") don't match"
        end
    end

    intensities=get_intensities(cel)
    if isempty(intensities)
        @error "intensity data not found in cel data"
        return (pm=pm,mm=mm)
    end

    sequential="1"
    index=findfirst(clf.header_tags .== "sequential")
    if index>0
        sequential=clf.header_values[index]
    else
        @error "tag <sequential> not found in clf data"
    end
    order="row_major"
    index=findfirst(clf.header_tags .== "order")
    if index>0
        order=clf.header_values[index]
    else
        @error "tag <order> not found in clf data"
    end
    if sequential != "1" || order != "row_major"
        @warn "unusual values for tags <sequential>="*sequential*" and/or <order>="*order*", expected values are 1 and row_major"
    end
    row_major=true
    if order != "row_major"
        row_major=false
    end
    cols=-1
    index=findfirst(clf.header_tags .== "cols")
    if index>0
        cols=parse(Int,clf.header_values[index])
    else
        @error "tag <cols> not found in clf data"
        cols=max_y+1
    end
    rows=-1
    index=findfirst(clf.header_tags .== "rows")
    if index>0
        rows=parse(Int,clf.header_values[index])
    else
        @error "tag <rows> not found in clf data"
        rows=max_x+1
    end

    check=true
    if sequential=="1" && length(clf.coordinates)>10
        seq=parse(Int,sequential)
        iter=keys(clf.coordinates)
        next=iterate(iter)
        for i = 1:10
            (probe_id, state) = next
            x=clf.coordinates[probe_id][1]
            y=clf.coordinates[probe_id][2]
            #https://media.affymetrix.com/support/developer/powertools/changelog/file-format-clf.html seems to be wrong here.
            #The other way round works:
            if row_major
                test_probe_id = y * cols + x + seq
            else
                test_probe_id = x * rows + y + seq
            end
            if parse(Int,probe_id) != test_probe_id
                check=false
            end
            next = iterate(iter, state)
        end
    end
    if ! check
        @error "sequential/deterministic check failed, double check results is recommended"
    end

    if use_meta_probeset_id
        id_list=eachindex(mps.meta)
    else
        id_list=eachindex(pgf.probe_ids)
    end
    for id in id_list
        pmIntensities=Vector{Float32}()
        mmIntensities=Vector{Float32}()
        if haskey(pm,id) && haskey(mm,id)
            pmIntensities=pm[id]
            mmIntensities=mm[id]
        end
        if use_meta_probeset_id
            probe_id_list=Vector{String}()
            for probeset_id in mps.meta[id]
                for probe_id in pgf.probe_ids[probeset_id]
                    push!(probe_id_list,probe_id)
                end
            end
        else
            probe_id_list=pgf.probe_ids[id]
        end
        for probe_id in probe_id_list
            (x,y)=clf.coordinates[probe_id]
            if row_major
                index = y * cols + x + 1
            else
                index = x * rows + y + 1
            end
            if haskey(pgf.pm,probe_id)
                push!(pmIntensities,intensities[index])
            else
                push!(mmIntensities,intensities[index])
            end
        end
        pm[id]=pmIntensities
        mm[id]=mmIntensities
    end
    return (pm=pm,mm=mm)
end

function intensities(cel::Cel, cdf::Cdf)
    pm=Dict{String,Vector{Float32}}()
    mm=Dict{String,Vector{Float32}}()

    cel_chiptype=AffyCelFiles.get_chiptype(cel)
    
    cdf_chiptype=cdf.chip["Name"]
    if cel_chiptype != cdf_chiptype
        @error "chiptype of cel data ("*cel_chiptype*") and cdf data ("*cdf_chiptype*") don't match"
        return (pm=pm,mm=mm)
    end

    intensities=get_intensities(cel)
    if isempty(intensities)
        @error "intensity data not found in cel data"
        return (pm=pm,mm=mm)
    end

    if ! isempty(cdf.definedUnitIndices)
        for unitIndex in eachindex(cdf.definedUnitIndices)
            unit=cdf.units[unitIndex]
            for blocks in unit.blocks
                for block in blocks
                    id=block.block["Name"]
                    pmIntensities=Vector{Float32}()
                    mmIntensities=Vector{Float32}()
                    if haskey(pm,id) && haskey(mm,id)
                        pmIntensities=pm[id]
                        mmIntensities=mm[id]
                    end
                    for cellsIndex in eachindex(block.cells)
                        cells=block.cells[cellsIndex]
                        pms=block.pm[cellsIndex]
                        for cellIndex in eachindex(cells)
                            index=cells[cellIndex]
                            if pms[cellIndex]
                                push!(pmIntensities,intensities[index+1])
                            else
                                push!(mmIntensities,intensities[index+1])
                            end
                        end
                    end
                    pm[id]=pmIntensities
                    mm[id]=mmIntensities
                end
            end
        end
    else
        @error "no valid cdf data"
    end
    return (pm=pm,mm=mm)
end

#parked code for later
function read_cel_with_cdf(cdf_file)
    cdf=Cdf()
    if isempty(cdf_file) && haskey(AffyCelFiles.cdfs,chiptype) && haskey(AffyCelFiles.cdfs_files,chiptype)
        cdf_file=cdfs_files[chiptype]
    end
    if ! isempty(cdf_file) && ( isdir(cdf_file) || isdirpath(cdf_file) )
        cdf_file = normpath(joinpath(cdf_file,chiptype*".cdf"))
    end
    if isempty(cdf_file)
        cdf_file=chiptype*".cdf"
    end
    if ! isfile(cdf_file)
        if haskey(AffyCelFiles.cdfs,chiptype) && haskey(AffyCelFiles.cdfs_files,chiptype)
            @warn "cdf file " * cdf_file * " does not exist or can't be opened, using "*AffyCelFiles.cdfs_files[chiptype]
            cdf_file=AffyCelFiles.cdfs_files[chiptype]
            cdf=cdf_read(cdf_file)
        else
            @error "cdf file " * cdf_file * " does not exist or can't be opened"
        end
    end
end




function mps_read(mps_file::AbstractString)
    mps_file=normpath(mps_file)
    if isfile(mps_file)
        new_crc=open(crc32c,mps_file)
        io=open(mps_file,"r")
    	mps=mps_read(io,new_crc)
    	close(io)
    
        if length(mps.chip_name)>0 && ! haskey(mpss,mps.chip_name)
            mpss[mps.chip_name]=mps
            mpss_files[mps.chip_name]=mps_file
            mpss_checksums[mps.chip_name]=open(crc32c,mps_file)
        end
    else
        @error "Mps file "*mps_file*" doesn't exist"
        mps=Mps()
    end
    return mps
end

function mps_read(io::IO, new_crc)::Mps
	seekstart(io)
    chip_name=""
    header_tags=Vector{String}()
    header_values=Vector{String}()
    other_ids=Dict{String, Vector{NamedTuple}}()
    meta=Dict{String,Vector{String}}()
    colnames=["probeset_id","transcript_cluster_id","probeset_list","probe_count"]
    transcript_cluster_id_index=2
    probeset_list_start=2   #2 tags before probeset_list
    probeset_list_end=1     #1 tag after probeset_list
    for line in eachline(io)
        if startswith(line,"#%")
            line=strip(line)
            (tag,value)=split(line,"=")
            tag=lstrip(tag,['#','%'])
            push!(header_tags,tag)
            push!(header_values,value)
            if tag=="lib_set_name"
                chip_name=value
                if haskey(mpss,value)
                    if new_crc == mpss_checksums[value] 
                        @info "mps for "*value*" already in cache, source "*mpss_files[value]
                        return mpss[value]
                    else
                        @info "mps for "*value*" already in cache, but source "*mpss_files[value]*" has changed"
                        delete!(mpss,value)
                        delete!(mpss_files,value)
                        delete!(mpss_checksums,value)
                    end
                end
            end
        else
            line=strip(line)
            if startswith(line,"probeset_id")
                colnames=split(line)
                probeset_list_start=findfirst("probeset_list".==colnames)-1
                probeset_list_end=length(colnames)-findfirst("probeset_list".==colnames)
                transcript_cluster_id_index=findfirst("transcript_cluster_id".==colnames)
                if ! (probeset_list_start>0 && probeset_list_end>0)
                    @error "probeset_list indices wrong, tags before="*string(probeset_list_start)*", tags after="*string(probeset_list_end)
                    probeset_list_start=2
                    probeset_list_end=1
                end
            else
                values=split(line)
                meta_probeset_id=values[1]
                if ! isnothing(transcript_cluster_id_index)
                    transcript_cluster_id=values[transcript_cluster_id_index]
                else
                    transcript_cluster_id=""
                end
                probeset_list=values[(probeset_list_start+1):(length(values)-probeset_list_end)]
                if haskey(meta,meta_probeset_id)
                    @error "ambiguous meta probeset_id in mps data"
                end
                meta[meta_probeset_id]=probeset_list
                for probeset_id in probeset_list
                    if haskey(other_ids,probeset_id)
                        #@error "ambiguous probeset_id in mps data: probeset_id part "*probeset_id*" of multiple meta probesets?"
                        other=other_ids[probeset_id]
                    else
                        other=Vector{NamedTuple}()
                    end
                    push!(other,(meta_probeset_id=meta_probeset_id,transcript_cluster_id=transcript_cluster_id))
                    other_ids[probeset_id]=other
                end
            end
        end
    end
    return Mps(chip_name,header_tags,header_values,other_ids,meta)
end

function clf_read(clf_file::AbstractString)
    clf_file=normpath(clf_file)
    if isfile(clf_file)
        new_crc=open(crc32c,clf_file)
        io=open(clf_file,"r")
    	clf=clf_read(io,new_crc)
    	close(io)
    
        if length(clf.chip_name)>0 && ! haskey(clfs,clf.chip_name)
            clfs[clf.chip_name]=clf
            clfs_files[clf.chip_name]=clf_file
            clfs_checksums[clf.chip_name]=open(crc32c,clf_file)
        end
    else
        @error "Clf file "*clf_file*" doesn't exist"
        clf=Clf()
    end
    return clf
end

function clf_read(io::IO, new_crc)::Clf
	seekstart(io)
    chip_name=""
    header_tags=Vector{String}()
    header_values=Vector{String}()
    coordinates=Dict{String,Tuple{Int,Int}}()
    x_first=true
    max_x=-1
    max_y=-1
    for line in eachline(io)
        if startswith(line,"#%")
            line=strip(line)
            (tag,value)=split(line,"=")
            tag=lstrip(tag,['#','%'])
            push!(header_tags,tag)
            push!(header_values,value)
            if tag=="lib_set_name"
                chip_name=value
                if haskey(clfs,value)
                    if new_crc == clfs_checksums[value] 
                        @info "clf for "*value*" already in cache, source "*clfs_files[value]
                        return clfs[value]
                    else
                        @info "clf for "*value*" already in cache, but source "*clfs_files[value]*" has changed"
                        delete!(clfs,value)
                        delete!(clfs_files,value)
                        delete!(clfs_checksums,value)
                    end
                end
            end
            if tag=="header0"
                values=split(value)
                if values[2]=="y" && values[3]=="x"
                    x_first=false
                    @warn "very unusual column order y,x found but expected x,y in clf file"
                end
            end
        else
            line=strip(line)
            values=split(line)
            if length(values)==3
                if x_first
                    x=parse(Int,values[2])
                    y=parse(Int,values[3])
                else
                    y=parse(Int,values[2])
                    x=parse(Int,values[3])
                end
                coordinates[values[1]]=(x,y)
                if x>max_x
                    max_x=x
                end
                if y>max_y
                    max_y=y
                end
            end
        end
    end

    return Clf(chip_name,header_tags,header_values,coordinates,max_x,max_y)
end

function pgf_read(pgf_file::AbstractString)
    pgf_file=normpath(pgf_file)
    if isfile(pgf_file)
        new_crc=open(crc32c,pgf_file)
        io=open(pgf_file,"r")
    	pgf=pgf_read(io,new_crc)
    	close(io)
    
        if length(pgf.chip_name)>0 && ! haskey(pgfs,pgf.chip_name)
            pgfs[pgf.chip_name]=pgf
            pgfs_files[pgf.chip_name]=pgf_file
            pgfs_checksums[pgf.chip_name]=open(crc32c,pgf_file)
        end
    else
        @error "Pgf file "*pgf_file*" doesn't exist"
        pgf=Pgf()
    end
    return pgf
end

function pgf_read(io::IO, new_crc)::Pgf
	seekstart(io)
    chip_name=""
    header_tags=Vector{String}()
    header_values=Vector{String}()
    probe_ids=Dict{String,Vector{String}}()
    pm=Dict{String,Bool}()
    probeset_id=""
    cur_probe_ids=Vector{String}()
    for line in eachline(io)
        if startswith(line,"#%")
            line=strip(line)
            (tag,value)=split(line,"=")
            tag=lstrip(tag,['#','%'])
            push!(header_tags,tag)
            push!(header_values,value)
            if tag=="lib_set_name"
                chip_name=value
                if haskey(pgfs,value)
                    if new_crc == pgfs_checksums[value] 
                        @info "pgf for "*value*" already in cache, source "*pgfs_files[value]
                        return pgfs[value]
                    else
                        @info "pgf for "*value*" already in cache, but source "*pgfs_files[value]*" has changed"
                        delete!(pgfs,value)
                        delete!(pgfs_files,value)
                        delete!(pgfs_checksums,value)
                    end
                end
            end
        elseif startswith(line,r"\d+\t")
            if ! isempty(cur_probe_ids) && ! isempty(probeset_id)
                probe_ids[probeset_id]=cur_probe_ids
            end
            line=strip(line)
            values=split(line)
            probeset_id=values[1]
            if haskey(probe_ids,probeset_id)
                cur_probe_ids=probe_ids[probeset_id]
            else
                cur_probe_ids=Vector{String}()
            end
        elseif startswith(line,r"\t\d+")
            line=strip(line)
            values=split(line)
            atom_id=values[1]
        elseif startswith(line,r"\t\t\d+")
            line=strip(line)
            values=split(line)
            probe_id=values[1]
            push!(cur_probe_ids,probe_id)
            type=values[2]
            if contains(type,"pm")
                pm[probe_id]=true
            end
        end
    end
    if ! isempty(cur_probe_ids) && ! isempty(probeset_id)
        probe_ids[probeset_id]=cur_probe_ids
    end
    return Pgf(chip_name,header_tags,header_values,probe_ids,pm)
end

function cdf_read(cdf_file::AbstractString)
    cdf_file=normpath(cdf_file)
    if isfile(cdf_file)
        new_crc=open(crc32c,cdf_file)
        io=open(cdf_file,"r")
    	cdf=cdf_read(io,new_crc)
    	close(io)
    
        if length(cdf.chip["Name"])>0 && ! haskey(cdfs,cdf.chip["Name"])
            cdfs[cdf.chip["Name"]]=cdf
            cdfs_files[cdf.chip["Name"]]=cdf_file
            cdfs_checksums[cdf.chip["Name"]]=open(crc32c,cdf_file)
        end
    else
        @error "Cdf file "*cdf_file*" doesn't exist"
        cdf=Cdf()
    end
    return cdf
end

function cdf_read(io::IO, new_crc=0x0)::Cdf
	seekstart(io)
    start_section=strip(readline(io))

    if start_section=="[CDF]"
        #ASCII text format is used by the MAS and GCOS 1.0 software. This was also known as the ASCII version.
        #  https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html#CDF_TEXT
        return cdf_read_ascii(io, new_crc)
    end

    seekstart(io)
    magic=read(io,Int8)
    @warn "Magic number is "*string(magic)
    if magic == 67
        #XDA format is used by the GCOS 1.2 and above software. This was also known as the binary or XDA version.
        #  https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html#CDF_XDA
        return cdf_read_xda(io)
    end

    @warn "Unknown CDF file format"
	return Cdf()
end

function cdf_read_xda(io::IO)
    #not implemented, need example file
    @warn "CDF file xda format not implemented, need example file"
	return Cdf()
end

function cdf_read_ascii(io::IO, new_crc)::Cdf
	seekstart(io)
    section=0
    currentQCindex=0
    currentUnitIndex=0
    currentBlockIndex=0
    indexOfINDEX=0
    indexOfPBASE=0
    indexOfTBASE=0
    cdf=Dict{String,String}()
    chip=Dict{String,String}()
    qc=Vector{Qc}()
    units=Vector{Unit}()
    definedUnitIndices=Dict{Int, Bool}()
    for line in eachline(io)
        line=strip(line)
        found=false
        if contains(line,"=")
            (tag,value)=split(line,"=")
            found=true
        end
        if found
            if section == 1   #[CDF]
                cdf[tag]=value
            elseif section == 2 #[Chip]
                chip[tag]=value
                if tag=="Name" && haskey(cdfs,value)
                    if new_crc == cdfs_checksums[value] 
                        @info "cdf for "*value*" already in cache, source "*cdfs_files[value]
                        return cdfs[value]
                    else
                        @info "cdf for "*value*" already in cache, but source "*cdfs_files[value]*" has changed"
                        delete!(cdfs,value)
                        delete!(cdfs_files,value)
                        delete!(cdfs_checksums,value)
                    end
                end
                if tag=="NumQCUnits"
                    numQCUnits=parse(Int,value)
                    qc=Vector{Qc}(undef,numQCUnits)
                end
                if tag=="MaxUnit"
                    maxUnits=parse(Int,value)
                    units=Vector{Unit}(undef,maxUnits)
                end
            elseif section == 3 #[QCn]
                qc[currentQCindex].qc[tag]=value
                if tag == "NumberCells"
                    numberCells=parse(Int,value)
                    qc[currentQCindex].cells[1]=Vector{Int}(undef,numberCells)
                end
                if startswith(tag,"CellHeader")
                    cols=split(value,"\t")
                    indexOfINDEX=findfirst(isequal("INDEX"),cols)
                end
                if startswith(tag,"Cell") && indexOfINDEX > 0
                    cell_index=-1
                    m=match(r"Cell(\d+)",tag)
                    if m !== nothing
                        cell_index=parse(Int,m.captures[1])
                        cols=split(value,"\t")
                        qc[currentQCindex].cells[1][cell_index]=parse(Int,cols[indexOfINDEX])
                    end                   
                end
            elseif section == 4 #[Unitn]
                units[currentUnitIndex].unit[tag]=value
                if tag == "NumberBlocks"
                    numberBlocks=parse(Int,value)
                    units[currentUnitIndex].blocks[1]=Vector{Block}(undef,numberBlocks)
                end
            elseif section == 5 #[Unitn_Blockm]
                if ! isassigned(units[currentUnitIndex].blocks[1],currentBlockIndex)
                    units[currentUnitIndex].blocks[1][currentBlockIndex]=Block() 
                end
                units[currentUnitIndex].blocks[1][currentBlockIndex].block[tag]=value
                if tag == "NumCells"
                    numCells=parse(Int,value)
                    units[currentUnitIndex].blocks[1][currentBlockIndex].cells[1]=Vector{Int}(undef,numCells)
                    units[currentUnitIndex].blocks[1][currentBlockIndex].pm[1]=Vector{Bool}(undef,numCells)
                end
                if startswith(tag,"CellHeader")
                    cols=split(value,"\t")
                    indexOfINDEX=findfirst(isequal("INDEX"),cols)
                    indexOfPBASE=findfirst(isequal("PBASE"),cols)
                    indexOfTBASE=findfirst(isequal("TBASE"),cols)
                end
                if startswith(tag,"Cell") && indexOfINDEX > 0 && indexOfPBASE > 0 && indexOfTBASE > 0
                    cell_index=-1
                    m=match(r"Cell(\d+)",tag)
                    if m !== nothing
                        cell_index=parse(Int,m.captures[1])
                        cols=split(value,"\t")
                        units[currentUnitIndex].blocks[1][currentBlockIndex].cells[1][cell_index]=parse(Int,cols[indexOfINDEX])
                        units[currentUnitIndex].blocks[1][currentBlockIndex].pm[1][cell_index]= cols[indexOfPBASE]!=cols[indexOfTBASE]
                    end                   
                end
            end
        end
        if startswith(line,"[")
            section=0
        end
        if line=="[CDF]"
            section=1
        elseif line=="[Chip]"
            section=2
        elseif startswith(line,"[QC")
            m=match(r"\[QC(\d+)\]",line)
            if m !== nothing
                indexOfINDEX=-1
                currentQCindex=parse(Int,m.captures[1])
                qc[currentQCindex]=Qc()
                section=3
            end
        elseif startswith(line,"[Unit") && ! contains(line,"Block")
            m=match(r"\[Unit(\d+)\]",line)
            if m !== nothing
                currentUnitIndex=parse(Int,m.captures[1])
                units[currentUnitIndex]=Unit()
                definedUnitIndices[currentUnitIndex]=true
                section=4
            end
        elseif startswith(line,"[Unit") && contains(line,"Block")
            m=match(r"\[Unit(\d+)_Block(\d+)\]",line)
            if m !== nothing
                currentUnitIndex=parse(Int,m.captures[1])
                currentBlockIndex=parse(Int,m.captures[2])
                section=5
            end
        end
    end
    return Cdf(cdf,chip,qc,units,definedUnitIndices)
end

function cel_read(cel_file::AbstractString)
    cel_file=normpath(cel_file)
    if isfile(cel_file)
        io=open(cel_file,"r")
        cel=cel_read(io)
        close(io)
    else
        @error "Cel file "*cel_file*" doesn't exist"
        cel=Cel()
    end
	return cel
end

function cel_read(io::IO)::Cel
	seekstart(io)
    magic=read(io,Int8)
    if magic == 59
        #Command Console generic data file format
        #  https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/generic.html
        return cel_read_generic(io)
    end
    
	seekstart(io)
    magic=read(io,Int32)
    if magic == 64
        #Version 4 Format
        #  https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html
        return cel_read_v4(io)
    end

	@warn "Unknown CEL file version"
	return Cel()
end

function cel_read_v4(io::IO)
    #not implemented, need example file
    @warn "CEL file version 4 not implemented, need example file"
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

    chiptype=""
    chiptype_index=findfirst(x->x=="affymetrix-array-type",dataHeader.names)
    if chiptype_index == Nothing
        @warn "chiptype is not part of headers"
        chiptype="unknown"
    else
        chiptype=dataHeader.values[chiptype_index]
    end
    
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
