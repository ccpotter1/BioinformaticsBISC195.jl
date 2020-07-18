module Assignment07

export normalizeDNA,
        composition,
        gc_content,
        complement,
        reverse_complement,
        parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end


# Your code here.
# Don't forget to export your functions!


function composition(seq)
    seq = normalizeDNA(seq) # make uppercase string, check invalid bases
    
    base_dict = Dict("A" =>  0, "C" => 0, "G" => 0, "T"=> 0, "N" => 0) #create dictionary w/ keys A,C,T,G,N
    for base in seq
         if base == 'A'
            base_dict["A"] = base_dict["A"] +1
         elseif base == 'C'
            base_dict["C"] = base_dict["C"]+1
         elseif base == 'G'
            base_dict["G"] = base_dict["G"]+1
        elseif base == 'T'
            base_dict["T"] = base_dict["T"]+1
         elseif base == 'N'
            base_dict["N"] = base_dict["N"]+1
         end
        end
    return base_dict
    
end


function gc_content(seq)
    base_dict=composition(seq)
    c=base_dict["C"]
    g=base_dict["G"]
    
    total=base_dict["N"] + base_dict["G"] + base_dict["T"] + base_dict["C"] + base_dict["A"]   #same as total=content[1]+content[2]+content[3]+content[4]
    #c in second position of the tuple 
    #g in third position
    #total is sum of all positions
    return (c+g)/total
end

 function complement(seq)

end # module Assignment07