
# AtomsBase interface is now in an extension
# these are the stubs that will be extended. 


# Define the _list_of_properties here for the extension to use
const _list_of_properties = [ 
   (:𝐫, :position),
   (:𝐯, :velocity),
   (:𝑚, :mass),
   (:𝑍, :atomic_number),
   (:S, :species),
   (:𝑞, :charge),
   (:μ, :dipole), 
   (:𝐩, :momentum),
   (:𝐸, :energy),
]

# Define function stubs that will be extended in AtomsBaseExt
function atom end
function symbol end
function aos_system end 
function soa_system end 


# Define setter function stubs that will be extended in AtomsBaseExt
for (sym, name) in _list_of_properties
    fname = Symbol("set_$name")
    fname_ip = Symbol("set_$(name)!")
    fnames_ip = Symbol("set_$(name)s!")
    @eval begin
        function $fname end
        function $fname_ip end
        function $fnames_ip end
        export $fname, $fname_ip, $fnames_ip
    end
end
