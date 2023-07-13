# ██████   █████   ██████ ██   ██  █████   ██████  ███████ ███████
# ██   ██ ██   ██ ██      ██  ██  ██   ██ ██       ██      ██
# ██████  ███████ ██      █████   ███████ ██   ███ █████   ███████
# ██      ██   ██ ██      ██  ██  ██   ██ ██    ██ ██           ██
# ██      ██   ██  ██████ ██   ██ ██   ██  ██████  ███████ ███████

using StatsBase
using Distributions
using LinearAlgebra
using JLD2



# ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
# ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
# █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
# ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
# ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████

include("1.1_basefunctions.jl")
include("1.2_ODEparameters_functions.jl")




#  ██████ ██████  ███████  █████  ████████ ███████     ██████   █████  ████████  █████
# ██      ██   ██ ██      ██   ██    ██    ██          ██   ██ ██   ██    ██    ██   ██
# ██      ██████  █████   ███████    ██    █████       ██   ██ ███████    ██    ███████
# ██      ██   ██ ██      ██   ██    ██    ██          ██   ██ ██   ██    ██    ██   ██
#  ██████ ██   ██ ███████ ██   ██    ██    ███████     ██████  ██   ██    ██    ██   ██

#       ███    ███ ███████ ████████  █████  ██     ██ ███████ ██████
#       ████  ████ ██         ██    ██   ██ ██     ██ ██      ██   ██
# █████ ██ ████ ██ █████      ██    ███████ ██  █  ██ █████   ██████
#       ██  ██  ██ ██         ██    ██   ██ ██ ███ ██ ██      ██   ██
#       ██      ██ ███████    ██    ██   ██  ███ ███  ███████ ██████

# create parametersets of meta-food webs 
# - create folder at working directory (check with pwd())
path_metaweb = "dat/0_metaweb"
mkpath(path_metaweb)


# - define IDs for parametersets 
metawebID = collect(1:20)


# - create and save random parametersets 
for i in eachindex(metawebID)
	parms_min = get_parms_min()
	metawebID_tmp = metawebID[i]

	name_out = "$path_metaweb/metaparms_$metawebID_tmp.jld2"
	@save name_out parms_min
end




#       ██████  ███████  █████  ██      ██     ██ ███████ ██████
#       ██   ██ ██      ██   ██ ██      ██     ██ ██      ██   ██
# █████ ██████  █████   ███████ ██      ██  █  ██ █████   ██████
#       ██   ██ ██      ██   ██ ██      ██ ███ ██ ██      ██   ██
#       ██   ██ ███████ ██   ██ ███████  ███ ███  ███████ ██████

# create parametersets across treatments 
# - create folder at working directory (check with pwd())
path_realweb = "dat/1_realweb"
mkpath(path_realweb)

# - get list of meta food web parameter sets
metafiles = readdir(path_metaweb)

# - initialize realwebID
realwebID = 1



# - create parameter sets across treatments
for metaparm in metafiles

	## load metaweb
	name_in = "$path_metaweb/$metaparm"
	@load name_in parms_min
	metaparms_min = deepcopy(parms_min)

	## modify treatments and save parameters
	for webtreatment in ["none", "non-nested", "nested"]
		for N_ikk in [1., 0.8, 0.6, 0.4, 0.2]

			# - 16sp mixtures
			parms_min = modify_parms_min(metaparms_min, N_ikk = N_ikk, webtreatment = webtreatment)
			name_out = "$path_realweb/realparms_$realwebID.jld2"
			@save name_out parms_min
			realwebID += 1

			# - 1sp mixtures
			for sp_focal in 1:16
				parms_min = modify_parms_min(metaparms_min, N_ikk = N_ikk, webtreatment = webtreatment, change_patchIDs = true, new_patchIDs = [sp_focal])
				name_out = "$path_realweb/realparms_$realwebID.jld2"
				@save name_out parms_min
				realwebID += 1	
			end
		end
	end

end


## sorting of treatments in parameter sets:
# - blocks of 17 are diversity treatment (1*16sp mixture, 16*1sp mixture), nested in ...
# - blocks of 5*17=85 are spatial resource overlap treatment (5 treatments for N_ikk),  nested in ...
# - blocks of 3*85=255 are food web treatments (3 scenarios), nested in ...
# - total of 20*255=5100 covering 20 meta food webs