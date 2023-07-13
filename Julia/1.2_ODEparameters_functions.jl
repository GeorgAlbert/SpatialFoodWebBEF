# ███████ ████████  █████  ██████  ████████     ██████   █████  ██████   █████  ███    ███
# ██         ██    ██   ██ ██   ██    ██        ██   ██ ██   ██ ██   ██ ██   ██ ████  ████
# ███████    ██    ███████ ██████     ██        ██████  ███████ ██████  ███████ ██ ████ ██
#      ██    ██    ██   ██ ██   ██    ██        ██      ██   ██ ██   ██ ██   ██ ██  ██  ██
# ███████    ██    ██   ██ ██   ██    ██        ██      ██   ██ ██   ██ ██   ██ ██      ██ ██

#       ███    ███ ██ ███    ██
#       ████  ████ ██ ████   ██
# █████ ██ ████ ██ ██ ██ ██  ██
#       ██  ██  ██ ██ ██  ██ ██
#       ██      ██ ██ ██   ████

# provide reduced parameter-set for (re-)doing simulations 
# - reduced to save hard disc space (can be 2+ GB)
# - focus on some key and all random parameters
# --> needs to be translated to full parameter set for simulations using get_parms_max()
# - specify treatments:
# 	- N_ikk, defines resource focus of plants on own patch, can have values 0.2, 0.4, 0.6, 0.8, 1.0
# 	- webtreatment, defines food web scenario, can have values none, non-nested, nested

function get_parms_min(; N_ikk = 1., webtreatment = "non-nested")

	## decode food web treatment
	if webtreatment == "none"
		# number of animals
		nA = 0
		# specify nested treatment 
		nested = false
	elseif webtreatment == "non-nested"
		# number of animals
		nA = 60
		# specify nested treatment 
		nested = false
	elseif webtreatment == "nested"
		# number of animals
		nA = 60
		# specify nested treatment 
		nested = true
	else
		@warn("Unknown value for webtreatment. Available values: none, non-nested, nested")
	end

	## Number of herbivores, omnivores, plants and resources
	# number of herbivores and omnivores (irrelevant when nA = 0)
	nH_pure = 16
	nH_omni = 16
	# number of plants
	nP = 16 
	# number of resources
	nN = 2


	## Animal body masses
	bodymass_animal = 10 .^ rand(Distributions.Uniform(log(10, 10^0), log(10, 10^8)), nA)


	## feeding efficiency
	feedeff_global = feedingeff(bodymass_animal, nH_omni, nH_pure, nP, nN, 100, 2)


	## assign plants to patches:
	# - randomly assigns even number of species
	plant_patchID = repeat(1:nP, Int(64 / nP)) |>
		x -> sample(x, 64, replace = false)


	## biomass and resource starting-densities:
	u0 = u_single(nA, nP, nN, S = [50, 25])


	## random parameters for capture coefficient and handling times
	beta = select_beta(nA, nP, nN, feedeff_global)

	eta_i = select_eta(nA, nP, nN, feedeff_global, -0.48, 0.03)
	eta_j = select_eta(nA, nP, nN, feedeff_global, -0.66, 0.02)


	## interference competition
	int = sigmalimit(0.8, 0.2; n = nA) |>
		x -> [x... zeros(nP+nN)...]


	## half saturation density of resource uptake
	K1_sp = reshape(rand(Distributions.Uniform(0.1, 0.2), nP), 1, :)
	K2_sp = reshape(rand(Distributions.Uniform(0.1, 0.2), nP), 1, :)


	## resource concentration in plant biomass - stochiometric parameters
	v1_sp = reshape(sigmalimit(2/3, 0.05, n = nP), 1, :)
	# v2 = 1 .- v1


	## interaction-specific q for plant-herbivore interactions
	q_plantherb = rand(Uniform(0, 1), nP, nA)

	## predator species-specific y_i for F_max
	y_i = sigmalimit(6, 1, n = nA)

	# NOTE: q_plantherb and y_i have more parameters than needed (parameters for carnivores and non-existing feeding interactions !) --> easier to combine in get_parms_max()



	## output vector
	p_min = (
		nA, 					# p_min[1]
		nH_pure,
		nH_omni,
		nP,
		nN,						# p_min[5]
		bodymass_animal,
		feedeff_global,
		int,
		beta,
		eta_i,					# p_min[10]
		eta_j,
		plant_patchID,
		K1_sp,
		K2_sp,
		v1_sp,					# p_min[15]
		nested,
		N_ikk,
		q_plantherb,
		y_i
		)

	return (u0, p_min)

end


#       ███    ███ ██ ███    ██     ███    ███  ██████  ██████  ██ ███████ ██    ██
#       ████  ████ ██ ████   ██     ████  ████ ██    ██ ██   ██ ██ ██       ██  ██
# █████ ██ ████ ██ ██ ██ ██  ██     ██ ████ ██ ██    ██ ██   ██ ██ █████     ████
#       ██  ██  ██ ██ ██  ██ ██     ██  ██  ██ ██    ██ ██   ██ ██ ██         ██
#       ██      ██ ██ ██   ████     ██      ██  ██████  ██████  ██ ██         ██

# modify reduced parameter sets from get_parms_min()
# - used to implement treatments while keeping everything else the same
# - plants can be redistributed in space by setting change_patchIDs to true
# - when change_patchIDs is true, new_patchIDs can be used --> allows implmenting of diversity treatment
#	- new_patchIDs needs to be Vector{Int} and number of patches (i.e. 64) needs to be multiple of length(new_patchIDs)
#	- NOTE: number of plants nP is not updated as it is required for proper sorting of parameters

function modify_parms_min(parms_min; N_ikk = 1., webtreatment = "non-nested", change_patchIDs = false, new_patchIDs = 1:16)
	
	## decode food web treatment
	if webtreatment == "none"
		# number of animals
		nA = 0
		# specify nested treatment 
		nested = false
	elseif webtreatment == "non-nested"
		# number of animals
		nA = 60
		# specify nested treatment 
		nested = false
	elseif webtreatment == "nested"
		# number of animals
		nA = 60
		# specify nested treatment 
		nested = true
	else
		@warn("Unknown value for webtreatment. Available values: none, non-nested, nested")
	end


	## define plant_patchID
	if change_patchIDs
		# check if plant_patchIDs has correct format
		if typeof(new_patchIDs) != Vector{Int}
			@warn("new_patchIDs needs to be type Vector{Int}")
		end
		if mod(64, length(new_patchIDs)) != 0
			@warn("number of patches (i.e. 64) is not multiple of length(new_patchIDs)")
		end

		plant_patchID = repeat(new_patchIDs, Int(64 / length(new_patchIDs))) |>
			x -> sample(x, 64, replace = false)
	else
		plant_patchID = parms_min[2][12]
	end


	## define new p_min
	p_min = (nA,
			parms_min[2][2:11]...,
			plant_patchID,
			parms_min[2][13:15]...,
			nested,
			N_ikk,
			parms_min[2][18:end]...)


	## define output, incl. removing animals from starting densities when required
	if nA == 0
		return (parms_min[1][:, parms_min[2][1]+1:end], p_min)
	else
		return (parms_min[1], p_min)
	end

end


#       ███    ███  █████  ██   ██
#       ████  ████ ██   ██  ██ ██
# █████ ██ ████ ██ ███████   ███
#       ██  ██  ██ ██   ██  ██ ██
#       ██      ██ ██   ██ ██   ██

# creates maximized parameter set 
# - provides starting densities and parameters necessary for simulation
# - expands minimized parameter set from get_parms_min() when provided for argument start_min

function get_parms_max(; start_min = get_parms_min())

	u0, p_min = start_min

	nA, nH_pure, nH_omni, nP, nN, bodymass_animal, feedeff_global, int,	beta, eta_i, eta_j, plant_patchID, K1_sp, K2_sp, v1_sp, nested, N_ikk, q_plantherb, y_i = p_min


	## spatial parameters 
	# number of spatial layers
	nLayer = 4
	# edge length of lowest layer
	nMaxedgelength = 8
	# maximum number of distinct patches (i.e. when nested)
	# - 64 + 16 + 4 + 1 (from plant to community scale)
	# - for food web scenarios non-nested and none, intermediate (16+4) and upper (16+4+1) patches are not populated
	nPatches = 	Int(sum(map(x-> (nMaxedgelength / 2^x)^2, 0:nLayer-1)))


	## Animal body masses
	bodymass_animals_long = [(nA > 0 ? bodymass_animal : [])..., zeros(nP+nN)...] |>
		x -> repeat(x, outer = nPatches) |>
		x -> reshape(x, 1, :)


	## feeding efficiency --> topology !
	# - randomized herbivore-plant interactions !
	# - autotrophic interactions set to zero (calculated differently)
	feedeff = feedeff_global |>
		# repeat to scale to size of filter_nested
		x -> vcat(map(y -> x, 1:nPatches)...) |>
		x -> hcat(map(y -> x, 1:nPatches)...)


	## assign species to spatial layers:
	spatiallayerID = sort_spatiallayer(nA, nP, bodymass_animal, nested = nested)


	## Identity filters
	# - animals
	# IDA = vcat(map(x -> (x-1)*(nA+nP+nN).+collect(1:nA), 1:nPatches)...)
	# - plants
	IDP = vcat(map(x -> (x-1)*(nA+nP+nN).+collect(nA+1:nA+nP), 1:nPatches)...)
	IDP_sel = (collect(1:64) .- 1) .* (nA+nP+nN) .+ nA .+ plant_patchID
	# - resources
	# IDN = vcat(map(x -> (x-1)*(nA+nP+nN).+collect(nA+nP+1:nA+nP+nN), 1:nPatches)...)
	IDN1 = vcat(map(x -> (x-1)*(nA+nP+nN).+collect(nA+nP+1:nA+nP+nN), 1:nPatches)...) |> x -> x[1:2:128]
	IDN2 = vcat(map(x -> (x-1)*(nA+nP+nN).+collect(nA+nP+1:nA+nP+nN), 1:nPatches)...) |> x -> x[2:2:128]


	## Matrix filters
	# preparation
	# - define spatial structure of overlapping patches
	xy_link = hcat(repeat(1:8, inner = 8), repeat(1:8, outer = 8)) |>
		x -> hcat(x, cld.(x[:, 1], 2), cld.(x[:, 2], 2)) |>
		x -> hcat(x, cld.(x[:, 3], 2), cld.(x[:, 4], 2)) |>
		x -> hcat(x, cld.(x[:, 5], 2), cld.(x[:, 6], 2))
	# - adjecency matrices to find neighboring patches
	adjacency_edge8 = xy_link[:, 1:2] |> x -> adjacency(x[:, 1], x[:, 2])
	

	# 1) filter to ensure interactions follow nested structure
	filter_nested = linklayer_patch(xy_link, nA, nP, nN)
	# 2) filter out impossible feeding interactions - based on feeding efficiencies
	filter_feedeff = feedeff .> 0
	# 3) filter out impossible feeding interactions - based on spatial layers occurence
	filter_u = u_initial(nA, nP, nN, plant_patchID, spatiallayerID, u_single = u0) |>
		x -> transpose(x .> 0) * (x .> 0)

	filter_master = filter_nested .* ifelse(nA > 0, filter_feedeff, 1.0) .* filter_u

	## Resource competition filter
	# resource 1
	N1 = IDN1 |> x -> rescomp(nA, nP, nN, nPatches, N_ikk, IDP_sel, x, adjacency_edge8)
	# resource 2
	N2 = IDN2 |> x -> rescomp(nA, nP, nN, nPatches, N_ikk, IDP_sel, x, adjacency_edge8)




	## biomass and resource starting densities:
	u0_new = u_initial(nA, nP, nN, plant_patchID, spatiallayerID, u_single = u0) |>
		x -> u_cleaner!(x, nA, nP, nN, feedeff, filter_master)

	if nA > 0
		## Parameters for feeding rates
		# - relative consumption rate
		omega = 1 ./ preypop(feedeff_global) |>
			x -> replace(x, Inf => 0.0) |>
			x -> hcat(map(y -> x, 1:nPatches)...)

		# - capture coefficient
		b_ij = capturecoef(select_b0(nA, nP, nN, feedeff),
						   repeat(beta, outer = Int(size(feedeff)[1] / (nA+nP+nN))),
						   feedeff, bodymass_animals_long)

		# - maximum feeding rate (for herbivory) based on individuals
	 	Fmax = repeat([y_i..., zeros(nP+nN)...], outer = nPatches) |>
	 		y_i_long -> reshape(y_i_long, 1, :) |>
	 		y_i_long -> y_i_long .* (0.141 .* bodymass_animals_long.^0.695) |>
			y_x -> map(x -> y_x, 1:(nA+nP+nN)*nPatches) |>
			y_x -> vcat(y_x...)

		# - handling time
	 	h_ij = handling(repeat(eta_i, outer = Int(size(feedeff)[1] / (nA+nP+nN))),
						repeat(eta_j, outer = Int(size(feedeff)[1] / (nA+nP+nN))),
						bodymass_animals_long,
						Fmax, IDP)

		# --> combine to avoid additional calculations during simulation
		omega_b_filtered = omega .* b_ij .* filter_master
		omega_b_h_filtered = omega .* b_ij .* h_ij .* filter_master

		# - Hill exponent
		qplus1 = Hillexp(bodymass_animal) |>
			# layer in plants:
			# - herbivorous interactions defined in q_plantherb; rest set to zero!
			x -> vcat(x, q_plantherb) |>
			x -> hcat(x, zeros(nA+nP, nP)) |>
			# layer in resources:
			x -> vcat(x, zeros(nN, nA+nP)) |>
			x -> hcat(x, zeros(nA+nP+nN, nN)) |>
			# repeat to scale to size
			x -> vcat(map(y -> x, 1:nPatches)...) |>
			x -> hcat(map(y -> x, 1:nPatches)...) |>
			# add 1 to avoid calculation in feedingrate()
			x -> x .+ 1

		# - interference competition
		int_long = hcat(map(y -> int, 1:nPatches)...)



		## additional ODE parameters
		# - conversion efficiency - carnivory, herbivory, autotroph (set to 0)
		e = [repeat([0.906], outer = nA)... repeat([0.545], outer = nP)... zeros(nN)...] |>
			x -> hcat(map(y -> x, 1:nPatches)...)

		# - metabolic demands per unit biomass - animals only !
		x_animals = [0.141 .* bodymass_animal.^-0.305... zeros(nP+nN)...] |>
			x -> hcat(map(y -> x, 1:nPatches)...)

	else
		# if no animals, create dummy variables
		omega_b_filtered = zeros(nPatches*(nP+nN), nPatches*(nP+nN))
		omega_b_h_filtered = zeros(nPatches*(nP+nN), nPatches*(nP+nN))
		qplus1 = zeros(nPatches*(nP+nN), nPatches*(nP+nN))
		int_long = zeros(1, nPatches*(nP+nN))
		e = zeros(1, nPatches*(nP+nN))
		x_animals = reshape(zeros(nPatches*(nP+nN)), 1, :)
	end


	## parametes for nutrient dynamics
	# - supply concentration
	S = [50, 25]

	# - turnover rate
	D = 0.25

	# - half saturation density of resource uptake
	K1 = reshape(K1_sp[plant_patchID], 1, :)
	K2 = reshape(K2_sp[plant_patchID], 1, :)

	# - resource concentration in plant biomass - stochiometric parameters
	v1 = reshape(v1_sp[plant_patchID], 1, :)
	v2 = 1 .- v1


	## extinction thresholds based on spatial layer !
	extinction_threshold = [
		repeat([1e-6], outer = (nA+nP+nN)*64)...,
		repeat([4e-6], outer = (nA+nP+nN)*16)...,
		repeat([1.6e-5], outer = (nA+nP+nN)*4)...,
		repeat([6.4e-5], outer = (nA+nP+nN)*1)...
		]
	extinction_threshold = reshape(extinction_threshold, 1, :)



	## combine and output parameters
	p = (
		bodymass_animals_long,	# p[1]
		qplus1,
		int_long,
		omega_b_filtered,
		omega_b_h_filtered,		# p[5]
		x_animals,
		e,
		N1,
		N2,
		K1,						# p[10]
		K2,
		S[1],
		S[2],
		D,
		v1,						# p[15]
		v2,
		IDP_sel,
		IDN1,
		IDN2,
		extinction_threshold	# p[20]
	)

	return (u0_new, p)

end


#       ███████ ███████ ██      ███████  ██████ ████████
#       ██      ██      ██      ██      ██         ██
# █████ ███████ █████   ██      █████   ██         ██
#            ██ ██      ██      ██      ██         ██
#       ███████ ███████ ███████ ███████  ██████    ██

# used to reduce size of parameters by removing non-existing interactions
# - can drastically reduce number of calculations
# - can also be used to remove animal populations that went extinct during simulation (i.e. if inital != true)

function sel_parms(u0, p, p_min; initial = true)

	## minimizing filter
	# - minimize size of parameters by removing species' values that are not needed
	# - never removes plants because they are essential for calculating nutrient dynamcis (even when extinct)
	filter_min = [findall(u0[1, :] .> 0)..., p[17]...] |>
		x -> sort(unique(x))

	## update ID vectors
	IDP_sel_min = findall(filter_min .∈ Ref(p[17]))
	IDN1_min = findall(filter_min .∈ Ref(p[18]))
	IDN2_min = findall(filter_min .∈ Ref(p[19]))


	## assemble new output parameters
	p = (
		p[1][:, filter_min],				# p[1]
		p[2][filter_min, filter_min],
		p[3][:, filter_min],
		p[4][filter_min, filter_min],
		p[5][filter_min, filter_min],		# p[5]
		p[6][:, filter_min],
		p[7][:, filter_min],
		p[8][filter_min,filter_min],
		p[9][filter_min, filter_min],
		p[10:16]...,						# p[10], p[15]
		IDP_sel_min,
		IDN1_min,
		IDN2_min,
		p[20][:, filter_min]				# p[20]
	)

	## define output
	# - if initial == true, p_min gets extended by filter_min (defines minimum number of densities tracked throughout simulation)
	if initial

		p_min = (
			p_min...,
			filter_min
		)

		return (u0[:, filter_min], p, p_min)
	
	# - if intial != true, p_min is not modified nor returned to avoid overwriting initial filter, but filter_min is added to output
	else

		return (u0[:, filter_min], p, filter_min)

	end
end