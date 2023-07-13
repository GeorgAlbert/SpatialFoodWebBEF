# ███████ ██  ██████  ███    ███  █████      ██      ██ ███    ███ ██ ████████
# ██      ██ ██       ████  ████ ██   ██     ██      ██ ████  ████ ██    ██
# ███████ ██ ██   ███ ██ ████ ██ ███████     ██      ██ ██ ████ ██ ██    ██
#      ██ ██ ██    ██ ██  ██  ██ ██   ██     ██      ██ ██  ██  ██ ██    ██
# ███████ ██  ██████  ██      ██ ██   ██     ███████ ██ ██      ██ ██    ██

# draw values from normal distribution around mean (mu) and redraw values outside of range of 3 standard deviations (sigma)
function sigmalimit(mu, sigma; n = 1)

	x = rand(Distributions.Normal(mu, sigma), n)

	while sum((x .< mu - 3 * sigma) .| (x .> mu + 3 * sigma)) > 0
		tmp = ((x .< mu - 3 * sigma) .| (x .> mu + 3 * sigma))

		x[tmp] = rand(Distributions.Normal(mu, sigma), sum(tmp))
	end

	return x
end


# ██      ██ ███    ██ ██   ██     ██████   █████  ████████  ██████ ██   ██
# ██      ██ ████   ██ ██  ██      ██   ██ ██   ██    ██    ██      ██   ██
# ██      ██ ██ ██  ██ █████       ██████  ███████    ██    ██      ███████
# ██      ██ ██  ██ ██ ██  ██      ██      ██   ██    ██    ██      ██   ██
# ███████ ██ ██   ████ ██   ██     ██      ██   ██    ██     ██████ ██   ██

#    ████████ ██████   ██████  ██████  ██   ██ ██  ██████
#       ██    ██   ██ ██    ██ ██   ██ ██   ██ ██ ██
# █████ ██    ██████  ██    ██ ██████  ███████ ██ ██
#       ██    ██   ██ ██    ██ ██      ██   ██ ██ ██
#       ██    ██   ██  ██████  ██      ██   ██ ██  ██████

# create matrix that determines between which patches across layers trophic interactions can occur
function linklayer_patch(xy_link, nA, nP, nN)

	n8 = size(unique(xy_link[:, 1:2], dims = 1), 1)
	n4 = size(unique(xy_link[:, 3:4], dims = 1), 1)
	n2 = size(unique(xy_link[:, 5:6], dims = 1), 1)
	n1 = size(unique(xy_link[:, 7:8], dims = 1), 1)

	# use large matrix to show between which patches across layers trophic interactions are possible
	# - initialize matrix without any interactions (i.e. all values = 0)
	linklayer_patch = zeros((n8+n4+n2+n1)*(nA+nP+nN), (n8+n4+n2+n1)*(nA+nP+nN))

	# - within patch, feeding interaction will always be possible
	patch_index = map(x -> collect(1:nA+nP+nN).+(x-1)*(nA+nP+nN), 1:(n8+n4+n2+n1))
	for i in patch_index
		linklayer_patch[i, i] .= ones()
	end

	# - between layers, interactions are defined by nested structure (following xy_link)
	for i in 1:n8
		# identify nesting patches with edge size 4
		tmp_patch_index_edge4 = intersect(
			findall(x -> x == xy_link[i, 3], unique(xy_link[:, 3:4], dims = 1)[:, 1]),
			findall(x -> x == xy_link[i, 4], unique(xy_link[:, 3:4], dims = 1)[:, 2])
			)
		# identify nesting patches with edge size 2
		tmp_patch_index_edge2 = intersect(
			findall(x -> x == xy_link[i, 5], unique(xy_link[:, 5:6], dims = 1)[:, 1]),
			findall(x -> x == xy_link[i, 6], unique(xy_link[:, 5:6], dims = 1)[:, 2])
			)

		# 8 to 4
		linklayer_patch[patch_index[i], patch_index[n8+tmp_patch_index_edge4[1]]] .= ones()
		# 8 to 2
		linklayer_patch[patch_index[i], patch_index[n8+n4+tmp_patch_index_edge2[1]]] .= ones()
		# 4 to 2
		linklayer_patch[patch_index[n8+tmp_patch_index_edge4[1]], patch_index[n8+n4+tmp_patch_index_edge2[1]]] .= ones()
	end
		# all to 1
	linklayer_patch[:, patch_index[end]] .= ones()

	return linklayer_patch
end



#       ███    ██ ██    ██ ████████      ██████  ██████  ███    ███ ██████  ███████ ████████ ███████
#       ████   ██ ██    ██    ██        ██      ██    ██ ████  ████ ██   ██ ██         ██    ██
# █████ ██ ██  ██ ██    ██    ██        ██      ██    ██ ██ ████ ██ ██████  █████      ██    █████
#       ██  ██ ██ ██    ██    ██        ██      ██    ██ ██  ██  ██ ██      ██         ██    ██
#       ██   ████  ██████     ██         ██████  ██████  ██      ██ ██      ███████    ██    ███████

# construct adjacency matrix based on x and y coordinates
# - only considers direct neighbors in x or y direction (no diagonal!)
# - uses periodic boundary conditions !!!
function adjacency(cell_x, cell_y)

	dist_x = [abs(i - j) for i in cell_x, j in cell_x] |>
		x -> convert(Matrix{Float64}, x) |>
		x -> replace(y -> y == maximum(x) ? 1 : y, x) |>
		x -> replace(y -> y > 1 ? NaN : y, x)

	dist_y = [abs(i - j) for i in cell_y, j in cell_y] |>
		x -> convert(Matrix{Float64}, x) |>
		x -> replace(y -> y == maximum(x) ? 1 : y, x) |>
		x -> replace(y -> y > 1 ? NaN : y, x)

	adjacency = ((dist_x .== 1) + (dist_y .== 0) .== 2) +
		((dist_x .== 0) + (dist_y .== 1) .== 2)

	if length(cell_x) == 1
		adjacency .= 1
	end

	return adjacency
end


# determines proportion of resources at patches accessible to plant at a given patch (according to N_ikk) 
function rescomp(nA, nP, nN, nPatches, N_ikk, IDP_sel, IDN_sel, adjacency_edge)

	patch_index = map(x -> collect(1:nA+nP+nN).+(x-1)*(nA+nP+nN), 1:nPatches)
	IDN_inpatch = IDN_sel[IDN_sel .∈ Ref(patch_index[1])]

	N_ident = zeros(nPatches*(nA+nP+nN), nPatches*(nA+nP+nN))
	for i in 1:size(adjacency_edge)[1]
		N_ident[IDN_sel[i], IDP_sel[i]] = 1
	end

	N_comp = zeros(nPatches*(nA+nP+nN), nPatches*(nA+nP+nN))
	findall(x -> x == 1, adjacency_edge) |>
		x -> map(y -> N_comp[patch_index[y[1]][IDN_inpatch], patch_index[y[2]]] .= ones(), x)
	N_comp[:, 1:end .∉ Ref(IDP_sel)] .= zeros()

	N = (N_ident .* N_ikk) + (N_comp .* ((1-N_ikk)/4))

	return N
end



# ███████  ██████  ██████  ████████      █████  ███    ██ ██ ███    ███  █████  ██      ███████
# ██      ██    ██ ██   ██    ██        ██   ██ ████   ██ ██ ████  ████ ██   ██ ██      ██
# ███████ ██    ██ ██████     ██        ███████ ██ ██  ██ ██ ██ ████ ██ ███████ ██      ███████
#      ██ ██    ██ ██   ██    ██        ██   ██ ██  ██ ██ ██ ██  ██  ██ ██   ██ ██           ██
# ███████  ██████  ██   ██    ██        ██   ██ ██   ████ ██ ██      ██ ██   ██ ███████ ███████

# assign animals to spatial layers based on bodymass

function sort_spatiallayer(nA, nP, bodymass_animal; nested = true)

	if nested
		# initialize (w/o resources!):
		spatiallayerID = Vector{String}(undef, nA+nP)
		# plants:
		spatiallayerID[1+nA:end] .= "edge8"

		if nA > 0
			# animals: sorted by body mass split every 2 orders of magnitude
			spatiallayerID[findall(x -> x <= 10^2, bodymass_animal)] .= "edge8"
			spatiallayerID[findall(x -> (x > 10^2) & (x <= 10^4), bodymass_animal)] .= "edge4"
			spatiallayerID[findall(x -> (x > 10^4) & (x <= 10^6), bodymass_animal)] .= "edge2"
			spatiallayerID[findall(x -> x > 10^6, bodymass_animal)] .= "edge1"
		end
	
	else
		spatiallayerID = repeat(["edge1"], nA+nP)
		spatiallayerID[1+nA:end] .= "edge8"
	end

	return spatiallayerID
end



# ███████ ███████ ███████ ██████  ██ ███    ██  ██████      ███████ ███████ ███████
# ██      ██      ██      ██   ██ ██ ████   ██ ██           ██      ██      ██
# █████   █████   █████   ██   ██ ██ ██ ██  ██ ██   ███     █████   █████   █████
# ██      ██      ██      ██   ██ ██ ██  ██ ██ ██    ██     ██      ██      ██
# ██      ███████ ███████ ██████  ██ ██   ████  ██████      ███████ ██      ██

#       ██████  ██       █████  ███    ██ ████████    ██   ██ ███████ ██████  ██████  ██
#       ██   ██ ██      ██   ██ ████   ██    ██       ██   ██ ██      ██   ██ ██   ██ ██
# █████ ██████  ██      ███████ ██ ██  ██    ██ █████ ███████ █████   ██████  ██████  ██
#       ██      ██      ██   ██ ██  ██ ██    ██       ██   ██ ██      ██   ██ ██   ██ ██
#       ██      ███████ ██   ██ ██   ████    ██       ██   ██ ███████ ██   ██ ██████  ██

# define plant-herbivore interactions following Thébault & Fontaine 2010:

# Rule1: Modularity
function Rule1(focalID, pmod, module_ident_1, module_ident_2)

	if rand(Binomial(1, pmod)) == 1
		setID = findall(module_ident_2 .== module_ident_1[focalID])
	else
		setID = 1:length(module_ident_2)
	end

	return setID
end

# Rule2: Nestedness
function Rule2(focalID, pnest, pinteract)

	if rand(Binomial(1, pnest)) == 1
		selID = sample(focalID, Weights(pinteract[focalID]./sum(pinteract[focalID])))
	else
		selID = sample(focalID)
	end

	return selID
end

# power law distribution as used by Thebault & Fontaine (2010)
function Powerlawdistribution(maxval)
	y = 1e100
	while y > maxval
		y = rand(Distributions.Uniform(0, 1)) |>
			x -> 0.5 * (1 - x)^-2 + 0.5 |>
			x -> round.(x)
	end
	return y
end


# determine plant-herbivore interactions 
function plantherbivore_interaction(nA, nP; connectance_goal = 0.2, n_module = 4, pnest = 0.2, pmod = 0.7)

	## check if nA and nP are multiple of n_module:
	if mod(nA, n_module) != 0
		@warn("'nA' mod 'n_module' != 0")
	end
	if mod(nP, n_module) != 0
		@warn("'nP' mod 'n_module' != 0")
	end

	## module identifier for each herbivore and plant
	module_ident_herbi = sample(repeat(1:n_module, Int(nA/n_module)), nA, replace = false)
	module_ident_plant = sample(repeat(1:n_module, Int(nP/n_module)), nP, replace = false)

	## probability of interaction
	pinteract_herbi = [Powerlawdistribution(nP) for i in 1:nA]
	pinteract_plant = [Powerlawdistribution(nA) for i in 1:nP]

	## initialize interaction matrix
	interaction = zeros(nP, nA)

	## fill interaction matrix
	while (sum(sum(interaction, dims = 1) .== 0.0) != 0) | (sum(sum(interaction, dims = 2) .== 0.0) != 0)
		# --> while wrapper to make sure there's no unconnected species !

		interaction = zeros(nP, nA)
		connectance = 0.0

		while connectance <= connectance_goal
			# --> while wrapper to add interactions until connectance reaches desired value

			# select plant species:
			plantID = Rule2(1:nP, pnest, pinteract_plant)
			# select herbivore species:
			herbiID = Rule1(plantID, pmod, module_ident_plant, module_ident_herbi) |>
				x -> Rule2(x, pnest, pinteract_herbi)
			# update interaction matrix:
			interaction[plantID, herbiID] = 1.0

			# update connectance:
			connectance = sum(interaction)/(nA*nP)
		end
	end

	return interaction
end


#        ██████  ██████  ███    ███ ██████  ██ ███    ██ ███████ ██████
#       ██      ██    ██ ████  ████ ██   ██ ██ ████   ██ ██      ██   ██
# █████ ██      ██    ██ ██ ████ ██ ██████  ██ ██ ██  ██ █████   ██   ██
#       ██      ██    ██ ██  ██  ██ ██   ██ ██ ██  ██ ██ ██      ██   ██
#        ██████  ██████  ██      ██ ██████  ██ ██   ████ ███████ ██████

# calculate feeding efficiencies for animals following Schneider et al. (2016) and link to plants following Thebault & Fontaine (2010)
# - selection of herbivorous species partially random

function feedingeff(bodymass_animal, nH_omni, nH_pure, nP, nN, optimalmassratio, gamma)

	# make sure bodymass is vector
	bodymass_animal = vec(bodymass_animal)

	nA = length(bodymass_animal)
	nH = nH_pure + nH_omni


	## calculate feeding efficiencies for animals:
	feedingeff = hcat(map(x -> (x ./ bodymass_animal) ./ optimalmassratio, bodymass_animal)...) |>
		x -> (x.*exp.(1 .-x)).^gamma
	# set values below threshold to zero:
	feedingeff[feedingeff .< 0.01] .= zeros()
	# add plants (no interactions):
	feedingeff = feedingeff  |>
		x -> vcat(x, zeros(nP, nA)) |>
		x -> hcat(x, zeros(nA+nP, nP))



	## select and connect nH herbivorous animals
	# select nH_pure pure herbivores (animals without resource species + randomly selected animals)
	IDH = collect(1:nA)[sum(feedingeff, dims = 1)[1:nA] .== 0]
	push!(IDH, sample(collect(1:nA)[1:end .∉ Ref(IDH)], nH_pure-length(IDH), replace = false)...)

	# modify feedingeff to remove carnivorous interactions of selected pure herbivores
	feedingeff[:, IDH] .= zeros()

	# select nH_omni omnivores (randomly choose among the remaining species)
	push!(IDH, sample(collect(1:nA)[1:end .∉ Ref(IDH)], nH-length(IDH), replace = false)...)

	# add herbivorous interactions:
	feedingeff[nA+1:end, IDH] = plantherbivore_interaction(nH, nP, connectance_goal = 0.2,
														   n_module = 4, pnest = 0.2, pmod = 0.7)


	## prepare output: layer in resources (set to zero)
	feedingeff = feedingeff |>
		x -> hcat(x, zeros(nA+nP, nN)) |>
		x -> vcat(x, zeros(nN, nA+nP+nN))

	return feedingeff
end




# ██████  ███████ ███    ██ ███████ ██ ████████ ██ ███████ ███████
# ██   ██ ██      ████   ██ ██      ██    ██    ██ ██      ██
# ██   ██ █████   ██ ██  ██ ███████ ██    ██    ██ █████   ███████
# ██   ██ ██      ██  ██ ██      ██ ██    ██    ██ ██           ██
# ██████  ███████ ██   ████ ███████ ██    ██    ██ ███████ ███████

#       ██ ███    ██ ██ ████████ ██  █████  ██          ██    ██
#       ██ ████   ██ ██    ██    ██ ██   ██ ██          ██    ██
# █████ ██ ██ ██  ██ ██    ██    ██ ███████ ██          ██    ██
#       ██ ██  ██ ██ ██    ██    ██ ██   ██ ██          ██    ██
#       ██ ██   ████ ██    ██    ██ ██   ██ ███████      ██████


# initializes densities for meta food web
function u_single(nA, nP, nN; S = [50, 25])
	# Animals:
	A = rand(Distributions.Uniform(0, 10), nA)
	while sum(A .== 0) > 0
		A[A .== 0] .= rand(Distributions.Uniform(0, 10), sum(A .== 0))
	end
	# Plants:
	P = rand(Distributions.Uniform(0, 10), nP)
	while sum(P .== 0) > 0
		P[P .== 0] .= rand(Distributions.Uniform(0, 10), sum(P .== 0))
	end
	# Resource:
	N = map(x -> rand(Distributions.Uniform(S[x]/2, S[x])), 1:nN)

	return [A... P... N...]
end

# spread densities across patches according to the spatial distribution of species and resources
function u_initial(nA, nP, nN, plant_patchID, spatiallayerID; u_single = u_single(nA, nP, nN, S = [50, 25]))

	# assign densities to patches 
	# - create initial tuples
	u_edge8 = map(x -> copy(u_single), 1:64)
	u_edge4 = map(x -> copy(u_single), 1:16)
	u_edge2 = map(x -> copy(u_single), 1:4)
	u_edge1 = map(x -> copy(u_single), 1:1)


	# - assign starting densities to patches
	for i in 1:64
		u_edge8[i][collect(1:nA+nP)[findall(x -> x != "edge8", spatiallayerID)]] .= zeros() # species
		u_edge8[i][nA.+collect(1:nP)[1:end .!= plant_patchID[i]]] .= zeros() # fix plants (only one per patch!)
	end
	for i in 1:16
		u_edge4[i][collect(1:nA+nP)[findall(x -> x != "edge4", spatiallayerID)]] .= zeros() # species
		u_edge4[i][end-1:end] .= zeros() # remove resources
	end
	for i in 1:4
		u_edge2[i][collect(1:nA+nP)[findall(x -> x != "edge2", spatiallayerID)]] .= zeros() # species
		u_edge2[i][end-1:end] .= zeros() # remove resources
	end
	for i in 1:1
		u_edge1[i][collect(1:nA+nP)[findall(x -> x != "edge1", spatiallayerID)]] .= zeros() # species
		u_edge1[i][end-1:end] .= zeros() # remove resources
	end

	u = hcat(u_edge8..., u_edge4..., u_edge2..., u_edge1...)

	return u
end


#        ██████ ██      ███████  █████  ███    ██     ██    ██
#       ██      ██      ██      ██   ██ ████   ██     ██    ██
# █████ ██      ██      █████   ███████ ██ ██  ██     ██    ██
#       ██      ██      ██      ██   ██ ██  ██ ██     ██    ██
#        ██████ ███████ ███████ ██   ██ ██   ████      ██████

# remove densities of animal species that have no resource species

# - calculate realized number of prey populations 
# 	- can consider matrix filters and densities when not used for meta food web
function preypop(feedeff; filter = ones(size(feedeff)), u = ones(1, size(feedeff)[1]))
	transpose(u .> 0) * (u .> 0) |>
		x -> sum((feedeff .* filter .* x) .> 0, dims = 1)
end

# - remove densities of populations without resource
function u_cleaner!(u, nA, nP, nN, feedeff, filter)

	if nA > 0
		# identifier for animals
		IDA = vcat(map(x -> (x-1)*(nA+nP+nN).+collect(1:nA), 1:85)...)

		while sum((preypop(feedeff, filter = filter, u = u) .> 0) .== (u .> 0)) != ((nA+nP+nN)*85-3*64)

			# identify which species (or resources) have densities and resource species
			tmp = (preypop(feedeff, filter = filter, u = u) .> 0) .== (u .> 0)
			# set densities of animals without resource species to zero
			u[IDA] = tmp[IDA] .* u[IDA]
		end
	end

	return u
end


#       ██       ██████  ███    ██  ██████      ██    ██
#       ██      ██    ██ ████   ██ ██           ██    ██
# █████ ██      ██    ██ ██ ██  ██ ██   ███     ██    ██
#       ██      ██    ██ ██  ██ ██ ██    ██     ██    ██
#       ███████  ██████  ██   ████  ██████       ██████

# expand reduced u to full size
# - this is either 6630 for nA = 60, or 1530 for nA = 0
function u_long(u, nA, filter_min)
	u_full = zeros(size(u)[1], nA == 0 ? 1530 : 6630)
	u_full[:, filter_min] = u

	u = copy(u_full)
end

# expand u to specified size
function u_long2(u, ncol, filter_min)
	u_full = zeros(size(u)[1], ncol)
	u_full[:, filter_min] = u

	return u_full
end



#       ████████ ██████   █████  ███    ██ ███████     ██    ██
#          ██    ██   ██ ██   ██ ████   ██ ██          ██    ██
# █████    ██    ██████  ███████ ██ ██  ██ ███████     ██    ██
#          ██    ██   ██ ██   ██ ██  ██ ██      ██     ██    ██
#          ██    ██   ██ ██   ██ ██   ████ ███████      ██████

# transform starting densities u (single row matrix) to z (multiple row matrix with u in first row and zeros in rest) 
# - used for preallocating memory in ODE calculation
# - handles u as vector or single-row matrix and brings it in correct shape

function u_addrow(u0; addrow = 1)
	z = [length(size(u0)) == 1 ? reshape(u0, 1, size(u0)[1]) : u0] |> 
		x -> vcat(x..., zeros(addrow, length(x...)))
	return z
end




# ███████ ███████ ███████ ██████  ██ ███    ██  ██████      ██████   █████  ████████ ███████ ███████
# ██      ██      ██      ██   ██ ██ ████   ██ ██           ██   ██ ██   ██    ██    ██      ██
# █████   █████   █████   ██   ██ ██ ██ ██  ██ ██   ███     ██████  ███████    ██    █████   ███████
# ██      ██      ██      ██   ██ ██ ██  ██ ██ ██    ██     ██   ██ ██   ██    ██    ██           ██
# ██      ███████ ███████ ██████  ██ ██   ████  ██████      ██   ██ ██   ██    ██    ███████ ███████

#        ██████  █████  ██████  ████████ ██    ██ ██████  ███████      ██████  ██████  ███████ ███████
#       ██      ██   ██ ██   ██    ██    ██    ██ ██   ██ ██          ██      ██    ██ ██      ██
# █████ ██      ███████ ██████     ██    ██    ██ ██████  █████       ██      ██    ██ █████   █████
#       ██      ██   ██ ██         ██    ██    ██ ██   ██ ██          ██      ██    ██ ██      ██
#        ██████ ██   ██ ██         ██     ██████  ██   ██ ███████      ██████  ██████  ███████ ██


# specify diettype based on topology 
function diettype(nA, nP, nN, feedeff)

	linkstructure = feedeff[1:nA+nP, 1:nA+nP]

	producer = [repeat([false], outer = nA)... , repeat([true], outer = nP)...]
	herbivorous = map(x -> sum(linkstructure[producer, x]) > 0, 1:nA+nP)
	carnivorous = map(x -> sum(linkstructure[.!producer, x]) > 0, 1:nA+nP)

	omni  = carnivorous .* herbivorous
	carni = carnivorous .* .!omni
	herbi = herbivorous .* .!omni

	diettype = Vector{String}(undef, nA+nP)

	if sum((producer .+ herbi .+ carni .+ omni) .!= 1) > 0
		@warn("Overlap between trophic groups or unlinked species!")
	else
		diettype[producer] .= "producer"
		diettype[herbi] .= "herbi"
		diettype[omni]  .= "omni"
		diettype[carni] .= "carni"
	end

	push!(diettype, "resource", "resource")
	diettype = repeat(diettype, outer = Int(size(feedeff)[1] / (nA+nP+nN)))

	return diettype
end


# determine scaling constant of capture coefficient based on diettype 
function select_b0(nA, nP, nN, feedeff)
	diettype_tmp = diettype(nA, nP, nN, feedeff) # diets: carni, omni, herbi, producer, resource

	b0 = Vector{Float64}(undef, length(diettype_tmp))
	b0[diettype_tmp .== "producer"] .= zeros()
	b0[diettype_tmp .== "resource"] .= zeros()
	b0[diettype_tmp .== "carni"] .= 50
	b0[diettype_tmp .== "omni"]  .= 100
	b0[diettype_tmp .== "herbi"] .= 400

	return b0
end

# determine scaling exponent of capture coefficient based on diettype 
function select_beta(nA, nP, nN, feedeff)
	diettype_tmp = diettype(nA, nP, nN, feedeff)[1:nA+nP+nN] # diets: carni, omni, herbi, producer, resource

	beta = Vector{Float64}(undef, length(diettype_tmp))
	beta[diettype_tmp .== "producer"] .= zeros()
	beta[diettype_tmp .== "resource"] .= zeros()
	beta[diettype_tmp .== "carni"] = sigmalimit(0.42, 0.05, n = sum(diettype_tmp .== "carni"))
	beta[diettype_tmp .== "omni"]  = sigmalimit(0.19, 0.04, n = sum(diettype_tmp .== "omni"))
	beta[diettype_tmp .== "herbi"] = sigmalimit(0.19, 0.04, n = sum(diettype_tmp .== "herbi"))

	beta = repeat(beta, outer = Int(size(feedeff)[1] / (nA+nP+nN)))

	return beta
end

# calculate capture coefficients
function capturecoef(b0, beta, feedeff, bodymass)
	vec(bodymass).^beta .* transpose(b0 .* vec(bodymass).^beta) .* feedeff
end



#       ██   ██  █████  ███    ██ ██████  ██      ██ ███    ██  ██████      ████████ ██ ███    ███ ███████
#       ██   ██ ██   ██ ████   ██ ██   ██ ██      ██ ████   ██ ██              ██    ██ ████  ████ ██
# █████ ███████ ███████ ██ ██  ██ ██   ██ ██      ██ ██ ██  ██ ██   ███        ██    ██ ██ ████ ██ █████
#       ██   ██ ██   ██ ██  ██ ██ ██   ██ ██      ██ ██  ██ ██ ██    ██        ██    ██ ██  ██  ██ ██
#       ██   ██ ██   ██ ██   ████ ██████  ███████ ██ ██   ████  ██████         ██    ██ ██      ██ ███████

# determine scaling exponent of handling time
function select_eta(nA, nP, nN, feedeff, mu, sigma)
	diettype_tmp = diettype(nA, nP, nN, feedeff)[1:nA+nP+nN] # diets: carni, omni, herbi, producer, resource

	eta = sigmalimit(mu, sigma, n = length(diettype_tmp))
	eta[diettype_tmp .== "producer"] .= zeros()
	eta[diettype_tmp .== "resource"] .= zeros()

	eta = repeat(eta, outer = Int(size(feedeff)[1] / (nA+nP+nN)))

	return eta
end

# calculate handling time
function handling(eta_i, eta_j, bodymass, Fmax, IDP)
	handling = vec(bodymass).^eta_j .* transpose(0.4 .* vec(bodymass).^eta_i)
	handling_Fmax = 1 ./ Fmax |> x -> replace(x, Inf => 0.0)

	handling[IDP, :] .= handling_Fmax[IDP, :]

	return handling
end

#       ██   ██ ██ ██      ██          ███████ ██   ██ ██████
#       ██   ██ ██ ██      ██          ██       ██ ██  ██   ██
# █████ ███████ ██ ██      ██          █████     ███   ██████
#       ██   ██ ██ ██      ██          ██       ██ ██  ██
#       ██   ██ ██ ███████ ███████     ███████ ██   ██ ██ ██

# calculate Hill exponent based on predator-prey bodymass ratios (carnivorous interaction)
# - to avoid calculating values for plants, use bodymass of animals only
function Hillexp(bodymass)
	q = hcat(map(x -> (x ./ vec(bodymass)), vec(bodymass))...) |>
		x -> 1 .* x.^2 ./ ((100)^2 .+ x.^2)
	return q
end



#        ██████  █████  ██       ██████        ██████   █████  ████████ ███████ ███████
#       ██      ██   ██ ██      ██             ██   ██ ██   ██    ██    ██      ██
# █████ ██      ███████ ██      ██             ██████  ███████    ██    █████   ███████
#       ██      ██   ██ ██      ██             ██   ██ ██   ██    ██    ██           ██
#        ██████ ██   ██ ███████  ██████ ██     ██   ██ ██   ██    ██    ███████ ███████

# calculate feeding rates
# - used in ODE
function feedingrate(u, qplus1, int, omega_b, omega_b_h, bodymass)
	transpose(u).^qplus1 |>
		x -> (omega_b .* x) ./
			 (1 .+ int .* u .+ sum(omega_b_h .* x, dims = 1)) ./
			 bodymass  |>
 		x -> replace(x, NaN => 0.0) # NaNs produced for herbi-/carnivorous plants with bodymass = 0
end
