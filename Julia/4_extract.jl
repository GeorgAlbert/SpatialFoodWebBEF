# ██████   █████   ██████ ██   ██  █████   ██████  ███████ ███████
# ██   ██ ██   ██ ██      ██  ██  ██   ██ ██       ██      ██
# ██████  ███████ ██      █████   ███████ ██   ███ █████   ███████
# ██      ██   ██ ██      ██  ██  ██   ██ ██    ██ ██           ██
# ██      ██   ██  ██████ ██   ██ ██   ██  ██████  ███████ ███████

using JLD2
using DataFrames
using CSV
using StatsBase
using Distributions
using LinearAlgebra



# ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
# ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
# █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
# ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
# ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████

include("1.1_basefunctions.jl")
include("1.2_ODEparameters_functions.jl")




# ██████   █████  ████████  █████
# ██   ██ ██   ██    ██    ██   ██
# ██   ██ ███████    ██    ███████
# ██   ██ ██   ██    ██    ██   ██
# ██████  ██   ██    ██    ██   ██

# assemble files that should be extracted
# - for each realwebID it takes the last savestepID

savepath = "dat/2_simweb"
files = readdir(savepath)

savestepID = files |>
	x -> map(y -> split(y, "_")[end], x) |>
	x -> map(y -> split(y, "-")[end], x) |>
	x -> map(y -> split(y, ".")[1], x) |>
	x -> parse.(Int, x)
realwebID = files |>
	x -> map(y -> split(y, "_")[end], x) |>
	x -> map(y -> split(y, "-")[2], x) |>
	x -> parse.(Int, x)

files_extract = Vector{Union{Missing, String}}(missing, length(unique(realwebID)))
for (i, j) in enumerate(unique(realwebID))
	savestepID_sel = findall(savestepID[realwebID .== j] .== maximum(savestepID[realwebID .== j]))
	files_extract[i] = files[realwebID .== j][savestepID_sel][1] |>
		x -> string(savepath, "/", x)
end



# ██ ███    ██ ██ ████████ ██  █████  ██          ██████  ███████
# ██ ████   ██ ██    ██    ██ ██   ██ ██          ██   ██ ██
# ██ ██ ██  ██ ██    ██    ██ ███████ ██          ██   ██ █████
# ██ ██  ██ ██ ██    ██    ██ ██   ██ ██          ██   ██ ██
# ██ ██   ████ ██    ██    ██ ██   ██ ███████     ██████  ██

# initialize dataframe
IDcolnames = [
	# ID
	:realwebID, :jobID, :calctime, :tmax,
	:nested, :N_ikk, :nP]
Pcolnames = [
	# Plants
	:densP_mean, :densP_sd, :densP_mean_surv, :densP_sd_surv, :densP_sum,
	:prodP_mean, :prodP_sd, :prodP_mean_surv, :prodP_sd_surv, :prodP_sum,
	:survP_ind, :survP_sp,
	:shannon_ent_dens, :shannon_div_dens, :shannon_ent_prod, :shannon_div_prod]
Acolnames = [
	# Animals
	:nA, :survA, :densA, :fluxHin, :fluxHout]

colnames = [IDcolnames..., Pcolnames..., Acolnames...]
n_col = length(colnames)

df = DataFrame(map(x -> Vector{Union{Missing, Float64}}(missing, length(files_extract)), 1:n_col), colnames)



# ███████ ██   ██ ████████ ██████   █████   ██████ ████████
# ██       ██ ██     ██    ██   ██ ██   ██ ██         ██
# █████     ███      ██    ██████  ███████ ██         ██
# ██       ██ ██     ██    ██   ██ ██   ██ ██         ██
# ███████ ██   ██    ██    ██   ██ ██   ██  ██████    ██


for i in 1:length(files_extract)

	global p_min

	@load files_extract[i] u t p_min calctime

	include("1.3_ODEdynamics_functions.jl")

	println(files_extract[i])


	## ID, meta-data and treatments:
	df.realwebID[i] = split(files_extract[i], "_")[end] |>
		x -> split(x, ".")[1] |>
		x -> split(x, "-")[2] |>
		x -> parse(Float64, x)
	df.jobID[i] = split(files_extract[i], "_")[end] |>
		x -> split(x, ".")[1] |>
		x -> split(x, "-")[1] |>
		x -> parse(Float64, x)

	df.calctime[i] = calctime
	df.tmax[i] = t[end]

	df.nested[i] = p_min[16]
	df.N_ikk[i] = p_min[17]
	df.nP[i] = length(unique(p_min[12]))



	## prepare extraction
	# - standardize u 
	u = u_long(u, p_min[1], p_min[end])
	u[u .< 1e-6] .= zeros()

	# - focus on last 1000 timesteps
	u1000 = u[t .> t[end]-1000, :]
	u1000_mean = mean(u1000, dims = 1)

	# - define filters
	IDA_all = vcat(map(x -> (x-1)*(+(p_min[[1, 4, 5]]...)).+collect(1:p_min[1]), 1:85)...)
	IDA = IDA_all |> x -> x[u[1, x] .> 0]
	IDP = vcat(map(x -> (x-1)*(+(p_min[[1, 4, 5]]...)).+collect(p_min[1]+1:p_min[1]+p_min[4]), 1:64)...) |>
		x -> x[u[1, x] .> 0]
	IDN1 = vcat(map(x -> (x-1)*(+(p_min[[1, 4, 5]]...)).+p_min[1]+p_min[4]+1, 1:64)...)
	IDN2 = vcat(map(x -> (x-1)*(+(p_min[[1, 4, 5]]...)).+p_min[1]+p_min[4]+2, 1:64)...)

	IDPN12 = sort([IDP..., IDN1..., IDN2...])



	## extract data if plants are alive at the end
	if sum(u[end, IDP]) != 0.0 

		## get parameters
		p = get_parms_max(start_min = (reshape(ones(+(p_min[[1, 4, 5]]...)), 1, :), p_min))[2]

		## Plants:
		df.densP_mean[i] = mean(u1000_mean[:, IDP])
		df.densP_sd[i] = std(u1000_mean[:, IDP])
		df.densP_mean_surv[i] = mean(u1000_mean[:, IDP] |> x -> x[x .> 0])
		df.densP_sd_surv[i] = std(u1000_mean[:, IDP] |> x -> x[x .> 0])
		df.densP_sum[i] = sum(u1000_mean[:, IDP])


		N1 = p[8][IDPN12, IDPN12]
		N2 = p[9][IDPN12, IDPN12]

		grow1000_mean = map(sel -> calcgrow(
					reshape(u1000[sel, IDPN12], 1, :), # u
					N1, # N1
					N2, # N2
					p[10], # K1
					p[11], # K2
					collect(1:3:192), collect(2:3:192), collect(3:3:192)),
				1:size(u1000)[1]) |>
			x -> vcat(x...) |>
			grow1000 -> mean(grow1000, dims = 1)

		df.prodP_mean[i] = mean(grow1000_mean)
		df.prodP_sd[i] = std(grow1000_mean)
		df.prodP_mean_surv[i] = mean(grow1000_mean |> x -> x[x .> 0])
		df.prodP_sd_surv[i] = std(grow1000_mean |> x -> x[x .> 0])
		df.prodP_sum[i] = sum(grow1000_mean)

		df.survP_ind[i] = count(u[end, IDP] .> 0)
		df.survP_sp[i] = count(map(x -> sum(u[end, IDP][p_min[12] .== x]) .> 0, unique(p_min[12])))


		df.shannon_ent_dens[i] = map(x -> sum(u1000_mean[:, IDP][p_min[12] .== x]), unique(p_min[12])) |>
			x -> x[x .> 0] |>
			x -> x ./ sum(x) |>
			x -> -sum(x .* log.(x))
		df.shannon_div_dens[i] = exp(df.shannon_ent_dens[i])

		df.shannon_ent_prod[i] = map(x -> sum(grow1000_mean[p_min[12] .== x]), unique(p_min[12])) |>
			x -> x[x .> 0] |>
			x -> x ./ sum(x) |>
			x -> -sum(x .* log.(x))
		df.shannon_div_prod[i] = exp(df.shannon_ent_prod[i])


		## Animals:
		if p_min[1] != 0

			df.nA[i] = reshape(u[1, IDA_all], p_min[1], :) |>
				x -> sum(x, dims = 2) |>
				x -> count(x .> 0)
			df.survA[i] = reshape(u[end, IDA_all], p_min[1], :) |>
				x -> sum(x, dims = 2) |>
				x -> count(x .> 0)
			df.densA[i] = reshape(u1000_mean[:, IDA_all], p_min[1], :) |>
				x -> map(indx -> mean(x[:, collect(indx)], dims = 2)[:], [1:64, 65:80, 81:84, 85]) |>
				x -> sum(sum(x))

			IDherbivorous_sel = findall([sum(p_min[7][61:76, 1:60], dims = 1)...] .> 0) |>
				y -> vcat(map(x -> (x-1)*(+(p_min[[1, 4, 5]]...)).+y, 1:85)...) |>
				x -> x[u1000[1, x] .> .0]
			IDherbivorous_out = findall([sum(p[5][IDherbivorous_sel, :], dims = 1)...] .> 0) |>
				x -> x[u1000[1, x] .> .0]
			IDherbivorous_in = findall([sum(p[5][:, IDherbivorous_sel], dims = 2)...] .> 0) |>
				x -> x[u1000[1, x] .> .0]

			IDfluxH = sort(unique([IDherbivorous_sel..., IDherbivorous_out..., IDherbivorous_in...]))

			feeding1000_mean = (p[2][IDfluxH, IDfluxH],
						   		p[3][:, IDfluxH],
						   		p[4][IDfluxH, IDfluxH],
						   		p[5][IDfluxH, IDfluxH],
						   		p[1][:, IDfluxH]) |>
				x -> map(sel -> reshape(u1000[sel, IDfluxH], 1, :) |>
					u_tmp -> u_tmp * transpose(feedingrate(reshape(u_tmp, 1, :), x...)), 1:size(u1000)[1]) |>
				x -> vcat(x...) |>
				feeding1000 -> mean(feeding1000, dims = 1)

			df.fluxHin[i] = sum(feeding1000_mean[findall(IDfluxH .∈ Ref(IDP))])
			df.fluxHout[i] = sum(feeding1000_mean[findall(IDfluxH .∈ Ref(IDherbivorous_sel))])

		else
			df[i, length(IDcolnames)+length(Pcolnames)+1:end] .= zeros()
		end
	else
		df[i, length(IDcolnames)+1:end] .= zeros()
	end
end


# ███████  █████  ██    ██ ███████
# ██      ██   ██ ██    ██ ██
# ███████ ███████ ██    ██ █████
#      ██ ██   ██  ██  ██  ██
# ███████ ██   ██   ████   ███████


df.file = files_extract
CSV.write("dat/dat.csv", df)
