# ██████   █████   ██████ ██   ██  █████   ██████  ███████ ███████
# ██   ██ ██   ██ ██      ██  ██  ██   ██ ██       ██      ██
# ██████  ███████ ██      █████   ███████ ██   ███ █████   ███████
# ██      ██   ██ ██      ██  ██  ██   ██ ██    ██ ██           ██
# ██      ██   ██  ██████ ██   ██ ██   ██  ██████  ███████ ███████

using DifferentialEquations
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





# ██████   █████  ██████   █████  ███    ███ ███████ ████████ ███████ ██████  ███████
# ██   ██ ██   ██ ██   ██ ██   ██ ████  ████ ██         ██    ██      ██   ██ ██
# ██████  ███████ ██████  ███████ ██ ████ ██ █████      ██    █████   ██████  ███████
# ██      ██   ██ ██   ██ ██   ██ ██  ██  ██ ██         ██    ██      ██   ██      ██
# ██      ██   ██ ██   ██ ██   ██ ██      ██ ███████    ██    ███████ ██   ██ ███████

#       ██ ███    ██ ██ ████████ ██  █████  ██
#       ██ ████   ██ ██    ██    ██ ██   ██ ██
# █████ ██ ██ ██  ██ ██    ██    ██ ███████ ██
#       ██ ██  ██ ██ ██    ██    ██ ██   ██ ██
#       ██ ██   ████ ██    ██    ██ ██   ██ ███████

# Initial parameters used to initiate simulations and for reference when simulations were interrupted

## ARGS
# - simulations use ARGS provided to Julia when script is called through command lines
# - there is two ARGS used. 
# 	- a job ID (SLURM_ARRAY_JOB_ID) that identifies the job, allowing to identify correct output/error files 
#	- a task ID (SLURM_ARRAY_TASK_ID) that identifies the task and corrosponds to the realwebID (see 2_create_parameters.jl) 
# --> should no ARGS be provided, script will use default job ID "xxx"; task ID can be changed for testing different scenarios
args = length(ARGS) > 0 ? ARGS : ["000", "1"]



## Parameters
# - load minimized parameters
name_startparms = string("dat/1_realweb/realparms_", args[2], ".jld2") 
@load name_startparms parms_min

# - expand and clean parameters
z0, p, p_min = parms_min |>
	x -> get_parms_max(start_min = x) |>
	x -> sel_parms(u_addrow(x[1]), x[2], parms_min[2])




#       ███████  █████  ██    ██ ███████ ██████
#       ██      ██   ██ ██    ██ ██      ██   ██
# █████ ███████ ███████ ██    ██ █████   ██   ██
#            ██ ██   ██  ██  ██  ██      ██   ██
#       ███████ ██   ██   ████   ███████ ██████

# If available, load previously saved densities, timestep, and calculation time; prepare parameters for simulation
# - is sensitive to file naming convention and may mess up if unexpected names occure (e.g. containing whitespace at wrong places)


# - specify where simulation output is saved (create folder if not done before) 
savepath = "dat/2_simweb"
mkpath(savepath) # NOTE that this may cause issues when multiple threats try to create folder

# - list all previously saved outputs
files = readdir(savepath)

# - retrieve IDs to select correct simulation output - utilizes file naming convention (sim_jobID-realwebID-savestepID.jld2)
# 	- savestepID: ID counting up for each time simulation of one parameterset is started (identify last saved simulation output)
savestepID = files |>
	x -> map(y -> split(y, "_")[end], x) |>
	x -> map(y -> split(y, "-")[end], x) |>
	x -> map(y -> split(y, ".")[1], x) |>
	x -> parse.(Int, x)
#	- realwebID: ID as used in "dat/1_realweb" (see also 2_create_parameters.jl)
realwebID = files |>
	x -> map(y -> split(y, "_")[end], x) |>
	x -> map(y -> split(y, "-")[2], x)

realwebID_sel = realwebID .== string(args[2])


# - load and prepare parameters for simulation
if sum(realwebID_sel) > 0
# 	- if available, load available data from previous simulations and adapt simulation parameters to potential extinctions

	name_saved = realwebID_sel |>
		x -> files[x][findall(savestepID[x] .== maximum(savestepID[x]))][1] |>
		x -> string(savepath, "/", x)
	@load name_saved u t calctime

	z02, p2, filter_min_tmp = u_addrow(u[end, :]) |>
		x -> sel_parms(x, p, p_min, initial = false)

else
# 	- if not available, initialize output and simulation parameters without changing them

	u = []
	t = []
	calctime = 0.

	z02, p2, filter_min_tmp = sel_parms(z0, p, p_min, initial = false)

end




# ███████ ██ ███    ███ ██    ██ ██       █████  ████████ ███████
# ██      ██ ████  ████ ██    ██ ██      ██   ██    ██    ██
# ███████ ██ ██ ████ ██ ██    ██ ██      ███████    ██    █████
#      ██ ██ ██  ██  ██ ██    ██ ██      ██   ██    ██    ██
# ███████ ██ ██      ██  ██████  ███████ ██   ██    ██    ███████

#       ███████ ███████ ████████ ██    ██ ██████
#       ██      ██         ██    ██    ██ ██   ██
# █████ ███████ █████      ██    ██    ██ ██████
#            ██ ██         ██    ██    ██ ██
#       ███████ ███████    ██     ██████  ██

# Setup for simulations

# - load functions used in simulations
include("1.3_ODEdynamics_functions.jl")

# - define last timestep 
tmax = 50000

# - define stepsize for saving simulations
tstep = 100



if length(t) != 0 
	# - if simulation picks up from previous simulation (i.e. there is values stored for t)
	if (t[end] < tmax) & (sum(u[end, p2[17]]) != 0.0)
		# - if simulation didn't yet reach last step or all plants went extinct 

		# - define time span for simulations until saving (unless extinctions occure) - aims at saving at regular time steps tstep
		tend = (t[end]-t[end]%tstep)+tstep
		tspan = (t[end], tend) 


		# - define at which timesteps solver records densities - reduces filesizes, especially for cyclic dynamics 
		# 	- initially uses smaller timesteps to allow investigating intial dynamics 
		#	- uses bigger timesteps in the middle to keep filesizes reasonable 
		#	- for last 1000 timesteps, save at every step to calculate unbiased means for analyses

		if t[end] < 1000 
			saveat_vec = (t[end]-t[end]%2.0)+2.0 |> tstart_new -> collect(tstart_new:2.0:tend)
		elseif t[end] < 3000
			saveat_vec = (t[end]-t[end]%10.0)+10.0 |> tstart_new -> collect(tstart_new:10.0:tend)
		elseif t[end] < tmax-1000
			saveat_vec = (t[end]-t[end]%100.0)+100.0 |> tstart_new -> collect(tstart_new:100.0:tend)
		else
			saveat_vec = (t[end]-t[end]%1.0)+1.0 |> tstart_new -> collect(tstart_new:1.0:tend)
		end
	else
		# - if simulation reached last step, take no further steps and give warning
		@warn("Simulation reached tmax or all plants died !!!")
		tspan = (tmax, tmax)
		saveat_vec = []
	end
else 
	# - if it's the first simulation for the specific parameter set
	tspan = (0, tstep) 
	saveat_vec = collect(0:2.0:tstep)

end

# - define extinction callback used to interrupt simulation when extinction occures 
condition2(u, t, integrator) = any((u[1, :] .<= integrator.p[20][1, :]) .& (u[1, :] .!= 0.0))
extinction_cb = DiscreteCallback(condition2, terminate!)



# - define problem to initialize simulation
prob = ODEProblem(foodwebsim, z02, tspan, p2)



#       ██ ███    ██ ██ ████████ ██  █████  ██
#       ██ ████   ██ ██    ██    ██ ██   ██ ██
# █████ ██ ██ ██  ██ ██    ██    ██ ███████ ██
#       ██ ██  ██ ██ ██    ██    ██ ██   ██ ██
#       ██ ██   ████ ██    ██    ██ ██   ██ ███████

# Initialize and save simulations

# - initial simulation
sol = @timed solve(prob, RK4(), abstol=1e-6, reltol=1e-3, maxiters = 1e8, saveat = saveat_vec, callback = extinction_cb)

# - prepare output
u = hcat(map(x -> x[1, :], sol[1].u)...) |>
	x -> collect(transpose(x)) |>
	x -> u_long2(x, length(p_min[end]), filter_min_tmp) |>
	x -> typeof(u) == Matrix{Float64} ? vcat(u, x) : x

t = vcat(t, sol[1].t)

calctime = calctime + sol[2]


# - save output
savestepID_new = savestepID[realwebID_sel] |> x -> length(x) > 0 ? maximum(x)+1 : 1
name_out = string(savepath, "/sim_", 
				  join([args..., savestepID_new], "-"),
				  ".jld2")

@save name_out u t p_min calctime





#        ██████  ██████  ███    ██ ████████ ██ ███    ██ ██    ██ ███████
#       ██      ██    ██ ████   ██    ██    ██ ████   ██ ██    ██ ██
# █████ ██      ██    ██ ██ ██  ██    ██    ██ ██ ██  ██ ██    ██ █████
#       ██      ██    ██ ██  ██ ██    ██    ██ ██  ██ ██ ██    ██ ██
#        ██████  ██████  ██   ████    ██    ██ ██   ████  ██████  ███████

# Continue simulations until done or all plants went extinct

# - to avoid getting sucked into while loops, we define a maximum number of steps the loop is allowed to take
#	- maximum remaining steps (without extinctions) + possible extinctions (without plants and resources)
#	--> can still be a lot of steps !
maxloop = (tmax/tstep - 1) + (size(z0)[2] - 64*3) 

for i in 1:maxloop

	global z02, p2, p_min, tspan, prob, saveat_vec, sol, u, t, calctime, filter_min_tmp, tend, extinction_cb

	if (t[end] < tmax) & (sum(u[end, p2[17]]) != 0.0)
		# - if we didn't reach tmax and not all plants are dead, continue calculating

		# - remove species below extinction threshold! 
		u[end, u[end, :] .<= p[20][1, :]] .= 0.

		z02, p2, filter_min_tmp = u_addrow(u[end, :]) |>
			x -> sel_parms(x, p, p_min, initial = false)

		tend = (t[end]-t[end]%tstep)+tstep
		tspan = (t[end], tend)

		if t[end] < 1000 
			saveat_vec = (t[end]-t[end]%2.0)+2.0 |> tstart_new -> collect(tstart_new:2.0:tend)
		elseif t[end] < 3000
			saveat_vec = (t[end]-t[end]%10.0)+10.0 |> tstart_new -> collect(tstart_new:10.0:tend)
		elseif t[end] < tmax-1000
			saveat_vec = (t[end]-t[end]%100.0)+100.0 |> tstart_new -> collect(tstart_new:100.0:tend)
		else
			saveat_vec = (t[end]-t[end]%1.0)+1.0 |> tstart_new -> collect(tstart_new:1.0:tend)
		end

		prob = ODEProblem(foodwebsim, z02, tspan, p2)

		sol = @timed solve(prob, RK4(), abstol=1e-6, reltol=1e-3, maxiters = 1e8, saveat = saveat_vec, callback = extinction_cb)

		u = hcat(map(x -> x[1, :], sol[1].u)...) |>
			x -> collect(transpose(x)) |>
			x -> u_long2(x, length(p_min[end]), filter_min_tmp) |>
			x -> typeof(u) == Matrix{Float64} ? vcat(u, x) : x

		t = vcat(t, sol[1].t)

		calctime = calctime + sol[2]

		@save name_out u t p_min calctime
	end
end



