
# ███████ ███████ ███████ ██████     ██    ███    ███ ███████ ████████  █████  ██████   ██████  ██
# ██      ██      ██      ██   ██    ██    ████  ████ ██         ██    ██   ██ ██   ██ ██    ██ ██
# █████   █████   █████   ██   ██ ████████ ██ ████ ██ █████      ██    ███████ ██████  ██    ██ ██
# ██      ██      ██      ██   ██ ██  ██   ██  ██  ██ ██         ██    ██   ██ ██   ██ ██    ██ ██
# ██      ███████ ███████ ██████  ██████   ██      ██ ███████    ██    ██   ██ ██████   ██████  ███████

# select how to calculate feeding and animal metabolism in ODE


if p_min[1] > 0 # p_min[1] is number of animals (nA)

	## function used when there is animals

	function feedmetabol(u, feed, x_all, e)
		## biomass increase due to feeding:
		(e * feed .* u) .- # NOTE: element-wise multiplication with u at end assures correct matrix-multiplication with e (sums up while multiplying)
		## biomass decrease due to predation:
		(u * transpose(feed)) .-
		## biomass decrease due to metabolism:
		(x_all .* u)
	end

else

	## function used whene there is no animals - only returns zeros ! 

	function feedmetabol(u, feed, x_all, e)
		zeros(size(u))
	end

end


#  ██████  ██████   ██████  ██     ██ ████████ ██   ██
# ██       ██   ██ ██    ██ ██     ██    ██    ██   ██
# ██   ███ ██████  ██    ██ ██  █  ██    ██    ███████
# ██    ██ ██   ██ ██    ██ ██ ███ ██    ██    ██   ██
#  ██████  ██   ██  ██████   ███ ███     ██    ██   ██

# select how to calculate growth in ODE


if p_min[17] != 1 # p_min[17] defines plant space-use scenario (N_ikk)

	## function used whene there is competition for resources

	function calcgrow(u, N1, N2, K1, K2, IDP_sel, IDN1, IDN2)
		min.((u * N1)[:, IDP_sel] |> x -> (x ./ (K1 .+ x)),
			 (u * N2)[:, IDP_sel] |> x -> (x ./ (K2 .+ x))) |>
			x -> u[:, IDP_sel].^0.75 .* x
	end

else

	## function used whene there is no competition for resources
	# - should be a bit quicker
	# - has unused arguments (N1, N2) to be used in place of alternative function above 

	function calcgrow(u, N1, N2, K1, K2, IDP_sel, IDN1, IDN2)
		min.(u[:, IDN1] |> x -> (x ./ (K1 .+ x)),
			 u[:, IDN2] |> x -> (x ./ (K2 .+ x))) |>
			x -> u[:, IDP_sel].^0.75 .* x
	end

end


# ███    ██ ██    ██ ████████     ██    ██ ███████ ███████
# ████   ██ ██    ██    ██        ██    ██ ██      ██
# ██ ██  ██ ██    ██    ██        ██    ██ ███████ █████
# ██  ██ ██ ██    ██    ██        ██    ██      ██ ██
# ██   ████  ██████     ██         ██████  ███████ ███████

# select how to calculate nutrient dynamics in ODE


if p_min[17] != 1 # p_min[17] defines plant space-use scenario (N_ikk)

	## function used whene there is competition for resources

	function calcnut(u, grow, Nx, Sx, D, vx, IDP_sel, IDNx)
		## nutrient turnover
		(D .* (Sx .- u[:, IDNx])) .-
		## nutrient loss to plant growth -- relative to nutrient contribution !
		(vx .* grow) * transpose(Nx[IDNx, IDP_sel] .* u[IDNx] |> x -> (x ./ sum(x, dims = 1)) |> x -> replace(x, NaN => 0)) # NOTE: division by zero is possible if all nutrients are zero for at least one plant -- replace with 0
	end

else

	## function used whene there is no competition for resources
	# - should be a bit quicker
	# - has unused arguments (Nx, IDP_sel) to be used in place of alternative function above 

	function calcnut(u, grow, Nx, Sx, D, vx, IDP_sel, IDNx)
		## nutrient turnover
		(D .* (Sx .- u[:, IDNx])) .-
		## nutrient loss to plant growth
		(vx .* grow)
	end

end



#  ██████  ██████  ███████
# ██    ██ ██   ██ ██     
# ██    ██ ██   ██ █████  
# ██    ██ ██   ██ ██     
#  ██████  ██████  ███████

function foodwebsim(dz, z, p, t)

	# define views to be used in calculations (pre-allocated memory)
	u = @view z[1, :]
	du = @view dz[1, :]
	grow = @view dz[2, 1:length(p[17])]

	# reshape data 
	u = reshape(max.(u, .0), 1, :) # max(...) prevents negative values (i.e. if extinction didn't trigger)
	du = reshape(du, 1, :)
	grow = reshape(grow, 1, :)



	# calculate plant growth
	grow .= calcgrow(u, p[8:11]..., p[17:19]...)



	# density changes
	# - animal feeding and metabolism
	du .= feedingrate(u, p[2:5]..., p[1]) |>
		feed -> feedmetabol(u, feed, p[6:7]...)
	# - plant growth and metabolism
	du[:, p[17]] .= du[:, p[17]] .+ grow .- (0.138 .* u[:, p[17]].^0.75)
	# - nutrient dynamics
	du[:, p[18]] .= calcnut(u, grow, p[8], p[12], p[14], p[15], p[17], p[18])
	du[:, p[19]] .= calcnut(u, grow, p[9], p[13], p[14], p[16], p[17], p[19])


	return dz # allows function to give output outside of ODEProblem(...)
end
