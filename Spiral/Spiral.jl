### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 19bd1463-ceb6-46ff-871f-3f0117ebeac9
begin
	using Pkg
	Pkg.activate("..")
	Pkg.develop("GridShielding")
	using GridShielding
	using Plots
	using PlutoUI
	using Measures
	using Unzip
	using ProgressLogging
	using StaticArrays
	using Distances
	include("../Shared Code/FlatUI.jl")
end

# ╔═╡ ac05e283-3e2f-4f8d-9ccc-c6c96bf2e7ee
md"""
# "Spiral" Example
"""

# ╔═╡ d7a85669-fd51-4730-ac8b-9542b82e279a
md"""
# Preliminaries
"""

# ╔═╡ 8cf97456-6c01-4cd4-b3ac-97e1115620a8
TableOfContents()

# ╔═╡ 8b405589-2ba9-4046-ad9d-e2e3ccd21c84
md"""
# Model of the System
"""

# ╔═╡ d355823b-cc1c-40b9-96e0-a09cd79be7ff
# system definition
const A = [0.0 -1; 1 0]  # dynamics matrix of the standard oscillator

# ╔═╡ 52be2607-df01-458b-898e-901406747e3e
const δ = 0.05  # time step

# ╔═╡ ec77c85b-f27e-406b-98c3-5efbb8fa860a
const speed = 0.2

# ╔═╡ cb7fdd55-db5c-40ae-938c-0f3dd2b26b59
# simulation
r0 = 1.5  # initial radius

# ╔═╡ 61fa2c7c-4c61-4b6c-9db1-f35c011df04a
@enum Action stay_course move_out move_in

# ╔═╡ 4f87cc27-d210-4273-a4fc-89622bb8ef6e
SpiralState = SVector{2, Float64}

# ╔═╡ 9c3302dc-c01c-46f4-b9fd-42b4f3f7dd61
x0 = SpiralState([r0, 0])  # initial state

# ╔═╡ b047d22d-f223-446e-b58e-daf6a77b898d
struct SpiralTrace
	states::Vector{SpiralState}
	times::Vector{Float64}
	actions::Vector{Action}
end

# ╔═╡ 88ffd438-a6b3-4b11-82d8-af9ce61ed222
function apply(x, a::Action)
    θ = atan(x...)
    r = sqrt(x[1]^2 + x[2]^2)
	if a == move_out
    	r′ = (1 + (speed * δ)) * r
	elseif a == move_in
    	r′ = (1 - (speed * δ)) * r
	elseif a == stay_course
		r′ = r
	else 
		error("Unexpected a $a")
	end
    return [r′ * sin(θ), r′ * cos(θ)]
end

# ╔═╡ bb1ff345-58ba-498a-8035-88ef1db7d917
const expAδ = exp(A*δ)

# ╔═╡ 4937d230-e39f-4460-8071-a4cadc7e8b6f
# successor computation
function successor(x, act)::SpiralState
    return expAδ * apply(x, act)
end

# ╔═╡ bf560bb3-2fcf-4cf3-b6a4-5bacb7cbc832
function simulate_sequence(x0, policy, duration)::SpiralTrace
	result = SpiralTrace([x0], [0.], [])
	for t in δ:δ:duration
		x = result.states[end]
		a = policy(x)
		x′ = successor(x, a)
		push!(result.states, x′)
		push!(result.times, t)
		push!(result.actions, a)
	end
	result
end

# ╔═╡ b8b4af3d-9f66-4e3a-8845-0a45703bc82a
random(_...) = rand([stay_course, move_out, move_in])

# ╔═╡ cbe5bb1f-4a78-4185-9067-ddc7aa332be5
trace = simulate_sequence(x0, random, 6)

# ╔═╡ f7ef34f2-533c-43bd-bdb5-f407671facb1
md"""
## Rocks

We will only be checking if we are hitting rocks as part of our safety property.
"""

# ╔═╡ fd8f218e-42d5-4432-b8ca-2ede5f746417
struct Rock
	position::SpiralState
	radius::Float64
end

# ╔═╡ 4e5dc265-42a2-4874-a5cf-ae05d3ce4fdf
const rocks = [
	Rock([0, 0], 0.1),
	Rock([1.2, 1.2], 0.1), 
	Rock([-0.5, 0.7], 0.1),
	Rock([1.4, -0.7], 0.1),
	Rock([-1.2, -0.4], 0.1),
]

# ╔═╡ 61f009cc-a578-49e4-963a-d5b995c2ae39
md"""
## Trace Visualisation
"""

# ╔═╡ ffcb7ecc-4ea8-4e66-b39a-809810e8f105
function circle_shape(position, radius)
	x, y = position
	θ = LinRange(0, 2*π, 500)
	x .+ radius*sin.(θ), y .+ radius*cos.(θ)
end

# ╔═╡ 87962272-e8a3-4687-91dd-b7f6c23af97c
begin
	🎈1 = "This variable is used to ensure correct reactivity."
	
	@recipe function plot_rock(rock::Rock)
		circle_shape(rock.position, rock.radius)
	end
	
	@recipe function plot_rock(rocks::Vector{Rock})
		[circle_shape(rock.position, rock.radius) for rock in rocks]
	end

	
end

# ╔═╡ 905533f4-bf2b-4540-b4f0-d18f7fe039d1
let 
	i = 30
	[1/(i - 2j) for j in 1:i]
end

# ╔═╡ 7557a054-9b6d-4858-ba82-1b613c362b4b
md"""
# Shielding
"""

# ╔═╡ 7dddbfa0-35f6-4cec-90b6-37ae41881a77
begin
	function is_safe(x::SpiralState)
		r = √(x[1]^2 + x[2]^2)
		if !(r < 2)
			return false
		end

		for rock in rocks
			if euclidean(rock.position, x) <= rock.radius
				return false
			end
		end
		return true
	end

	function is_safe(bounds::Bounds)
		for x in SupportingPoints(3, bounds)
			if !is_safe(SpiralState(x))
				return false
			end
		end
		return true
	end
end

# ╔═╡ 9ff397e2-9aa7-433f-9a59-a4d3cb38a9bb
function plot_trace(trace::SpiralTrace, i=nothing)
	i = something(i, length(trace.states))
	path_colors = [is_safe(x) ? colors.PETER_RIVER : colors.POMEGRANATE 
		for x in trace.states[1:i]]

	alphas = [10/(i - j) for j in 1:i]
	
	plot(rocks, 
		label=nothing,
		seriestype=:shape,
		linewidth=0,
		color=colors.CONCRETE)
	
	plot!([(x[1], x[2]) for x in trace.states[1:i]],
		xlabel="x1",
		ylabel="x2",
		xlim=(-2.2, 2.2),
		ylim=(-2.2, 2.2),
		ratio=1,
		size=(300, 300),
		alpha=alphas,
		color=path_colors,
		marker=:circle,
		markersize=3,
		markerstrokewidth=0,
		label=nothing)
	
	scatter!([trace.states[i][1]], [trace.states[i][2]], 
		color=colors.NEPHRITIS,
		marker=:circle,
		markersize=4,
		markerstrokewidth=0,
		label=nothing)
end

# ╔═╡ 0a46d16c-e86d-443b-9001-6663a65ccccc
function animate_sequence(trace::SpiralTrace)
	🎈1
	🎥 = @animate for i in 1:length(trace.states)
		
		plot_trace(trace, i)
	end
	gif(🎥, fps=(1/δ)*2, show_msg=false)
end

# ╔═╡ f0b44048-4c9c-475f-87a4-9ff331885c32
animate_sequence(trace)

# ╔═╡ b67bbfee-633d-4d8b-846f-6a63ec778f93
plot_trace(trace)

# ╔═╡ e8a5370b-e18d-4957-b695-d6c50fec1182
begin
	any_action = actions_to_int(instances(Action))
	no_action = actions_to_int([])
	(;no_action, any_action)
end

# ╔═╡ 01c9a1be-4258-4aef-9855-db55874f6ae0
granularity = 0.01

# ╔═╡ 0f42e2d4-6446-401e-a062-b5aa893e9ac5
begin
	grid = Grid(granularity, Bounds(Float64[-2.1, -2.1], Float64[2.1, 2.1]))
	initialize!(grid, x -> is_safe(x) ? any_action : no_action)
	grid
end

# ╔═╡ 1c429eaf-f72f-4197-8242-12f41db29f81
size(grid)

# ╔═╡ 3c15b3dc-d3f0-4a0c-89cd-da2bbd6e4865
size(grid) |> prod

# ╔═╡ 74bebc78-c912-4ac8-9c94-1f85838e836b
begin
	shieldlabels = 	[
		"{$(join(int_to_actions(Action, i), ", "))}"
		for i in no_action:any_action]
	shieldcolors = 
		[colors.WET_ASPHALT, colors.AMETHYST, colors.SUNFLOWER, colors.PETER_RIVER, colors.CARROT, colors.EMERALD, colors.PUMPKIN, colors.CLOUDS]
	(;shieldcolors, shieldlabels)
end

# ╔═╡ 332ac25d-b07f-40ee-ab7d-e6cdfb2d075f
begin
	draw(grid;
		aspectratio=:equal,
		legend=:outerright,
		size=(800, 600),
		clim=(no_action, any_action),
		colors=shieldcolors, 
		color_labels=shieldlabels,)

	#draw_barbaric_transition!(simulation_model, partition, action)
	
	
	plot!(rocks, 
		seriestype=:shape,
		color=colors.CONCRETE,
		linewidth=0,
		label=nothing)
	
end

# ╔═╡ 5084de0e-e5e8-435b-9264-d858362957fb
action = stay_course

# ╔═╡ 7d5ba76c-e4cb-4c12-928e-2b9c30435288
x = [1, 1.1]

# ╔═╡ 1ce94535-be9e-4f07-b0f3-062509540516
partition = box(grid, x)

# ╔═╡ 0d0aac2a-15ae-422a-b33f-ef2810920d15
simulation_model = SimulationModel((x, a, _) -> successor(SpiralState(x), a), Bounds([], []), 3, 1)

# ╔═╡ 1c586fc6-1ec6-4f23-910b-4a9efb47e9fa
reachability_function = get_barbaric_reachability_function(simulation_model)

# ╔═╡ d652b956-7473-4042-91b1-b508067e9a62
md"""
## Mainmatter
"""

# ╔═╡ d5d84407-294d-4247-b00f-291530edf099
@bind max_steps NumberField(1:1000, default=1000)

# ╔═╡ af490584-bd8c-4e03-a64b-32bab91afc33
reachability_function_precomputed = get_transitions(reachability_function, Action, grid);

# ╔═╡ 398ce525-47fc-494a-ad4c-312b6c49c566
shield, max_steps_reached = make_shield(reachability_function_precomputed, Action, grid; max_steps)

# ╔═╡ 9ea2e239-00f8-4517-aeac-566415a5fa8a
begin
	draw(shield;
			xlabel="x1",
			ylabel="x2",
			aspectratio=:equal,
			legend=:outerright,
			size=(800, 600),
			clim=(no_action, any_action),
			colors=shieldcolors, 
			color_labels=shieldlabels,)
	
		
		plot!(rocks, 
			seriestype=:shape,
			color=colors.CONCRETE,
			linewidth=0,
			label=nothing)
		
end

# ╔═╡ 4bf2b834-15a3-4d1a-8314-8b851a9ca33b
if max_steps_reached
md"""
!!! warning "Shield not done"
	Increase `max_steps`
"""
else
md"""
!!! success "Shield finished 👍"
	Fixed-point iteration completed.
"""
end

# ╔═╡ 8bb0f2ad-5802-492d-b209-158d35b66a18
is_safe(x0)

# ╔═╡ 7dadc8b9-c52b-4578-8da4-36862f4c59f0
function apply_shield(shield::Grid, policy)
    return (s) -> begin
		a = policy(s)
		if s ∉ shield
			return a
		end
        allowed = int_to_actions(Action, get_value(box(shield, s)))
        if a ∈ allowed
            return a
        elseif length(allowed) > 0
			a′ = rand(allowed)
            return a′
        else
            return a
        end
    end
end

# ╔═╡ abc8c87f-df92-4320-91e9-7f1d5ad4462f
md"""
## Example trace
"""

# ╔═╡ b538485b-0c21-49c7-b9c0-bd7f13e2da3b
shielded_random = apply_shield(shield, random)

# ╔═╡ f161fcf5-83f2-4a44-87a5-58d0a052d01b
shielded_trace = simulate_sequence(x0, shielded_random, 60)

# ╔═╡ 14013a0c-9a27-44c5-b9ac-5d844fd3fe30
animate_sequence(shielded_trace)

# ╔═╡ 531281a1-6f2d-49fc-b608-cb1f8fc966ae
plot_trace(shielded_trace)

# ╔═╡ ba2fadc7-7677-4a50-a07b-5a542beb5b8a
md"""
# Shield with Altered State Space

I'm trying a different code structure in this notebook, letting the shields for the  regular state-space and the altered state space exist together. As such, there are a lot of variables and functions that are copied over. **These variables are pre-fixed with `a_` when they are related to the altered state-space.** 
"""

# ╔═╡ eae9e329-e241-46b3-aaf3-378d3067cd5a
md"""
## Transformation Functions
"""

# ╔═╡ 096bd467-1eb3-44f7-9b39-fc6490cedd57
AlteredState = SVector{2, Float64}

# ╔═╡ 508f0e5e-b668-400b-95d3-c484477fc855
function f(x)::AlteredState
	θ = atan(x...)
	r = euclidean(zeros(2), x)
	θ, r
end

# ╔═╡ 6c07703e-e945-4b7c-8db6-8790bc9fe9a1
f(x0)

# ╔═╡ e2cd27b9-79bf-400f-a8b5-3dc895d98ff1
function f⁻¹(x)::SpiralState
	θ, r = x
	r*sin(θ), r*cos(θ)
end

# ╔═╡ a8c248b4-e1b5-45f8-9ed3-1d81239f063f
begin
	function a_is_safe(bounds::Bounds)
		for x in SupportingPoints(6, bounds)
			if !is_safe(f⁻¹(x))
				return false
			end
		end
		return true
	end
end

# ╔═╡ 87df944d-74bb-443d-aa98-81209817b5f1
f⁻¹(f(x0)) ≈ x0

# ╔═╡ ca22b676-f646-41b1-85fe-c8f469f8fc06
begin
	a_grid = Grid([0.1, 0.005], Bounds([-pi - 0.1, 0.], [pi + 0.1, 3.5]))
	
	initialize!(a_grid, bounds -> a_is_safe(bounds) ? any_action : no_action)
	a_grid
end

# ╔═╡ f73c8ecc-d145-4b55-9ac1-5127633d50cb
a_is_safe.([Bounds(box(a_grid, f(rock.position))) for rock in rocks])

# ╔═╡ 172f1142-89d2-4646-881d-fb3d064a2c68
a_simulation_model = SimulationModel((x, a, _) -> f(successor(SpiralState(f⁻¹(x)), a)), Bounds([], []), 3, 1)

# ╔═╡ aa73fe1b-3c3e-45dc-8be4-332ebc529851
a_reachability_function = get_barbaric_reachability_function(a_simulation_model)

# ╔═╡ f5143b1b-cb09-4df0-9871-f975ba3b5b97
md"""
## Mainmatter
"""

# ╔═╡ 2607e37b-8e9e-4658-b99a-72ecf1643321
a_reachability_function_precomputed = get_transitions(a_reachability_function, Action, a_grid);

# ╔═╡ 22501f27-54f7-452c-b253-9d4612667c54
@bind a_max_steps NumberField(1:1000, default=1000)

# ╔═╡ 4459115d-27b5-47ae-b53e-aee8d4c6d6cf
a_shield, a_max_steps_reached = make_shield(a_reachability_function_precomputed, Action, a_grid; max_steps=a_max_steps)

# ╔═╡ 4663c31a-bede-45af-876d-8650f2a5125a
size(a_grid)

# ╔═╡ cc8226c6-1a47-4050-af50-5c8ca7f7a3a9
length(grid), length(a_grid)

# ╔═╡ 351e2fb7-70d9-4c5e-81d5-6ba651c490e8
draw(a_shield;
		xlabel="θ",
		ylabel="r",
		aspectratio=:equal,
		legend=:outerright,
		size=(800, 600),
		clim=(no_action, any_action),
		colors=shieldcolors, 
		color_labels=shieldlabels,)

# ╔═╡ 98822937-709b-4654-b3dc-e23342d28f0f
if a_max_steps_reached
md"""
!!! warning "Shield not done"
	Increase `max_steps`
"""
else
md"""
!!! success "Shield finished 👍"
	Enjoy your new shield
"""
end

# ╔═╡ e1a91b83-07d3-48f0-9c25-e27f3ce0258f
function draw_function(policy::Function, x_min, x_max, y_min, y_max, G; plotargs...)
	size_x, size_y = Int((x_max - x_min)/G), Int((y_max - y_min)/G)
	matrix = Matrix(undef, size_x, size_y)
	for i in 1:size_x
		for j in 1:size_y
			x, y = i*G - G + x_min, j*G - G + y_min

			matrix[i, j] = policy([x, y])
		end
	end
	x_tics = G+x_min:G:x_max
	y_tics = G+y_min:G:y_max
	middle_x, middle_y = [(x_max - x_min)/2 + x_min], [(y_max - y_min)/2 + y_min]
	plot(;plotargs...)
	heatmap!(x_tics, y_tics, transpose(matrix);
			plotargs...)
end

# ╔═╡ ee3c081b-9bc1-4dfc-b056-620eb305b4cd
let
	l, u = grid.bounds.lower, grid.bounds.upper
	
	draw_function(x -> get_value(box(a_shield, f(x))), l[1], u[1], l[2], u[2], 0.01,
			color=shieldcolors,
			xlabel="x1",
			ylabel="x2",
			ratio=1,
			colorbar=nothing)
	
	plot!(rocks, 
		seriestype=:shape,
		color=colors.CONCRETE,
		linewidth=0,
		label=nothing)
	
end

# ╔═╡ a5b2aa32-8267-45da-8664-98ea3afe4671
md"""
## Example trace
"""

# ╔═╡ c8b232d3-fe05-453e-9bed-837c83a81a6e
function a_apply_shield(shield::Grid, policy)
    return (s) -> begin
		a = policy(s)
		if f(s) ∉ shield
			return a
		end
        allowed = int_to_actions(Action, get_value(box(shield, f(s))))
        if a ∈ allowed
            return a
        elseif length(allowed) > 0
			a′ = rand(allowed)
            return a′
        else
            return a
        end
    end
end

# ╔═╡ aa8ebc72-8300-405b-b6da-704e39f19506
a_shielded_random = a_apply_shield(a_shield, random)

# ╔═╡ fe72c613-014f-4c07-bee5-45efeb2bb770
a_shielded_trace = simulate_sequence(x0, a_shielded_random, 100)

# ╔═╡ 13646ad3-c39b-46d0-bddf-c16206dbb7eb
animate_sequence(a_shielded_trace)

# ╔═╡ 87356f08-c963-4928-8a0f-48ca9b8d82b2
plot_trace(a_shielded_trace)

# ╔═╡ d2427d85-7e20-4d2e-865c-2c582f65fe87
function a_plot_trace!(trace::SpiralTrace, i=nothing)
	i = something(i, length(trace.states))
	path_colors = [is_safe(x) ? colors.PETER_RIVER : colors.POMEGRANATE 
		for x in trace.states[1:i]]

	alphas = [10/(i - j) for j in 1:i]
	
	plot!([(f(x)[1], f(x)[2]) for x in trace.states[1:i]],
		xlabel="θ",
		ylabel="r",
		ratio=1,
		size=(300, 300),
		alpha=alphas,
		color=path_colors,
		marker=:circle,
		markersize=3,
		markerstrokewidth=0,
		label=nothing)
	
	scatter!([f(trace.states[i])[1]], [f(trace.states[i])[2]], 
		color=colors.NEPHRITIS,
		marker=:circle,
		markersize=4,
		markerstrokewidth=0,
		label=nothing)
end

# ╔═╡ 9ce546c6-361d-4186-b585-2534d38614b6
function a_animate_sequence(trace::SpiralTrace)
	🎈1
	background = draw(a_shield;
				xlabel="θ",
				ylabel="r",
				clim=(no_action, any_action),
				colors=shieldcolors, 
				color_labels=shieldlabels,
				legend=:none)
		
	🎥 = @animate for i in 1:length(trace.states)
		plot(background)
		a_plot_trace!(trace, i)
	end
	gif(🎥, fps=(1/δ)*2, show_msg=false)
end

# ╔═╡ 875fe159-8fb4-437b-a327-3e2129957484
begin
	
	a_animate_sequence(a_shielded_trace)

end

# ╔═╡ Cell order:
# ╟─ac05e283-3e2f-4f8d-9ccc-c6c96bf2e7ee
# ╟─d7a85669-fd51-4730-ac8b-9542b82e279a
# ╠═19bd1463-ceb6-46ff-871f-3f0117ebeac9
# ╠═8cf97456-6c01-4cd4-b3ac-97e1115620a8
# ╟─8b405589-2ba9-4046-ad9d-e2e3ccd21c84
# ╠═d355823b-cc1c-40b9-96e0-a09cd79be7ff
# ╠═52be2607-df01-458b-898e-901406747e3e
# ╠═ec77c85b-f27e-406b-98c3-5efbb8fa860a
# ╠═cb7fdd55-db5c-40ae-938c-0f3dd2b26b59
# ╠═9c3302dc-c01c-46f4-b9fd-42b4f3f7dd61
# ╠═61fa2c7c-4c61-4b6c-9db1-f35c011df04a
# ╠═4f87cc27-d210-4273-a4fc-89622bb8ef6e
# ╠═b047d22d-f223-446e-b58e-daf6a77b898d
# ╠═88ffd438-a6b3-4b11-82d8-af9ce61ed222
# ╠═bb1ff345-58ba-498a-8035-88ef1db7d917
# ╠═4937d230-e39f-4460-8071-a4cadc7e8b6f
# ╠═bf560bb3-2fcf-4cf3-b6a4-5bacb7cbc832
# ╠═b8b4af3d-9f66-4e3a-8845-0a45703bc82a
# ╠═cbe5bb1f-4a78-4185-9067-ddc7aa332be5
# ╟─f7ef34f2-533c-43bd-bdb5-f407671facb1
# ╠═fd8f218e-42d5-4432-b8ca-2ede5f746417
# ╠═4e5dc265-42a2-4874-a5cf-ae05d3ce4fdf
# ╟─61f009cc-a578-49e4-963a-d5b995c2ae39
# ╠═ffcb7ecc-4ea8-4e66-b39a-809810e8f105
# ╠═87962272-e8a3-4687-91dd-b7f6c23af97c
# ╠═9ff397e2-9aa7-433f-9a59-a4d3cb38a9bb
# ╠═905533f4-bf2b-4540-b4f0-d18f7fe039d1
# ╠═0a46d16c-e86d-443b-9001-6663a65ccccc
# ╠═f0b44048-4c9c-475f-87a4-9ff331885c32
# ╠═b67bbfee-633d-4d8b-846f-6a63ec778f93
# ╟─7557a054-9b6d-4858-ba82-1b613c362b4b
# ╠═7dddbfa0-35f6-4cec-90b6-37ae41881a77
# ╠═e8a5370b-e18d-4957-b695-d6c50fec1182
# ╠═01c9a1be-4258-4aef-9855-db55874f6ae0
# ╠═0f42e2d4-6446-401e-a062-b5aa893e9ac5
# ╠═1c429eaf-f72f-4197-8242-12f41db29f81
# ╠═3c15b3dc-d3f0-4a0c-89cd-da2bbd6e4865
# ╠═74bebc78-c912-4ac8-9c94-1f85838e836b
# ╟─332ac25d-b07f-40ee-ab7d-e6cdfb2d075f
# ╠═5084de0e-e5e8-435b-9264-d858362957fb
# ╠═1ce94535-be9e-4f07-b0f3-062509540516
# ╠═7d5ba76c-e4cb-4c12-928e-2b9c30435288
# ╠═0d0aac2a-15ae-422a-b33f-ef2810920d15
# ╠═1c586fc6-1ec6-4f23-910b-4a9efb47e9fa
# ╟─d652b956-7473-4042-91b1-b508067e9a62
# ╠═d5d84407-294d-4247-b00f-291530edf099
# ╠═af490584-bd8c-4e03-a64b-32bab91afc33
# ╠═398ce525-47fc-494a-ad4c-312b6c49c566
# ╟─9ea2e239-00f8-4517-aeac-566415a5fa8a
# ╟─4bf2b834-15a3-4d1a-8314-8b851a9ca33b
# ╠═8bb0f2ad-5802-492d-b209-158d35b66a18
# ╠═7dadc8b9-c52b-4578-8da4-36862f4c59f0
# ╟─abc8c87f-df92-4320-91e9-7f1d5ad4462f
# ╠═b538485b-0c21-49c7-b9c0-bd7f13e2da3b
# ╠═f161fcf5-83f2-4a44-87a5-58d0a052d01b
# ╠═14013a0c-9a27-44c5-b9ac-5d844fd3fe30
# ╠═531281a1-6f2d-49fc-b608-cb1f8fc966ae
# ╟─ba2fadc7-7677-4a50-a07b-5a542beb5b8a
# ╟─eae9e329-e241-46b3-aaf3-378d3067cd5a
# ╠═096bd467-1eb3-44f7-9b39-fc6490cedd57
# ╠═508f0e5e-b668-400b-95d3-c484477fc855
# ╠═6c07703e-e945-4b7c-8db6-8790bc9fe9a1
# ╠═e2cd27b9-79bf-400f-a8b5-3dc895d98ff1
# ╠═a8c248b4-e1b5-45f8-9ed3-1d81239f063f
# ╠═f73c8ecc-d145-4b55-9ac1-5127633d50cb
# ╠═87df944d-74bb-443d-aa98-81209817b5f1
# ╠═ca22b676-f646-41b1-85fe-c8f469f8fc06
# ╠═172f1142-89d2-4646-881d-fb3d064a2c68
# ╠═aa73fe1b-3c3e-45dc-8be4-332ebc529851
# ╟─f5143b1b-cb09-4df0-9871-f975ba3b5b97
# ╠═2607e37b-8e9e-4658-b99a-72ecf1643321
# ╠═4459115d-27b5-47ae-b53e-aee8d4c6d6cf
# ╠═22501f27-54f7-452c-b253-9d4612667c54
# ╠═4663c31a-bede-45af-876d-8650f2a5125a
# ╠═cc8226c6-1a47-4050-af50-5c8ca7f7a3a9
# ╟─351e2fb7-70d9-4c5e-81d5-6ba651c490e8
# ╟─98822937-709b-4654-b3dc-e23342d28f0f
# ╠═e1a91b83-07d3-48f0-9c25-e27f3ce0258f
# ╠═ee3c081b-9bc1-4dfc-b056-620eb305b4cd
# ╟─a5b2aa32-8267-45da-8664-98ea3afe4671
# ╠═c8b232d3-fe05-453e-9bed-837c83a81a6e
# ╠═aa8ebc72-8300-405b-b6da-704e39f19506
# ╠═fe72c613-014f-4c07-bee5-45efeb2bb770
# ╠═13646ad3-c39b-46d0-bddf-c16206dbb7eb
# ╠═87356f08-c963-4928-8a0f-48ca9b8d82b2
# ╠═d2427d85-7e20-4d2e-865c-2c582f65fe87
# ╠═9ce546c6-361d-4186-b585-2534d38614b6
# ╠═875fe159-8fb4-437b-a327-3e2129957484
