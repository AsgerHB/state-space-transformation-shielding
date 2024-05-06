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
	Pkg.develop("UppaalTraceParser")
	using GridShielding
	using UppaalTraceParser
	using Plots
	using PlutoUI
	using Measures
	using Unzip
	using ProgressLogging
	using StaticArrays
	using PlutoLinks
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

# ╔═╡ a115e214-31a5-4a62-a798-07ebbc67caa1
function multiline(str)
	HTML("""
	<pre style='max-height:30em; margin:8pt 0 8pt 0; overflow-y:scroll'>
	$str
	</pre>
	""")
end

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

# ╔═╡ acd8f481-ba99-42cf-b95a-8c175b996e1d
1/δ

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

# ╔═╡ b8b4af3d-9f66-4e3a-8845-0a45703bc82a
random(_...) = rand([move_in, stay_course, move_out])

# ╔═╡ cbe5bb1f-4a78-4185-9067-ddc7aa332be5
trace = simulate_sequence(x0, random, 6)

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
function plot_trace(trace::SpiralTrace, i=nothing; background=plot())
	i = something(i, length(trace.states))
	path_colors = [is_safe(x) ? colors.PETER_RIVER : colors.POMEGRANATE 
		for x in trace.states[1:i]]

	alphas = [10/(i - j) for j in 1:i]

	plot(background)
	
	plot!(rocks, 
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
function animate_sequence(trace::SpiralTrace; background=plot())
	🎈1
	🎥 = @animate for i in 1:length(trace.states)
		
		plot_trace(trace, i; background)
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

# ╔═╡ 5084de0e-e5e8-435b-9264-d858362957fb
action = stay_course

# ╔═╡ 7d5ba76c-e4cb-4c12-928e-2b9c30435288
x = [1, 1.1]

# ╔═╡ 1ce94535-be9e-4f07-b0f3-062509540516
partition = box(grid, x)

# ╔═╡ 44db0342-b1f1-4e0c-a422-b95bf4c78dbe
samples_per_axis = 6

# ╔═╡ 0d0aac2a-15ae-422a-b33f-ef2810920d15
simulation_model = SimulationModel((x, a, _) -> successor(SpiralState(x), a), Bounds([], []), samples_per_axis, 1)

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

# ╔═╡ 74bebc78-c912-4ac8-9c94-1f85838e836b
begin
	shieldlabels = 	[
		"{$(join(int_to_actions(Action, i), ", "))}"
		for i in no_action:any_action]
	shieldcolors = 
		[	colors.WET_ASPHALT, 
			colors.CONCRETE, 
			colors.POMEGRANATE, 
			colors.ORANGE, 
			colors.PETER_RIVER,
			colors.WISTERIA,
			colors.SUNFLOWER,
			colors.CLOUDS]
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
	
	scatter!([],
		markershape=:o,
		markerstrokewidth=0,
		color=colors.CONCRETE, 
		label="Rock")
end

# ╔═╡ 9ea2e239-00f8-4517-aeac-566415a5fa8a
begin
	p1 = draw(shield;
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
	
	scatter!([],
		markershape=:o,
		markerstrokewidth=0,
		color=colors.CONCRETE, 
		label="Rock")
		
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
shielded_trace = simulate_sequence(x0, shielded_random, 25)

# ╔═╡ 14013a0c-9a27-44c5-b9ac-5d844fd3fe30
animate_sequence(shielded_trace)

# ╔═╡ 2f658797-0581-40c4-9054-4ed968ad143b
let
	📈 = plot(p1, legend=nothing)
	animate_sequence(shielded_trace, background=📈)
end

# ╔═╡ 531281a1-6f2d-49fc-b608-cb1f8fc966ae
plot_trace(shielded_trace)

# ╔═╡ 7d1b7c12-9bf3-4b7a-891a-7139acefffd3
md"""
## Make UPPAAL-friendly
Same deal as with the bouncing-ball: Since UPPAAL's simulation is not accurate, the radius can decrease by an infinitesimal amount, even if the action is `stay_course`. 

This is not part of the model, this is just UPPAAL being inaccurate. To address it, the shield will be made robust with an ad-hoc method, redefining the values of some squares at the border that are supposed to be unreachable.
"""

# ╔═╡ 6e5eee38-2267-45a2-9a0c-a2663ca207a8
"""
	adjacent_to(partition::Partition, value)

Returns true if `partition` is adjacent to (including corners) another partition with value `value`.
"""
function adjacent_to(partition::Partition, value)
	if length(partition.indices) != 2 
		error("Not implemented: Only 2d grids supported.")
	end
	grid_size = size(partition.grid)
	for offset_x in [-1, 0, 1]
		for offset_y in [-1, 0, 1]
			if offset_x == offset_y == 0
				continue
			end
			indices′ = partition.indices .+ [offset_x, offset_y]
			if indices′[1] < 1 || indices′[2] < 1
				continue
			elseif indices′[1] >  grid_size[1] || indices′[2] >  grid_size[2]
				continue
			end
			partition′ = Partition(partition.grid, indices′)
			if get_value(partition′) == value
				return true
			end
		end
	end
	return false
end

# ╔═╡ 3ba16e83-5a6f-4bff-a2ad-dd20785f33dc
let
	result = Partition[]
	must_move_in = actions_to_int([move_in])
	for partition in shield
		if get_value(partition) != no_action continue end
		if adjacent_to(partition, must_move_in)
			push!(result, partition)
		end
	end
	plot(p1)
	
	scatter!([Tuple(Bounds(p).lower) for p in result], 
		title="Testing adjacent_to",
		legend=nothing, 
		markersize=3,
		size=(300, 300))
end

# ╔═╡ 744dcdcf-1db6-4544-b4f0-d3d49149b1d7
function make_uppaal_friendly!(shield::Grid)
	must_move_out = actions_to_int([move_out])
	must_move_in = actions_to_int([move_in])
	stay_or_move_out = actions_to_int([stay_course, move_out])
	stay_or_move_in = actions_to_int([stay_course, move_in])
	updates = []
	for partition in shield
		if get_value(partition) != no_action continue end
		if adjacent_to(partition, must_move_out)
			push!(updates, (partition, must_move_out))
			continue
		end
		if adjacent_to(partition, must_move_in)
			push!(updates, (partition, must_move_in))
			continue
		end
		if adjacent_to(partition, stay_or_move_out)
			push!(updates, (partition, stay_or_move_out))
			continue
		end
		if adjacent_to(partition, stay_or_move_in)
			push!(updates, (partition, stay_or_move_in))
			continue
		end
	end
	for update in updates
		set_value!(update...)
	end
end

# ╔═╡ d5101b8b-e65d-4a61-a616-803c9415bdcc
@bind make_uppaal_friendly_button CounterButton("Make UPPAAL-friendly")

# ╔═╡ d9506fdc-4bd4-40f6-9686-0403628eecf7
made_uppaal_friendly = if make_uppaal_friendly_button == 1
	make_uppaal_friendly!(shield)
	"done."
elseif make_uppaal_friendly_button > 1
	"already done."
end

# ╔═╡ 19e227b0-df36-4835-9a4d-514490c67062
md"""
## Check Safety
"""

# ╔═╡ 60a01054-37be-4bcd-98fb-023350073ecf
runs = 1000

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
	a_grid = Grid([0.1, 0.005], Bounds([-pi - 0.1, 0.], [pi + 0.1, 2.1]))
	
	initialize!(a_grid, bounds -> a_is_safe(bounds) ? any_action : no_action)
	a_grid
end

# ╔═╡ f73c8ecc-d145-4b55-9ac1-5127633d50cb
a_is_safe.([Bounds(box(a_grid, f(rock.position))) for rock in rocks])

# ╔═╡ 172f1142-89d2-4646-881d-fb3d064a2c68
a_simulation_model = SimulationModel((x, a, _) -> clamp(Vector(f(successor(SpiralState(f⁻¹(x)), a))), a_grid.bounds), Bounds([], []), samples_per_axis, 1)

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
begin
	p2 = draw(a_shield;
		xlabel="θ",
		ylabel="r",
		aspectratio=:equal,
		legend=:outerright,
		size=(800, 600),
		clim=(no_action, any_action),
		colors=shieldcolors, 
		color_labels=shieldlabels,)

	a_rocks = [Rock(f(rock.position), rock.radius) for rock in rocks
		if rock.position != [0, 0]]

	plot!(a_rocks, 
		seriestype=:shape,
		color=colors.CONCRETE,
		linewidth=0,
		label=nothing)
	
	scatter!([],
		markershape=:o,
		markerstrokewidth=0,
		color=colors.CONCRETE, 
		label="Rock")
end

# ╔═╡ c5570c5b-5d4c-44bf-9dde-2d5e9a8e8ed4
begin
	plot(p2, projection=:polar)
end

# ╔═╡ 6c3550be-a891-4019-a97d-d4f28e3b1912
get_value(box(a_shield, 1.5, 1.5))

# ╔═╡ 98822937-709b-4654-b3dc-e23342d28f0f
if a_max_steps_reached
md"""
!!! warning "Shield not done"
	Increase `max_steps`
"""
else
md"""
!!! success "Shield finished 👍"
	Fixed point reached.
"""
end

# ╔═╡ e1a91b83-07d3-48f0-9c25-e27f3ce0258f
function draw_function(policy::Function, x_min, x_max, y_min, y_max, G; 
	colors=[],
	color_labels=[],
	plotargs...)
	
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

	# Show labels
	if length(color_labels) > 0
		if length(color_labels) != length(colors)
			throw(ArgumentError("Length of argument color_labels does not match  number of colors."))
		end
		for (color, label) in zip(colors, color_labels)
			# Apparently shapes are added to the legend even if the list is empty
		    plot!(Float64[], Float64[], seriestype=:shape, 
		        label=label, color=color)
		end
	end
end

# ╔═╡ ee3c081b-9bc1-4dfc-b056-620eb305b4cd
p3 = let
	l, u = grid.bounds.lower, grid.bounds.upper
	
	📈  = draw_function(x -> get_value(box(a_shield, 
					clamp(Vector(f(x)), a_shield.bounds) )), 
		
			l[1], u[1], l[2], u[2], 0.005;
			color=shieldcolors,
			xticks=-2:1:2,
			yticks=-2:1:2,
			xlabel="x1",
			ylabel="x2",
			ratio=1,
			size=(800, 600),
			legend=:outerright,
			colors=shieldcolors, 
			color_labels=shieldlabels,
			colorbar=nothing)
	
	plot!(rocks, 
		seriestype=:shape,
		color=colors.CONCRETE,
		linewidth=0,
		label=nothing)

	scatter!([],
		markershape=:o,
		markerstrokewidth=0,
		color=colors.CONCRETE, 
		label="Rock")
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
a_shielded_trace = simulate_sequence(x0, a_shielded_random, 10)

# ╔═╡ 13646ad3-c39b-46d0-bddf-c16206dbb7eb
let
	📈 = plot(p3, legend=nothing)
	animate_sequence(a_shielded_trace, background=📈)
end

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

# ╔═╡ 319cebe6-33ce-4a6f-be54-e9fb6dd059b0
md"""
## Make UPPAAL-friendly
Same deal as above.
"""

# ╔═╡ c8867810-099d-4dfc-909a-a5e0a0773bb8
@bind a_make_uppaal_friendly_button CounterButton("Make UPPAAL-friendly")

# ╔═╡ a7b0f796-6532-4a87-979c-f9c8759813f6
a_made_uppaal_friendly = if a_make_uppaal_friendly_button == 1
	make_uppaal_friendly!(a_shield)
	"done."
elseif a_make_uppaal_friendly_button > 1
	"already done."
end

# ╔═╡ ddf10c36-f525-4384-b5ac-731abc9d2b1e
md"""
## Check Safety
"""

# ╔═╡ cda1ca0f-b02b-42cf-85fb-20175578396b
a_runs = 1000

# ╔═╡ 3061a25a-bae8-4763-b84b-59e3536fd63c
function check_safety(policy, duration; runs=1000)
	deaths = 0
	example_trace = nothing
	@progress for run in 1:runs
		trace = simulate_sequence((r0, 0), policy, duration)
		for s in trace.states
			if !is_safe(s)
				deaths += 1
				example_trace = trace
				break
			end
		end
		example_trace = something(example_trace, trace)
	end
	deaths, example_trace
end

# ╔═╡ 4c6ca97b-4871-4f4c-85e8-ee02719a2433
deaths, example_trace = check_safety(shielded_random, 100; runs)

# ╔═╡ c780a26f-d6c8-4573-b919-2500b9e4d8f4
let
	header = if deaths > 0
		"""!!! danger "Shield unsafe"

		"""
	else
		"""!!! success "Shield safe"

		"""
	end

	Markdown.parse("""$header
		Out of **$runs** runs, **$deaths** of them contained a safety violation.
	""")
end

# ╔═╡ 150a3ed4-3a97-4a75-b836-3c990d127f8f
a_deaths, a_example_trace = check_safety(a_shielded_random, 100; runs=a_runs)

# ╔═╡ c9412832-1653-449b-bff3-99147ec7f3a6
let
	header = if a_deaths > 0
		"""!!! danger "Shield unsafe"

		"""
	else
		"""!!! success "Shield safe"

		"""
	end

	Markdown.parse("""$header
		Out of **$runs** runs, **$a_deaths** of them contained a safety violation.
	""")
end

# ╔═╡ 4ec52885-8315-403d-8eeb-d49a4eac8f16


# ╔═╡ f47335e6-999d-449d-bf95-1f184b898042
md"""
# Exporting the Shield
"""

# ╔═╡ e514a65b-eb40-4193-92cd-406273e43d9b
@bind target_dir TextField(95, default=mktempdir())

# ╔═╡ 3b0d8736-713c-4fde-85db-7cfb239c45fe
target_dir; @bind open_folder_button CounterButton("Open Folder")

# ╔═╡ 551adaa4-3c3c-47a1-a263-db3fefcaf4f0
if open_folder_button > 0
	run(`nautilus $target_dir`, wait=false)
end; "This cell opens `$target_dir` in nautilus"

# ╔═╡ 04d7c2d6-fb7c-4329-b349-0d8563c8cde8
md"""
### Export as serialized julia-tuple

Easy export and import between julia code.
"""

# ╔═╡ ab11ba27-c11f-41cb-8266-9b8c9e39008a
let
	made_uppaal_friendly # Reactivity
	filename = "Spiral - Standard State Space.shield"
	
	robust_grid_serialization(joinpath(target_dir, filename), shield)
	
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ c58c5add-cca1-4102-a71a-af82e3115993
let
	a_made_uppaal_friendly # Reactivity
	filename = "Spiral - Altered State Space.shield"
	
	robust_grid_serialization(joinpath(target_dir, filename), a_shield)
	
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ fe28308c-829d-4b9d-a821-9f2aba4204aa
md"""
### Export as a function in a shared-object library

Use this library to access the shield from C and C++ code.

The shield is compiled into a shared-object binary, which exports the function `int get_value(double v, double p)`. It takes the state-variables as input and returns the bit-encoded list of allowed actions. (See `int_to_actions`.)
"""

# ╔═╡ 35079353-ac6d-4c09-b25f-cdcaf526c914
begin
	made_uppaal_friendly # Reactivity
	shield_so = "spiral_shield_standard_state_space.so"
	shield_so = joinpath(target_dir, shield_so)
	
	get_libshield(shield; destination=shield_so, force=true)
	
	"Exported `'$shield_so'`." |> Markdown.parse
end

# ╔═╡ abefb2a8-19fd-4c92-a81a-b0513572b756
begin
	a_made_uppaal_friendly # Reactivity
	a_shield_so = "spiral_shield_altered_state_space.so"
	a_shield_so = joinpath(target_dir, a_shield_so)
	
	get_libshield(a_shield; destination=a_shield_so, force=true)
	
	"Exported `'$a_shield_so'`." |> Markdown.parse
end

# ╔═╡ 206bb5c4-f1fd-4093-b7d4-726a39a2bb16
const foo = shield_so

# ╔═╡ c17af70b-2bc5-4949-8133-499c8ea3c0e2
c_get_value = (x1, x2) -> @ccall foo.get_value(x1::Cdouble, x2::Cdouble)::Cint

# ╔═╡ e5f19917-868f-41f0-8873-6eef2d40bfaf
md"""
### Export to Numpy

Exports a zip-file containing a serialized numpy-array along with a JSON file with details on how to read it.
"""

# ╔═╡ 46e15850-9632-4232-9d8f-05f3a9cb4832
let
	made_uppaal_friendly # Reactivity

	meta_info = (;variables=["x1", "x2"], 
		actions=Action,
		env_id="Spiral")
	
	filename = "Spiral Shield - Standard State Space.zip"
	
	numpy_zip_file(shield, joinpath(target_dir, filename); meta_info...)
	
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ 852cb308-c800-42ee-9412-412fe9b9a05e
let
	a_made_uppaal_friendly # Reactivity
	meta_info = (;variables=["angle", "radius"], 
		actions=Action,
		env_id="Bouncing Ball")
	
	filename = "Spiral Shield - Altered State Space.zip"
	
	numpy_zip_file(a_shield, joinpath(target_dir, filename); meta_info...)
	
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ 641033dd-aa48-43fd-b81b-c3ad84b58327
md"""
# Sampling plot

The plot that explains how sampling is done.
"""

# ╔═╡ 80934c04-60dd-4fbc-9563-0954aa72eec2
let
	# Overwrite successor function with custom time step
	function apply(x, a::Action, δ)
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
	successor(x, a, δ) = (exp(A*δ))*apply(x, a, δ*2)

	action = move_in
	δ = 0.9
	
	a_grid = Grid(0.5, Bounds([-pi - 0.1, 0.], [pi + 0.1, 3.5]))
	a_partition = box(a_grid, 0, 1.5)
	a_samples = SupportingPoints(2, a_partition)

	samples = [Tuple(f⁻¹(s)) for s in a_samples]
	reached = [Tuple(successor(x, action, δ)) for x in samples]
	
	trajectories = [[Tuple(successor(x, action, δ′)) for δ′ in 0:0.01:δ] 
		for x in samples]

	a_reached = [Tuple(f(x)) for x in reached]

	a_reached_partitions = [box(a_grid, s) for s in a_reached] |> unique

	margin = 1
	
	xlim = (min([s[1] for s in [samples..., reached...]]...) - margin,
			max([s[1] for s in [samples..., reached...]]...) + margin)
	
	ylim = (min([s[2] for s in [samples..., reached...]]...) - margin,
			max([s[2] for s in [samples..., reached...]]...) + margin)

	a_margin = 1
	
	a_xlim = (min([s[1] for s in [a_samples..., a_reached...]]...) - a_margin,
			max([s[1] for s in [a_samples..., a_reached...]]...) + a_margin)
	
	a_ylim = (min([s[2] for s in [a_samples..., a_reached...]]...) - a_margin,
			max([s[2] for s in [a_samples..., a_reached...]]...) + a_margin)
	
	📈1 = draw(a_grid, 
		xlabel="θ",
		ylabel="r",
		xlim = a_xlim,
		ylim = a_ylim,
		legend=:outertop,
		show_grid=true,
		ratio=1)
	
	plot!(Bounds(a_partition), 
		color=colors.EMERALD, 
		linewidth=0,
		opacity=0.5,
		label="initial partition")

	for a_reached_partition in a_reached_partitions
		plot!(Bounds(a_reached_partition), 
			color=colors.PETER_RIVER, 
			linewidth=0,
			opacity=0.5,
			label=nothing) # label added later
	end
	
	scatter!(unzip(a_samples), 
		markersize=5,
		markerstrokewidth=0,
		label="initial samples",
		markercolor=colors.EMERALD)
	
	scatter!(unzip(a_reached), 
		markersize=5,
		markerstrokewidth=0,
		label="successors translated to altered state space",
		markercolor=colors.PETER_RIVER)

	plot!([], seriestype=:shape,
		color=colors.PETER_RIVER, 
		linewidth=0,
		opacity=0.5,
		label="partition marked as reached")
	
	📈2 = plot(;
		legend=:outertop,
		ratio=1,
		xlabel="x1",
		ylabel="x2",
		grid=false,
		xlim,
		ylim)
	
	scatter!(unzip(samples), 
		markersize=5,
		markerstrokewidth=0,
		label="samples in standard state space",
		markercolor=colors.EMERALD)
	
	
	scatter!(unzip(reached), 
		markersize=5,
		markerstrokewidth=0,
		label="successor states",
		markercolor=colors.PETER_RIVER)

	for trajectory in trajectories
		
		plot!(trajectory, 
			linecolor=colors.CONCRETE,
			linestyle=:dash,
			arrow=:cap,
			label=nothing)
	end
	
	
	plot(📈1, 📈2,
		size=(600, 300))
end

# ╔═╡ b65755fc-1623-4f67-9aec-54dcd3872ae0
md"""
# Run in UPPAAL
"""

# ╔═╡ 52ac166d-61f6-4a4f-b65d-7b8e1db22261
@bind use_altered_state_space_shield_in_uppaal CheckBox(default=true)

# ╔═╡ b015bc8a-bc0e-436e-8edd-52c6de7cb528
@bind uppaal_model TextField(80, joinpath(pwd(), "Spiral.xml"))

# ╔═╡ bfdd4c8f-b980-439b-b774-598a77e8ac08
if use_altered_state_space_shield_in_uppaal
	replacements = Dict(
		"SHIELD_IMPORT" => "import \"$a_shield_so\" ",
		"SHIELD_SIGNATURE" => "    int get_value(double angle, double radius);",
		"SHIELD_LOOKUP" => "    int value = get_value(angle, radius);",
	)
else
	replacements = Dict(
		"SHIELD_IMPORT" => "import \"$shield_so\" ",
		"SHIELD_SIGNATURE" => "    int get_value(double x1, double x2);",
		"SHIELD_LOOKUP" => "    int value = get_value(x1, x2);",
	)
end

# ╔═╡ 24730e1a-b691-4fb0-8558-12724fe07697


# ╔═╡ c2a9ad89-308c-4291-9efe-cb3a3cd36187
@bind query TextField((80, 4), "simulate[<=10;100] {x1, x2, collisions, unsafe_state_entered}")

# ╔═╡ 01bac781-c4b0-4378-8e42-8dcc6eb5ac6e
uppaal_output = run_model(uppaal_model, query, replacements, 
	working_dir=target_dir,
	discretization=0.001);

# ╔═╡ 53cee99d-cb92-4d43-98ac-8e86b04e298a
uppaal_output |> multiline

# ╔═╡ 2e94277a-f16d-47ab-b215-9bb17e0537df
uppaal_traces = parse_traces(uppaal_output, δ);

# ╔═╡ d629a401-e2b3-4020-9bdd-5f2101aba751
md"""
## Inspect trace (code no work)
All thsi code is bust and I'm not fixing it.
"""

# ╔═╡ 877c0c7b-563c-4b6b-8195-14577d04ef14
uppaal_example_trace′ = let
	trace = uppaal_traces[1]
	for i in 1:length(uppaal_traces)
		if uppaal_traces[i]["collisions"][end] > 0
			trace = uppaal_traces[i]
			break
		end
	end
	trace
end

# ╔═╡ bb4d99ba-9752-430c-9898-edf372c2ee41
if uppaal_example_trace′["collisions"][end] > 0
	"""
	!!! danger "Found an unsafe UPPAAL trace."
		Collisions: $(uppaal_example_trace′["collisions"][end])
	""" |> Markdown.parse
else
	"""
	!!! success "All $(length(uppaal_traces)) traces safe"
		Collisions: $(uppaal_example_trace′["collisions"][end])
	""" |> Markdown.parse
end

# ╔═╡ 4bea0ab7-3254-46fb-b48f-a04b5cf90819
function convert_uppaal_action(a::Number)::Action
	isapprox(a, 0.0; atol=1e-8) ? stay_course :
	isapprox(a, 1.0; atol=1e-8) ? move_out :
	isapprox(a, 2.0; atol=1e-8) ? move_in : error("unexpected action: $a")
end

# ╔═╡ 21cdf64b-ba00-447d-9b04-aac3a7a0269e
uppaal_example_trace = let
	x1 = uppaal_example_trace′["x1"]
	x2 = uppaal_example_trace′["x2"]
	actions = uppaal_example_trace′["action"]
	trace = SpiralTrace(
		collect(zip(x1, x2)), 
		collect(0:δ:100), 
		[convert_uppaal_action(a) for a in actions][1:end])

end;

# ╔═╡ 3e0a5ab6-7158-46fb-9684-2027fee6c340
@bind i NumberField(1:length(uppaal_example_trace′["unsafe_state_entered"]), 
	default=findfirst((>)(0), uppaal_example_trace′["unsafe_state_entered"]))

# ╔═╡ 620bb25c-2bf9-4fc1-9543-1254d8eeb3db
begin
	plot_trace(uppaal_example_trace, i, background=p3)
	plot!(legend=nothing)
end

# ╔═╡ 410f184a-2055-419a-a986-65b817bdb124
uppaal_example_trace.states[i]

# ╔═╡ 5e29d540-3ad9-4d89-b829-ae6fe61bdf36
uppaal_example_trace′["collisions"][i]

# ╔═╡ 12844f71-adeb-41d2-afd3-b6bc6c53fa92
i_action = uppaal_example_trace.actions[i]

# ╔═╡ feabe121-3aae-492d-8f68-245256a7a0d8
i_partition = box(a_shield, f(uppaal_example_trace.states[i]))

# ╔═╡ 75110d13-3374-4497-89c5-f82c9d726782
(uppaal_example_trace.states[i]), f(uppaal_example_trace.states[i])

# ╔═╡ 24bceb53-25cd-4a60-9a84-cfeecec75dc7
# Allowed
i_allowed = int_to_actions(Action, get_value(i_partition))

# ╔═╡ e72dd985-c655-4ac5-aad4-761a01832738
let
	header = if i_action ∈ i_allowed
	"""
	!!! success "Took shielded action"
	"""
	else
	"""
	!!! danger "Shield ignored :-("
	"""
	end

	"""
	$header
		`action: $(uppaal_example_trace.actions[i])`
		`action: $(uppaal_example_trace′["action"][i])`
	
		`x1: $(uppaal_example_trace.states[i][1])`
		`x2: $(uppaal_example_trace.states[i][2])`
		
		`f(x): $(f(uppaal_example_trace.states[i]))`
	
		`step: $(uppaal_example_trace′["step"][i])`
	
		`unsafe_state_entered: $(uppaal_example_trace′["unsafe_state_entered"][i])`
	
		`allowed.stay_course: $(uppaal_example_trace′["allowed.stay_course"][i])`
	
		`allowed.move_in: $(uppaal_example_trace′["allowed.move_in"][i])`
	
		`allowed.move_out: $(uppaal_example_trace′["allowed.move_out"][i])`

		`i_allowed: $([string(a) for a in i_allowed])`
	""" |> Markdown.parse
end

# ╔═╡ 4bd9a1a7-d1b1-415a-a238-315d13521299
shieldcolors[get_value(i_partition) + 1]

# ╔═╡ 62aa856d-a577-4f8a-b57e-85206fe3347c
# Reachable partitions
let
	round_4(a) = round(a, digits=4)
	prettyprint(b::Bounds) = "Bounds($(round_4.(b.lower)), $(round_4.(b.upper))"
	reachable = a_reachability_function(i_partition, i_action)
	reachable = [Partition(a_shield, indices) for indices in reachable]
	reachable = [(get_value(partition), Bounds(partition)) for partition in reachable]
	
	reachable = [(shieldcolors[value + 1], prettyprint(bounds)) 
		for (value, bounds) in reachable]
end

# ╔═╡ aaac12e2-3ece-407d-839a-b4676d70efb3
md"""
## Check reachability uppaal etc
"""

# ╔═╡ 3a408283-6921-4ebb-bfca-b547922486a9
let
	# Reachability in altered state space.
	# Paste in states from uppaal and write the correct action
	action = stay_course
	
	# before
x1 = -1.3120041233945752
x2 = 0.18522502698172874
	x = (x1, x2)
	# after
x1 = -1.3196053587033947
x2 = 0.11941916379615652
	x′ = (x1, x2)

	partition = box(a_shield, f(x))
	partition′ = box(a_shield, f(x′))

	reachable = a_reachability_function(partition, action)
	@info partition′.indices ∈ reachable
	@info get_value(partition′)
	reachable = [Partition(a_shield, i) for i in reachable]

	plot(Bounds(partition), color=colors.NEPHRITIS, label=nothing)
	for partition″ in reachable
		plot!(Bounds(partition″), color=colors.BELIZE_HOLE, label=nothing)
	end
	plot!(Bounds(partition′), color=colors.ALIZARIN, opacity=0.5, label=nothing)

	supporting_points = SupportingPoints(6, partition)
	successors = [f(successor(f⁻¹(x), action)) for x in supporting_points]
	scatter!([Tuple(x) for x in supporting_points], color=colors.EMERALD, label=nothing)
	scatter!([Tuple(x) for x in successors], color=colors.PETER_RIVER, label=nothing)
	
	plot!([Tuple(f(x)), Tuple(f(x′))],
		linewidth=3,
		linestyle=:dot,
		color=colors.WET_ASPHALT,
		marker=:circle,
		ratio=1, 
		xlabel="θ",
		ylabel="r",
		label=nothing)
end

# ╔═╡ 140850b6-a3b4-4e07-a17a-5ad77e31df79
let
	
	# Overwrite successor function with higher time step
	successor′(x, δ) = (exp(A*δ))*[x...]
	
	# Reachability in original state space.
	# Paste in states from uppaal and write the correct action
	action = move_out
	
	# before
x1 = -0.5310075417723606
x2 = 0.8102798596338565
	x = (x1, x2)
	# after
x1 = -0.5651185227456723
x2 = 0.7748812344580948
	x′ = (x1, x2)

	partition = box(shield, x)
	partition′ = box(shield, x′)

	reachable = reachability_function(partition, action)
	@info "partition′ in reachable?" partition′.indices ∈ reachable
	@info "actoins partition′" int_to_actions(Action, get_value(partition′))
	@info "actoins partition" int_to_actions(Action, get_value(partition))
	reachable = [Partition(shield, i) for i in reachable]

	plot(xlim=(x[1] - 0.1, x[1] + 0.1), ylim=(x[2] - 0.1, x[2] + 0.1),)
	plot!(Bounds(partition), color=colors.NEPHRITIS, label=nothing)
		
	plot!(rocks;
		seriestype=:shape,
		color=colors.CONCRETE,
		linewidth=0,
		label=nothing)
	
	for partition″ in reachable
		plot!(Bounds(partition″), color=colors.BELIZE_HOLE, label=nothing)
	end
	plot!(Bounds(partition′), color=colors.ALIZARIN, opacity=0.5, label=nothing)

	supporting_points = SupportingPoints(6, partition)
	successors = [successor(x, action) for x in supporting_points]
	scatter!([Tuple(x) for x in supporting_points], color=colors.EMERALD, label=nothing)
	scatter!([Tuple(x) for x in successors], color=colors.PETER_RIVER, label=nothing)
	
	plot!([Tuple(x), Tuple(x′)],
		linewidth=3,
		linestyle=:dot,
		color=colors.WET_ASPHALT,
		marker=:circle,
		ratio=1, 
		xlabel="x1",
		ylabel="x2",
		label=nothing)	

	plot!([Tuple(successor′(x′, δ′)) for δ′ in 0:0.01:δ],
		label=nothing,
		linestyle=:dash,
		color=colors.WET_ASPHALT)
end

# ╔═╡ Cell order:
# ╟─ac05e283-3e2f-4f8d-9ccc-c6c96bf2e7ee
# ╟─d7a85669-fd51-4730-ac8b-9542b82e279a
# ╠═19bd1463-ceb6-46ff-871f-3f0117ebeac9
# ╠═8cf97456-6c01-4cd4-b3ac-97e1115620a8
# ╟─a115e214-31a5-4a62-a798-07ebbc67caa1
# ╟─8b405589-2ba9-4046-ad9d-e2e3ccd21c84
# ╠═d355823b-cc1c-40b9-96e0-a09cd79be7ff
# ╠═52be2607-df01-458b-898e-901406747e3e
# ╠═ec77c85b-f27e-406b-98c3-5efbb8fa860a
# ╠═acd8f481-ba99-42cf-b95a-8c175b996e1d
# ╠═cb7fdd55-db5c-40ae-938c-0f3dd2b26b59
# ╠═9c3302dc-c01c-46f4-b9fd-42b4f3f7dd61
# ╠═61fa2c7c-4c61-4b6c-9db1-f35c011df04a
# ╠═4f87cc27-d210-4273-a4fc-89622bb8ef6e
# ╠═b047d22d-f223-446e-b58e-daf6a77b898d
# ╠═88ffd438-a6b3-4b11-82d8-af9ce61ed222
# ╠═bb1ff345-58ba-498a-8035-88ef1db7d917
# ╠═4937d230-e39f-4460-8071-a4cadc7e8b6f
# ╠═bf560bb3-2fcf-4cf3-b6a4-5bacb7cbc832
# ╟─f7ef34f2-533c-43bd-bdb5-f407671facb1
# ╠═fd8f218e-42d5-4432-b8ca-2ede5f746417
# ╠═4e5dc265-42a2-4874-a5cf-ae05d3ce4fdf
# ╟─61f009cc-a578-49e4-963a-d5b995c2ae39
# ╠═ffcb7ecc-4ea8-4e66-b39a-809810e8f105
# ╠═87962272-e8a3-4687-91dd-b7f6c23af97c
# ╠═9ff397e2-9aa7-433f-9a59-a4d3cb38a9bb
# ╠═905533f4-bf2b-4540-b4f0-d18f7fe039d1
# ╠═0a46d16c-e86d-443b-9001-6663a65ccccc
# ╠═b8b4af3d-9f66-4e3a-8845-0a45703bc82a
# ╠═cbe5bb1f-4a78-4185-9067-ddc7aa332be5
# ╠═f0b44048-4c9c-475f-87a4-9ff331885c32
# ╠═b67bbfee-633d-4d8b-846f-6a63ec778f93
# ╟─7557a054-9b6d-4858-ba82-1b613c362b4b
# ╠═7dddbfa0-35f6-4cec-90b6-37ae41881a77
# ╠═e8a5370b-e18d-4957-b695-d6c50fec1182
# ╠═01c9a1be-4258-4aef-9855-db55874f6ae0
# ╠═0f42e2d4-6446-401e-a062-b5aa893e9ac5
# ╠═1c429eaf-f72f-4197-8242-12f41db29f81
# ╠═3c15b3dc-d3f0-4a0c-89cd-da2bbd6e4865
# ╟─332ac25d-b07f-40ee-ab7d-e6cdfb2d075f
# ╠═5084de0e-e5e8-435b-9264-d858362957fb
# ╠═1ce94535-be9e-4f07-b0f3-062509540516
# ╠═7d5ba76c-e4cb-4c12-928e-2b9c30435288
# ╠═44db0342-b1f1-4e0c-a422-b95bf4c78dbe
# ╠═0d0aac2a-15ae-422a-b33f-ef2810920d15
# ╠═1c586fc6-1ec6-4f23-910b-4a9efb47e9fa
# ╟─d652b956-7473-4042-91b1-b508067e9a62
# ╠═d5d84407-294d-4247-b00f-291530edf099
# ╠═af490584-bd8c-4e03-a64b-32bab91afc33
# ╠═398ce525-47fc-494a-ad4c-312b6c49c566
# ╠═74bebc78-c912-4ac8-9c94-1f85838e836b
# ╟─9ea2e239-00f8-4517-aeac-566415a5fa8a
# ╟─4bf2b834-15a3-4d1a-8314-8b851a9ca33b
# ╠═8bb0f2ad-5802-492d-b209-158d35b66a18
# ╠═7dadc8b9-c52b-4578-8da4-36862f4c59f0
# ╟─abc8c87f-df92-4320-91e9-7f1d5ad4462f
# ╠═b538485b-0c21-49c7-b9c0-bd7f13e2da3b
# ╠═f161fcf5-83f2-4a44-87a5-58d0a052d01b
# ╠═14013a0c-9a27-44c5-b9ac-5d844fd3fe30
# ╠═2f658797-0581-40c4-9054-4ed968ad143b
# ╠═531281a1-6f2d-49fc-b608-cb1f8fc966ae
# ╟─7d1b7c12-9bf3-4b7a-891a-7139acefffd3
# ╠═6e5eee38-2267-45a2-9a0c-a2663ca207a8
# ╟─3ba16e83-5a6f-4bff-a2ad-dd20785f33dc
# ╠═744dcdcf-1db6-4544-b4f0-d3d49149b1d7
# ╟─d5101b8b-e65d-4a61-a616-803c9415bdcc
# ╠═d9506fdc-4bd4-40f6-9686-0403628eecf7
# ╟─19e227b0-df36-4835-9a4d-514490c67062
# ╠═60a01054-37be-4bcd-98fb-023350073ecf
# ╠═4c6ca97b-4871-4f4c-85e8-ee02719a2433
# ╟─c780a26f-d6c8-4573-b919-2500b9e4d8f4
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
# ╠═c5570c5b-5d4c-44bf-9dde-2d5e9a8e8ed4
# ╠═6c3550be-a891-4019-a97d-d4f28e3b1912
# ╟─98822937-709b-4654-b3dc-e23342d28f0f
# ╠═e1a91b83-07d3-48f0-9c25-e27f3ce0258f
# ╟─ee3c081b-9bc1-4dfc-b056-620eb305b4cd
# ╟─a5b2aa32-8267-45da-8664-98ea3afe4671
# ╠═c8b232d3-fe05-453e-9bed-837c83a81a6e
# ╠═aa8ebc72-8300-405b-b6da-704e39f19506
# ╠═fe72c613-014f-4c07-bee5-45efeb2bb770
# ╠═13646ad3-c39b-46d0-bddf-c16206dbb7eb
# ╠═87356f08-c963-4928-8a0f-48ca9b8d82b2
# ╠═d2427d85-7e20-4d2e-865c-2c582f65fe87
# ╠═9ce546c6-361d-4186-b585-2534d38614b6
# ╠═875fe159-8fb4-437b-a327-3e2129957484
# ╟─319cebe6-33ce-4a6f-be54-e9fb6dd059b0
# ╠═c8867810-099d-4dfc-909a-a5e0a0773bb8
# ╠═a7b0f796-6532-4a87-979c-f9c8759813f6
# ╟─ddf10c36-f525-4384-b5ac-731abc9d2b1e
# ╠═cda1ca0f-b02b-42cf-85fb-20175578396b
# ╠═3061a25a-bae8-4763-b84b-59e3536fd63c
# ╠═150a3ed4-3a97-4a75-b836-3c990d127f8f
# ╟─c9412832-1653-449b-bff3-99147ec7f3a6
# ╠═4ec52885-8315-403d-8eeb-d49a4eac8f16
# ╟─f47335e6-999d-449d-bf95-1f184b898042
# ╠═e514a65b-eb40-4193-92cd-406273e43d9b
# ╟─3b0d8736-713c-4fde-85db-7cfb239c45fe
# ╟─551adaa4-3c3c-47a1-a263-db3fefcaf4f0
# ╟─04d7c2d6-fb7c-4329-b349-0d8563c8cde8
# ╠═ab11ba27-c11f-41cb-8266-9b8c9e39008a
# ╠═c58c5add-cca1-4102-a71a-af82e3115993
# ╟─fe28308c-829d-4b9d-a821-9f2aba4204aa
# ╠═35079353-ac6d-4c09-b25f-cdcaf526c914
# ╠═abefb2a8-19fd-4c92-a81a-b0513572b756
# ╠═206bb5c4-f1fd-4093-b7d4-726a39a2bb16
# ╠═c17af70b-2bc5-4949-8133-499c8ea3c0e2
# ╟─e5f19917-868f-41f0-8873-6eef2d40bfaf
# ╠═46e15850-9632-4232-9d8f-05f3a9cb4832
# ╠═852cb308-c800-42ee-9412-412fe9b9a05e
# ╟─641033dd-aa48-43fd-b81b-c3ad84b58327
# ╟─80934c04-60dd-4fbc-9563-0954aa72eec2
# ╟─b65755fc-1623-4f67-9aec-54dcd3872ae0
# ╠═52ac166d-61f6-4a4f-b65d-7b8e1db22261
# ╠═b015bc8a-bc0e-436e-8edd-52c6de7cb528
# ╠═bfdd4c8f-b980-439b-b774-598a77e8ac08
# ╠═24730e1a-b691-4fb0-8558-12724fe07697
# ╠═c2a9ad89-308c-4291-9efe-cb3a3cd36187
# ╠═01bac781-c4b0-4378-8e42-8dcc6eb5ac6e
# ╠═53cee99d-cb92-4d43-98ac-8e86b04e298a
# ╠═2e94277a-f16d-47ab-b215-9bb17e0537df
# ╟─bb4d99ba-9752-430c-9898-edf372c2ee41
# ╟─d629a401-e2b3-4020-9bdd-5f2101aba751
# ╠═877c0c7b-563c-4b6b-8195-14577d04ef14
# ╠═21cdf64b-ba00-447d-9b04-aac3a7a0269e
# ╠═4bea0ab7-3254-46fb-b48f-a04b5cf90819
# ╠═620bb25c-2bf9-4fc1-9543-1254d8eeb3db
# ╠═3e0a5ab6-7158-46fb-9684-2027fee6c340
# ╠═410f184a-2055-419a-a986-65b817bdb124
# ╠═5e29d540-3ad9-4d89-b829-ae6fe61bdf36
# ╠═12844f71-adeb-41d2-afd3-b6bc6c53fa92
# ╟─e72dd985-c655-4ac5-aad4-761a01832738
# ╠═feabe121-3aae-492d-8f68-245256a7a0d8
# ╠═75110d13-3374-4497-89c5-f82c9d726782
# ╠═24bceb53-25cd-4a60-9a84-cfeecec75dc7
# ╠═4bd9a1a7-d1b1-415a-a238-315d13521299
# ╠═62aa856d-a577-4f8a-b57e-85206fe3347c
# ╟─aaac12e2-3ece-407d-839a-b4676d70efb3
# ╠═3a408283-6921-4ebb-bfca-b547922486a9
# ╠═140850b6-a3b4-4e07-a17a-5ad77e31df79
