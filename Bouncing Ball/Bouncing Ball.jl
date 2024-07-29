### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ 9c8abfbc-a5f0-11ec-3a9b-9bfd0b447638
begin
	using Pkg
	Pkg.activate("..")
	using GridShielding
	using Plots
	using PlutoUI
	using Measures
	using Unzip
	using PyCall
	using JSON
	using ProgressLogging
	include("../Shared Code/FlatUI.jl")
end

# ╔═╡ c663a860-4562-4de0-9b08-edc041cde9e6
md"""
# Preamble
"""

# ╔═╡ bffbac67-8a3b-4155-a665-0c39f93d3dd7
TableOfContents()

# ╔═╡ 43b8a393-b906-4c46-b8fb-d277954510b3
md"""
## Make paper-friendly figures? 
"""

# ╔═╡ a0eeeee8-04b3-4cfd-a75f-03c53a6a91e0
@bind make_paper_friendly_figures CheckBox(default=false)

# ╔═╡ 0ddc1964-874a-4b47-ad13-0a070d42a51b
begin
	default_font = default(:fontfamily)
	default_size = default(:size)
	default_margin = default(:margin)
	
	paper_font = "serif-roman" 	# https://gr-framework.org/fonts.html
	paper_size = (300, 220)
	paper_margin = 0mm
end;

# ╔═╡ 28d80bda-7ced-4ac1-b538-9f1c5187a4a4
theme_type = if make_paper_friendly_figures
	Plots.default(fontfamily=paper_font)
	Plots.default(size=paper_size)
	Plots.default(margin=paper_margin)
	"Paper-firendly it is!"
else
	Plots.default(fontfamily=default_font)
	Plots.default(size=default_size)
	Plots.default(margin=default_margin)
	"Using Julia Plots defaults :-)"
end

# ╔═╡ e59ef411-1779-4ef6-8d82-37892f1387e8
md"""
## Some BB-mechanics
"""

# ╔═╡ 6fee7dcf-a0ee-431a-a5b7-d31c54ffa1a6
BB = GridShielding.BB

# ╔═╡ 67d30df1-b60a-4835-a331-94957908ae4a
const m = BB.bbmechanics

# ╔═╡ 19281cd7-e79d-4535-9c53-9e9de9882eb0
const g = m.g

# ╔═╡ dc918da4-5ab4-4795-a220-67ffbccb97d1
position_from_e_pot(e) = e/g

# ╔═╡ 0a785d98-a03f-4cb1-8d5a-1f61fda1ab4b
md"""

# A New State Space

$(E_{mek},~~ v)$

"""

# ╔═╡ d30c6363-75e6-42f8-b3a0-4a032df063ef
begin
	π_xlabel = "\$E_{m}\$"
	π_ylabel = "\$v\$"
	
	π_xlabel_simple = "E_m"
	π_ylabel_simple = "v"
end;

# ╔═╡ 61fa5c6f-2f91-4429-a479-10266f6332c8
md"""
## $E_{mek}$, Mechanical Energy

$E_{mek} = E_{pot} + E_{kin}$

$E_{pot} =  g p$

$E_{kin} = {1 \over 2}  v^2$
"""

# ╔═╡ f93a65f3-9bdf-493b-994c-a26f34818a96
e_pot(g, p) = abs(g)*p

# ╔═╡ a3af719b-3b92-4c39-a95e-478d5b3179a2
e_kin(g, v) = 0.5*v^2

# ╔═╡ b775a061-3279-4121-806c-e99d211c36b0
e_mek(g, v, p) = e_kin(g, v) + e_pot(g, p)

# ╔═╡ 5a251063-64e8-4ced-8e28-34cb4812f931
position_from_e_pot(e_mek(g, 0, 4))

# ╔═╡ 3544f929-e518-485f-bdec-eaf1506f3226
md"""
`v_0 =` $(@bind v_0 NumberField(-13:1:13, default = 4))
`p_0 =` $(@bind p_0 NumberField(0:1:13, default = 4))

"""

# ╔═╡ b60a9495-7d59-4faa-a399-ac83a83d934d
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

# ╔═╡ c1a7ffdd-767d-418d-96af-f13b357e980e
@bind gradient Select([[:black, :deepskyblue, :white], :heat, :matter, :curl, :dense, :phase, :algae])

# ╔═╡ cad96c13-e9fa-45ae-b046-f976ae2ee901
theme_type; p2 = draw_function((vp) -> e_mek(g, vp...), -15, 15, 0, 10, 0.01,
	color=cgrad(gradient, 10, categorical=false),
	xlabel="Velocity (m/s)",
	ylabel="Position (m)",
	colorbar_title="Mechanical Energy (J)")

# ╔═╡ 490abcb1-80ea-4bbe-9b4f-b8133d22d9dd
theme_type; p1 = draw_function((vp) -> e_mek(g, vp...), -15, 15, 0, 10, 0.01 ,
	color=cgrad(gradient, 20, categorical=true),
	xlabel="Velocity (m/s)",
	ylabel="Position (m)",
	colorbar_title="Mechanical Energy (J)")

# ╔═╡ 2556f5de-5e22-4f88-b4bf-f3f4c87d06be
p1, p2; @bind ExportButton CounterButton("Export")

# ╔═╡ 3cfe65d2-7f6b-47c9-9c5f-ebc09229a2e6
if ExportButton > 0 let
	Plots.svg(p1, "../Graphics/Bouncing Ball/BB Mechanical Energy Grouped.svg")
	Plots.svg(p2, "../Graphics/Bouncing Ball/BB Mechanical Energy Smooth.svg")
end end

# ╔═╡ 0ad2d3c9-3a81-4582-b7dc-52225e0c99e9
velocity_from_e_kin(e) = sqrt(2*e)

# ╔═╡ 72d3376b-c91a-43c6-b237-53e83970fd4f
velocity_from_e_kin(e_kin(g, 4))

# ╔═╡ a92f10f7-46c7-4e66-b406-094e66da8cfe
velocity_from_e_kin(e_kin(g, -4))

# ╔═╡ 9cbb038a-62d3-43bf-9281-e672db8fb19e
md"""
## Transformation Function

This function is called **f** and it's inverse **f⁻¹** everywhere else. I just chose $\pi$ here  because I thought I was doing a projection initially.
"""

# ╔═╡ 78cb48d3-bedf-48e9-9479-8c71bcc10f6f
function π(v, p) 
	
	return e_mek(g, v, p), v
	
end

# ╔═╡ 236497fd-670b-45b1-8ca3-41b204a4d287
function π⁻¹(e_mek, v)
	p = (e_mek - 1/2*v^2)/-g
	#if p < 0 return nothing end
	return v, p
end

# ╔═╡ b527f190-ff38-48d3-97ae-aeeed8fdd273
md"""
### How traces look in the new projection
"""

# ╔═╡ aad4b9e6-2fbb-46a9-9311-f9e534a17002
md"""
# Shielding in the Transformed State Space
"""

# ╔═╡ f363e7ad-ad45-4fca-83c3-7b04ffdf48eb
bbshieldlabels = 	[
	"{$(join(int_to_actions(BB.Action, i), ", "))}"
	for i in 0:3]

# ╔═╡ 8f0f7850-c149-4735-a2d5-f58182251d34
bbshieldcolors = [colors.WET_ASPHALT, colors.AMETHYST, colors.SUNFLOWER, colors.CLOUDS];

# ╔═╡ 5c3ed2b7-81c0-42e5-b157-d65e25537791
vp_bounds = Bounds((-13, 0), (13, 10))

# ╔═╡ 98663ba0-e134-45eb-b700-8fb78dcc6971
md"""
## The Transformed Grid
"""

# ╔═╡ 2408a96c-8634-4fe9-91aa-af32ac2c7dec
const π_bounds = let
	#e_mek_upper = e_mek(g, vp_bounds.upper...)
	e_mek_upper = 100.
	v_lower, v_upper = vp_bounds.lower[1], vp_bounds.upper[1]
	
	Bounds((0., v_lower), ceil.((e_mek_upper, v_upper)))
end

# ╔═╡ 7b3c0c32-6845-46a6-9507-a4351c14dd47
md"""
### Figures Showing how the Transformed Grid Covers the Original State Space
"""

# ╔═╡ 443301cb-ef1c-40b3-a552-f86e46e0cbe8
# ╠═╡ disabled = true
#=╠═╡
let
	valid_partitions = let
		result = Partition[]
		for partition in Grid([8, 8, 1], π_grid.bounds)
			points = BruteForceSampler(partition, global_supporting_points) |> collect
			if length(points) > 0
				push!(result, partition)
			end
		end
		result
	end;
	
	# This is the plot of reachable partitions. It's super slow.
	p1 = plot([], 
		title=π_zlabel,
		seriestype=:shape, 
		color=colors.PETER_RIVER, 
		xlabel=π_xlabel,
		ylabel=π_ylabel,
		label=nothing)
	
	plot!(π_bounds, color=:white, linecolor=:white, label=nothing)

	for partition in valid_partitions
		bounds = Bounds(partition)
		if bounds.lower[3] == 0 continue end
		plot!(bounds, color=colors.PETER_RIVER, label=nothing, lw=1)
	end
	plot!(legend=:bottomright)

	p2 = plot([], 
		seriestype=:shape, 
		title="not $π_zlabel",
		color=colors.PETER_RIVER, 
		xlabel=π_xlabel,
		ylabel=π_ylabel,
		label="Ball can actually have this state")

	
	plot!(π_bounds, color=:white, linecolor=:white, label=nothing)

	for partition in valid_partitions
		bounds = Bounds(partition)
		if bounds.lower[3] == 1 continue end
		plot!(bounds, color=colors.PETER_RIVER, label=nothing, lw=1)
	end
	plot!(legend=:bottomright)

	plot(p1, p2, size=(800, 400))
end
  ╠═╡ =#

# ╔═╡ 898fcb77-2f6d-42b9-93c7-dce396664174
md"""
## Reachability Function
"""

# ╔═╡ 206a65db-e953-4216-9689-31966739c88d
const samples_per_random_axis = [3]

# ╔═╡ c2d118ff-daaa-4649-8937-76f6f4de684b
samples_per_axis = [10, 10, 1]

# ╔═╡ f0612487-06c4-4330-a0f0-fc4dd367d083
prod([samples_per_axis..., samples_per_random_axis...])

# ╔═╡ 5bf69f54-8ec2-4561-b696-7199ce83c839
round_8(x) = round(x, RoundNearestTiesUp, digits=8)

# ╔═╡ 7a76b511-0462-4889-bfe6-939a89ca4702
md"""
!!! danger "wtf?"
	I have no idea why there are two different reachability function variants. 
	The changes seem arbitrary.
	I remember I was messing around to get a better result, but I can't for the life of me remember what.

	The prime-variant does not round the result afterwards, and they have different conditions for skipping a sample.
"""

# ╔═╡ 104f1f24-44c8-4ea8-9d6a-732984a96e91
@bind spa NumberField(1:1000, default=samples_per_axis[1]/2)

# ╔═╡ 6327ed76-cf69-4389-8ce2-e0e9c42eb11f
md"""
### Deprecated: Brute-force sampling

In cases where there is no easy $π⁻¹$ function available: Just iterate through samples that cover the whole state space and check for every sample if it is a member of the π-partition. Horribly inefficient, but it works for this example.
"""

# ╔═╡ 01190c0f-b8bb-403f-8eed-57d683ad302a
# Number of samples per unit
const global_sampling_resolution = 5

# ╔═╡ c98583c9-3105-46b3-80b4-06b84d6e1db6
global_supporting_points::Vector{Tuple{Float64, Float64}} = SupportingPoints([
		(u - l)*global_sampling_resolution
		for (l, u) in zip(vp_bounds.lower, vp_bounds.upper)
	], 
	vp_bounds) |> collect;

# ╔═╡ f351c1ed-89d0-495c-8720-7f1ffa9ddd93
begin
	# Brute-force approach to sample generation
	# by checking membership of every single sample in the state space
	# for some resolution of samples
	struct BruteForceSampler
		partition::Partition
		points::Vector{Tuple{Float64, Float64}}
	end

	
	Base.IteratorSize(::BruteForceSampler) = Base.SizeUnknown()
	
	Base.iterate(sampler::BruteForceSampler) = begin
		Base.iterate(sampler, 1)
	end

	Base.iterate(sampler::BruteForceSampler, i::Int64) = begin
		l = length(sampler.points)
		while π(sampler.points[i]...) ∉ sampler.partition && i < l
			i += 1
		end
		if i == l
			nothing
		else
			(sampler.points[i], i + 1)
		end
	end

	BruteForceSampler
end

# ╔═╡ f9fbc97a-3b7d-45f4-a784-45d2cf515fa1
scatter(global_supporting_points,
	marker=(2, :red, :circle),
	markerstrokewidth=0,
	xlabel="v",
	size=(300, 300),
	ylabel="p",
	label="global supporting points")

# ╔═╡ 94ced2a5-7ad8-49e3-b1d1-e1d5b8ee9868
md"""
## Synthesize a Shield
"""

# ╔═╡ 7dad96ac-3c70-4b75-86a1-3ab374d631fa
@bind max_steps NumberField(0:1000, default=1000)

# ╔═╡ 3f4d6c5b-b4a9-42c2-909e-569061590af7
@bind show_point CheckBox()

# ╔═╡ 60401048-7e4a-45c8-a0aa-4fb9338714ab
#=╠═╡
v = shielded_trace[1][i]
  ╠═╡ =#

# ╔═╡ a31a8a05-c145-43a9-b844-ccfaf9f49645
#=╠═╡
p = shielded_trace[2][i]
  ╠═╡ =#

# ╔═╡ 8790b998-d96e-4437-b9bb-d77571d4bd1b
# ╠═╡ disabled = true
#=╠═╡
@bind i NumberField(1:length(shielded_trace[1]), default=21)
  ╠═╡ =#

# ╔═╡ fd928206-accf-44fc-8762-599fe34c26b6
@bind action Select(BB.Action |> instances |> collect, default=BB.nohit)

# ╔═╡ 7802329e-9ef1-40a5-8d5f-79010fa6ac1f
BB.simulate_point(m, (v_0, p_0), action)

# ╔═╡ 22d05a23-bcad-4281-8303-5082a3d8e785
@bind v NumberField(-15:0.2:15, default=0)

# ╔═╡ 2a4c1d40-bd6d-4e83-94d8-c6a3cfa8aee0
@bind p NumberField(0:0.1:8, default=7)

# ╔═╡ 2565bc31-e784-4c58-b741-1ccf41813dca
md"""
## Visualise shield in original state-space
"""

# ╔═╡ 8c03490e-0c85-4287-b8a3-03c82ad05f07
# Returns samples per axis but only on the edges.
function get_edge_points(samples_per_axis, bounds)
	samples = SupportingPoints(samples_per_axis, bounds)
	# Have to do it like this to get a clockwise order
	lower1 = [s for s in samples if s[1] ≈ bounds.lower[1]] |> sort
	upper1 = [s for s in samples if s[1] ≈ bounds.upper[1]] |> sort |> reverse
	lower2 = [s for s in samples if s[2] ≈ bounds.lower[2]] |> sort |> reverse
	upper2 = [s for s in samples if s[2] ≈ bounds.upper[2]] |> sort
	vcat(lower1, upper2, upper1, lower2)
end

# ╔═╡ 06ab9ae5-fe2a-4734-8929-ca5babe0d9f5
function transformed_grid(grid, f; samples=4, outer_bounds=nothing)
	result = Shape[]
	for partition in grid
		bounds = Bounds(partition)
		points = get_edge_points(samples, bounds)
		points = [f(p) for p in points]
		if !isnothing(outer_bounds) && !any(p ∈ outer_bounds for p in points)
			continue
		end
		shape = Shape(points)
		push!(result, shape)
	end
	result
end

# ╔═╡ 080a4374-104e-4c30-b946-313475fb0c11
any_action, no_action = actions_to_int([BB.hit BB.nohit]), actions_to_int([])

# ╔═╡ 26092473-69d3-4777-9890-48fa928ccc94
function initial_value_of_vp_partition(bounds::Bounds)
	vl, pl = bounds.lower
	vu, pu = bounds.upper

	if e_mek(g, vl, pl) < 0.5 || e_mek(g, vu, pl) < 0.5
		no_action
	else
		any_action
	end
end

# ╔═╡ fc8619b5-8dbc-47b3-b66e-24ceeeb45f7f
begin
	vp_grid = Grid(1, vp_bounds)
	initialize!(vp_grid, initial_value_of_vp_partition)
	vp_grid
end

# ╔═╡ f4364c08-d09b-4dcc-89ea-e3a58490d901
function reachability_function(partition, action)::Vector{Vector{Int64}}
	result = Vector{Int64}[]
	grid = partition.grid
	for π_point in SupportingPoints(samples_per_axis, partition)
		if π_point ∉ Bounds(partition) continue end
		point = π⁻¹(π_point...)
		if isnothing(point) continue end
		if point ∉ vp_grid.bounds continue end
		for r in SupportingPoints(samples_per_random_axis, Bounds((-1,), (1,)))
			point′ = BB.simulate_point(m, point, r, action)
			π_point′ = π(point′...)
			π_point′ = round_8.(π_point′)
			if π_point′ ∉ grid
				continue
			end
			partition′ = box(grid, π_point′)
			if partition′.indices ∈ result
				continue
			end
			push!(result, partition′.indices)
		end
	end
	result
end

# ╔═╡ 0335457d-5081-4f34-b086-7f597413c9f7
function reachability_function′(partition, action)::Vector{Vector{Int64}}
	result = Vector{Int64}[]
	grid = partition.grid
	for π_point in SupportingPoints([spa, spa, 3], partition)
		#if π_point ∉ Bounds(partition) continue end
		point = π⁻¹(π_point...)
		if isnothing(point) continue end
		for r in SupportingPoints(samples_per_random_axis, Bounds((-1,), (1,)))
			point′ = BB.simulate_point(m, point, r, action)
			if point′ ∉ vp_grid.bounds continue end
			π_point′ = π(point′...)
			if π_point′ ∉ grid
				continue
			end
			partition′ = box(grid, π_point′)
			if partition′.indices ∈ result
				continue
			end
			push!(result, partition′.indices)
		end
	end
	result
end

# ╔═╡ 814e17fe-4824-410d-a46f-da73729d6e8c
function initial_value_of_π_partition(bounds::Bounds)::Int64
	e_mek_lower, _ = bounds.lower
	if e_mek_lower < 0.5
		no_action
	else
		any_action
	end
end

# ╔═╡ 3e00e758-2e2e-42da-9152-fff188f75875
begin
	π_grid = Grid([4, 1], π_bounds)
	initialize!(π_grid, initial_value_of_π_partition)
	π_grid
end

# ╔═╡ d9e38cab-da45-4367-ad44-0c02ced097b7
π_grid.bounds

# ╔═╡ cc239362-e3a9-4e7b-bef7-737233e2d338
size(π_grid), length(π_grid)

# ╔═╡ 670639a2-dc12-45af-bb38-5d197ff41fd4
let	
	theme_type
	p1 = draw(π_grid,
		xlabel=π_xlabel,
		ylabel=π_ylabel,
		colors=[:white, :white], 
		#color_labels=bbshieldlabels,
		size=(400, 400),
		margin=4mm,
		show_grid=true,
		legend=:topright)
end

# ╔═╡ c2c5207f-ee2e-46fa-877a-18a9bc691e11
let
	theme_type
	π_points = [p for p in SupportingPoints([20, 20, 3], π_grid.bounds) 
		if p ∈ π_grid.bounds]
	
	# π_points = [p for p in π_points if p[3] == 0]
	
	vp_points = [π⁻¹(p...) for p in π_points]
	vp_points = [p for p in vp_points if !isnothing(p)]

	
	
	p1 = scatter(π_points |> unzip,
		xlabel=π_xlabel,
		ylabel=π_ylabel,
		label="samples")
	
	p2 = scatter(vp_points |> unzip,
		xlabel="v",
		ylabel="p",
		label="translated samples")
	
	plot!(vp_grid.bounds,
		label="bounds of (v,p) grid",
		alpha=0.8)

	plot(p1, p2, 
		size=(600, 300), 
		legend=:outertop)
end

# ╔═╡ 1f1c79cb-d4d4-4e1b-9a34-b958ed864a7d
let
	theme_type
	π_grid = Grid([6, 6], π_grid.bounds)
	
	p1  = plot(π_grid.bounds,
		color=:white,
		line=nothing,
		label=nothing,
		legend=:outertop,
		xlabel="\$E_{mek}\$",
		ylabel="v",)

	π_point = π(v, p)

	π_partition = box(π_grid, π_point)
	bounds = Bounds(π_partition)
	bounds = Bounds(bounds.lower[1:2], bounds.upper[1:2])

	plot!(bounds, 
		linewidth=0,
		color=:red, 
		label="partition")

	scatter!([π_point[1]], [π_point[2]], 
		marker=(4, :green, :circle),
		markerstrokewidth=1,
		label="π(v, p)")
	
	p2 = plot(vp_grid.bounds,
		color=:white,
		line=nothing,
		label=nothing,
		legend=:outertop,
		xlabel="v",
		ylabel="p")
	
	π_partition = box(π_grid, π(v, p))
	sp = global_supporting_points
		
	sp = BruteForceSampler(box(π_grid, π(v, p)), global_supporting_points) |> collect 
	if length(sp) > 0
		scatter!(sp |> unzip, 
			marker=(2, :red, :circle),
			markerstrokewidth=0,
			label="BruteForce sampled")
	else
		@info "???"
	end

	sp = SupportingPoints([8, 8, 3], π_partition)
	sp = [p for p in sp if p ∈ π_grid.bounds]
	sp = [π⁻¹(p...) for p in sp]
	sp = [p for p in sp if !isnothing(p)]
	sp = [p for p in sp if p ∈ vp_grid.bounds]
	
	scatter!(sp, 
		marker=(2, :blue, :circle),
		markerstrokewidth=0,
		label="π⁻¹ sampled",
		alpha=0.1)
	
	scatter!([v], [p], 
		label="(v,p)",
		marker=(4, :green, :circle),
		markerstrokewidth=1,)

	plot(p1, p2, size=(800, 400), margin=3mm)
end

# ╔═╡ 966304ab-8d5e-452b-9d47-c234a14626e6
begin
	boundsify(indices::Vector{Int64}) = Partition(π_grid, indices) |> Bounds
end

# ╔═╡ 77750ddb-f774-4c95-8963-a4fd45806bb6
π; π_grid; @bind do_it_button CounterButton("Do it")

# ╔═╡ e762cebe-cea0-48ea-952b-55d14fbba5bb
if do_it_button >= 0
	reachability_function_precomputed = 
		get_transitions(reachability_function, BB.Action, π_grid);
end

# ╔═╡ af696d4b-aa09-4339-b471-d9c91f065364
shield, max_steps_reached = make_shield(reachability_function_precomputed, BB.Action, π_grid; max_steps)

# ╔═╡ d0dd5ad2-97b6-4d7a-a97b-cb33b29230e6
function animate_trace(states, times; 
		left_background=nothing, 
		right_background=nothing)
	
	vs = [v for (v, p) in states]
	ps = [p for (v, p) in states]
	e_kins = [e_kin(g, v) for v in vs]
	e_pots = [e_pot(g, p) for (v, p) in zip(vs, ps)]
	e_meks = [e_pot(g, p) + e_kin(g, v) for (v, p) in zip(vs, ps)]
	πs = [π(v, p) for (v, p) in zip(vs, ps)]
	layout = 2
	
	x1, y1 = [s[1] for s in πs], [s[2] for s in πs]
	x1label=π_xlabel
	y1label=π_ylabel

	if isnothing(shield)
		x1lims=(minimum(x1) - 3, maximum(x1) + 3)
		y1lims=(minimum(y1) - 3, maximum(y1) + 3)
	else
		x1lims=(shield.bounds.lower[1], shield.bounds.upper[1])
		y1lims=(shield.bounds.lower[2], shield.bounds.upper[2])
	end
	
	x2, y2 = vs, ps
	x2label="t"
	y2label="p"
	x2lims=(minimum(x2) - 3, maximum(x2) + 3)
	y2lims=(minimum(y2) - 0, maximum(y2) + 3)

	animation = @animate for (i, _) in enumerate(times)
		
		p1 = isnothing(left_background) ? plot() : plot(left_background)
		
		plot!(x1[1:i], y1[1:i],
			xlims=x1lims,
			ylims=y1lims,
			xlabel=x1label,
			ylabel=y1label,
			color=colors.WET_ASPHALT,
			linewidth=2,
			markersize=2,
			markeralpha=1,
			markershape=:circle)
		
		scatter!([x1[i]], [y1[i]], marker=(3, :circle, :red))
		
		
		p2 = isnothing(right_background) ? plot() : plot(right_background)
			
		plot!(x2[1:i], y2[1:i],
			xlims=x2lims,
			ylims=y2lims,
			xlabel=x2label,
			ylabel=y2label,
			color=colors.WET_ASPHALT,
			linewidth=2,
			markersize=2,
			markeralpha=1,
			markershape=:circle)

		#hline!([4], label=nothing, color=colors.WET_ASPHALT)
		scatter!([x2[i]], [y2[i]], marker=(3, :circle, :red))
		
		plot(p1, p2, 
			layout=layout, 
			size=(800, 400), 
			legend=nothing)
	end
	
	gif(animation, joinpath(tempdir(), "trace.gif"), fps=10, show_msg=false)
end

# ╔═╡ a3e566e8-6b31-4d07-a2b9-b3b90f178d63
Bounds(box(shield, π(7, 0)))

# ╔═╡ 74c833b3-f3b9-4f3e-b579-b34ff2b17220
shield.array |> unique

# ╔═╡ ccbdaed6-0e4d-4bef-941c-9ad831b52e61
plot(transformed_grid(shield, x -> π⁻¹(x...)), 
	label=nothing, 
	fill=nothing, 
	linecolor=colors.SILVER,
	size=(200, 200))

# ╔═╡ e247dfa7-6000-4df1-8a28-328463e32c49
length(shield)

# ╔═╡ 702172e9-59d7-4a77-b663-a89f66132a1f
partition = box(shield, π(v, p))

# ╔═╡ 3fdb6a5a-81e6-43ab-b3f5-4118fe2275c7
reachable = boundsify.(reachability_function(partition, action))

# ╔═╡ 7f4b10fe-bed4-4f0a-bc4e-0a6f0d0ca8f1
reachable′ = boundsify.(reachability_function′(partition, action))

# ╔═╡ 5b65f23f-ecd1-4911-98e8-57a582cdb4d3
bounds = Bounds(partition)

# ╔═╡ 24350838-772a-4357-b4fd-5275d6a70393
length(π_grid)

# ╔═╡ 3961c068-f268-48c5-926c-99cd5c501018
let
	partition = box(π_grid, π(v, p))
	reachability_function(partition, action)
end

# ╔═╡ e494556c-1106-49ce-85b4-729136b9b0b3
md"""
## Apply the shield
"""

# ╔═╡ efef17e1-8cd7-4d5b-a805-3d4a7345cf9d
function apply_shield(shield::Grid, policy)
    return (s) -> begin
		a = policy(s)
		if π(s...) ∉ shield
			return a
		end
        allowed = int_to_actions(BB.Action, get_value(box(shield, π(s...))))
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

# ╔═╡ f5bd346f-ba38-42c5-8920-7ec127f8c547
random(s...) = if (rand(1:10) == 1) BB.hit else BB.nohit end

# ╔═╡ ff60b015-12cf-478b-9a60-93a9b93d0f5f
trace = BB.simulate_sequence(m, (0, 10), random, 10)

# ╔═╡ 87651747-c606-4f15-b335-649492faedd9
plot(); BB.animate_trace(trace.states, trace.times)

# ╔═╡ 937afb55-7775-482d-8674-260c8de29614
animate_trace(trace.states, trace.times)

# ╔═╡ 087cbfb4-9f42-4f9a-85cd-e92ff2004cc8
shielded_random = apply_shield(shield, random)

# ╔═╡ d4cbae79-3a44-4f1f-839e-3b652bf83a42
shielded_random((v, p))

# ╔═╡ 92f2f097-02e7-4c7b-a8f0-9d0be416444f
length(shield)

# ╔═╡ d7b1d3d3-4ced-47f0-918a-3c3aa8cae5ed
md"""
# Evaluation
"""

# ╔═╡ cf85b021-ab85-4926-86c4-16854fbbe545
let
	unsafe_initial_state = nothing
	v = 0
	for p in 7:0.1:10
		partition = box(shield, π(v, p))
		if get_value(partition) == no_action
			unsafe_initial_state = (v, p)
		end
	end
	unsafe_initial_state

	if isnothing(unsafe_initial_state)
		md"""!!! success "Initial states safe."
			🖒
		"""
	else
		Markdown.parse("""!!! warning "Initial state unsafe."
			The shield is too conservative and considers some of the initial states unsafe. 
			Example: **$(unsafe_initial_state)**
		""")
	end
end

# ╔═╡ cb82c3f3-a4c7-435d-af89-7dc8785d4c51
lazy(state) = BB.nohit

# ╔═╡ 568bbecc-0726-43d2-ba8e-cc2c468c44b2
shielded_lazy = apply_shield(shield, lazy)

# ╔═╡ 2a7637b1-4a00-43e5-b601-e1f1bf7144bb
function generate_shielded_trace()
	trace = BB.simulate_sequence(m, (0, rand(7:0.01:10)), shielded_lazy, 120)
	trace.states, trace.actions
end

# ╔═╡ e4dd2de5-f959-4bcf-b67e-a12fe3c9d188
function is_safe(state)
	v, p = state
	return !(p == 0 && abs(v) < 1)
end

# ╔═╡ 05b5e4d4-9bea-49b5-ae51-0daa2fb8478d
@bind checks NumberField(1:10^20, default=1000)

# ╔═╡ 10782c73-e1ef-4576-b139-485b559f6ba6
safety_report = evaluate_safety(generate_shielded_trace, is_safe, checks)

# ╔═╡ 777c1274-e8f8-441d-a79c-17248c85dd39
shielded_trace = safety_report.example_trace[1]

# ╔═╡ 4f238ba9-035c-4046-ba45-604bf514be67
begin
	trace_to = findfirst(s -> !is_safe(s), shielded_trace)

	trace_to = something(trace_to, length(shielded_trace))

	trace_from = max(1, trace_to - 100)
	(;trace_from, trace_to)
end

# ╔═╡ 62a83575-a88a-4ce3-8d03-0a88a4864762
length(collect(1:BB.bbmechanics.t_hit:120))

# ╔═╡ 048a3a18-64c5-4ef3-9b25-23306e100dd2
md"""
# Exporting the Shield
"""

# ╔═╡ 9b4e88ad-2de3-4b8d-8a69-bbb8660cc293
@bind target_dir TextField(95, default=mktempdir())

# ╔═╡ 1a3ebfb8-47b8-41a3-b63d-875b03187a4e
target_dir; @bind open_folder_button CounterButton("Open Folder")

# ╔═╡ 5789cc7e-4a58-4a83-9f3a-87c982028c59
if open_folder_button > 0
	run(`nautilus $target_dir`, wait=false)
end; "This cell opens `$target_dir` in nautilus"

# ╔═╡ d97db0c1-6bba-4c92-b8a6-0f60a9ea6357
md"""
### Export as serialized julia-tuple

Easy export and import between julia code.
"""

# ╔═╡ 2813bdf7-9530-40a0-bdb5-ab1213f54b31
let
	filename = "$π_xlabel_simple, $π_ylabel_simple.shield"
	
	robust_grid_serialization(joinpath(target_dir, filename), shield)
	
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ 15e42007-54fd-4af2-b7ea-b86ca62a9f19
md"""
### Export as a function in a shared-object library

Use this library to access the shield from C and C++ code.

The shield is compiled into a shared-object binary, which exports the function `int get_value(double v, double p)`. It takes the state-variables as input and returns the bit-encoded list of allowed actions. (See `int_to_actions`.)
"""

# ╔═╡ 75b00a19-f89b-4b9b-b85e-f84a7c5f2e9b
md"""
### Export to Numpy

Exports a zip-file containing a serialized numpy-array along with a JSON file with details on how to read it.
"""

# ╔═╡ b62de837-53d8-4e61-97a3-d629d9388165
md"""
# Make the Shield 'UPPAAL friendly'

Alright, so. UPPAAL uses a Runge-Kutta implementation to approximate the differential equations. This is fine  for a lot of things, but it slightly violates the mechanical energy invariant. 

Due to imprecision in the calculation of $v$ and $p$, the ball will sometimes appear to lose or gain energy. In some instances, this will bring it just under the threshold between having a safe amount of mechanical energy, and an unsafe amount.

For this reason, the purple _must hit_ area is made to "drip downwards" to ensure that if the ball passes slightly below a _must hit_ partition, the partition below will also force the ball to be hit.

This is a very specific hack to address an innacurate simulation of the state-space. I note that the shield is perfectly safe in this notebook, and that inspecting safety violations in UPPAAL yields traces that all go ever so slightly below the partition it was supposed to actually pass through.
"""

# ╔═╡ b55d3379-4a1d-40e4-a663-ec07c119df33
shield; @bind make_uppaal_friendly_button CounterButton("Make UPPAAL Firendly")

# ╔═╡ e1cc2426-a7fe-41c9-bfdc-6bbd64126123
function make_uppaal_friendly!(shield::Grid)
	must_hit = actions_to_int([BB.hit])
	updates = []
	for partition in shield
		if get_value(partition) != must_hit continue end

		# partition below
		indices′ = partition.indices |> collect
		indices′[1] -= 1
		partition′ = Partition(shield, indices′)
		
		set_value!(partition′, must_hit)
	end
	return shield
end

# ╔═╡ 1bd93c52-5484-4364-82e3-0407fb6d0779
begin
	🎈1 = "reactivity"
	if make_uppaal_friendly_button == 0
		md"This cell makes the shield 'uppaal friendly'"
	elseif make_uppaal_friendly_button == 1
		make_uppaal_friendly!(shield)
		md"Made the shield 'uppaal friendly'!"
	elseif make_uppaal_friendly_button > 1
		md"Already made the shield 'uppaal friendly'"
	end
end

# ╔═╡ 021e2fb4-1760-4421-916b-fb2ef306cb13
shield_plot_new_statespace = let
	🎈1, theme_type
	
	partition = box(shield, π(v, p))

	slice::Vector{Any} = partition.indices
	slice[1] = slice[2] = Colon()
	
	p1 = draw(shield,
		show_grid=true,
		xlabel=π_xlabel,
		ylabel=π_ylabel,
		yflip=true,
		colors=[colors.WET_ASPHALT, colors.AMETHYST,  colors.SUNFLOWER, colors.CLOUDS],
		legend=:bottomright)

	plot!([], seriestype=:shape, color=colors.WET_ASPHALT, label="{}")
	plot!([], seriestype=:shape, color=colors.AMETHYST, label="{hit}")
	plot!([], seriestype=:shape, color=colors.CLOUDS, label="{hit, nohit}")

	if show_point
		## Reachable partitions ##
		π_point = π(v, p)

		## Barbaric sampling endpoints ##
		points = BruteForceSampler(partition, global_supporting_points)
		foo = []
		for p in points
			for r in SupportingPoints(samples_per_random_axis, Bounds((0,), (1,)))
				p′ = π(BB.simulate_point(m, p, r, action)...)
				push!(foo, (p′[1], p′[2]))
			end
		end
		if length(foo) > 0
			p1 = plot(p1)
			
			scatter!(foo |> unzip,
				marker=(1, colors.PETER_RIVER, :circle),
				markerstrokewidth=0,
				label="Barbaric sample endpoints")
		end

		plot!(Bounds(partition),
			linewidth=0,
			color=colors.EMERALD,
			opacity=0.5)
		
		## The actual point ##		
		scatter!([π_point[1]], [π_point[2]],
			marker=(3, colors.EMERALD, :circle),
			markerstrokewidth=0,
			label="(v, p)")
	end

	plot(p1)
end

# ╔═╡ a566b33b-7005-43c3-afce-b8793447f615
shield_plot_old_statespace = let
	🎈1, theme_type

	xlim = -13, 13
	ylim = 0, 8
	
	draw_function(
		s -> (box(shield, clamp(shield.bounds, π(s...) |> collect)) |> get_value), 
		xlim..., ylim..., (make_paper_friendly_figures ? 0.01 : 0.05),
		color=cgrad([colors.WET_ASPHALT, colors.AMETHYST, colors.SUNFLOWER, colors.CLOUDS], 10, categorical=true),
		xlabel="v",
		ylabel="p",
		colorbar=nothing)

	plot!(transformed_grid(shield, x -> π⁻¹(x...); samples=3, outer_bounds=vp_bounds);
		xlim, 
		ylim,
		label=nothing, 
		fill=nothing, 
		linecolor=colors.SILVER)

	#=
	plot!([], seriestype=:shape, color=colors.WET_ASPHALT, label="{}")
	plot!([], seriestype=:shape, color=colors.AMETHYST, label="{hit}")
	plot!([], seriestype=:shape, color=colors.CLOUDS, label="{hit, nohit}")
	=#
	
	if show_point
		points = BruteForceSampler(partition, global_supporting_points)
		barbaric_sample_endpoints = []
		for p in points
			for r in SupportingPoints(samples_per_random_axis, Bounds((0,), (1,)))
				push!(barbaric_sample_endpoints, BB.simulate_point(m, p, r, action))
			end
		end
		scatter!(barbaric_sample_endpoints |> unzip,
			marker=(1, colors.PETER_RIVER, :circle),
			markerstrokewidth=0,
			label="Barbaric sample endpoints")
		
		scatter!([v], [p],
			marker=(5, colors.EMERALD, :circle),
			markerstrokewidth=0,
			label="(v, p)")
	end
	plot!()
end

# ╔═╡ b097e128-a1df-44f0-8fb7-347d9317abfc
animate_trace(shielded_trace[trace_from:trace_to], 
	[i*BB.bbmechanics.t_hit for (i, _) in enumerate(shielded_trace)][trace_from:trace_to], 
	left_background=shield_plot_new_statespace,
	right_background=shield_plot_old_statespace)

# ╔═╡ fbe0e11d-a06b-41fe-b349-ccbcc66ffd3f
begin
	🎈1
	shield_so = "shield.so"
	shield_so = joinpath(target_dir, shield_so)
	
	get_libshield(shield; destination=shield_so, force=true)
	
	"Exported `'$shield_so'`." |> Markdown.parse
end

# ╔═╡ 1cce35be-253e-4e75-8f0f-fdf1aed9799d
let
	🎈1
	filename = "shield.zip"
	
	numpy_zip_file(shield, joinpath(target_dir, filename); 
		variables=[π_xlabel_simple, π_ylabel_simple],
		binary_variables=[3], 
		actions=BB.Action, 
		env_id="Bouncing Ball")
	
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ 700c196c-dafe-4116-bac8-1024acee9642
md"""
# Analyze safety violation in UPPAAL

!!! warn "TODO"
	All this code for parsing UPPAAL traces is now available as a package: 

	[UppaalTraceParser.jl](https://github.com/AsgerHB/UppaalTraceParser.jl/)
"""

# ╔═╡ aefb7e86-276d-4536-9ea0-33487a5015a8
⨝ = joinpath

# ╔═╡ 3071c7f7-4d77-45fb-866f-a27f76270284
function multiline(str)
	HTML("""
	<pre style='max-height:30em; margin:8pt 0 8pt 0; overflow-y:scroll'>
	$str
	</pre>
	""")
end

# ╔═╡ a29a9de0-0fee-4887-bf9a-3780725e418f
← = push!

# ╔═╡ a683c891-05bf-4119-af89-002f6def9ef9
function parse_pair(str)
	left = match(r"\(([-0-9.e]+),", str)
	if isnothing(left) error("left side of pair not found in $str") end
	left = left[1]
	left = parse(Float64, left)
	right = match(r",([-0-9.e]+)\)", str)
	if isnothing(right) error("right side of pair not found in $str") end
	right = right[1]
	right = parse(Float64, right)
	(left, right)
end

# ╔═╡ 8a7697e1-ef70-40b0-9432-4fe2bc302bad
function parse_pairs(str)
	[parse_pair(s) for s in split(str, " ") if s != ""]
end

# ╔═╡ 4f368fec-9cca-4236-99fe-ffd35cc9335d
function get_pairs_for_array(output::String, keyword::String)::Vector{Tuple{Float64, Float64}}
	re_trace = Regex("\\Q$keyword\\E:\\n\\[0\\]:(?<values>.*)", "m")
	m = match(re_trace, output)
	return parse_pairs(m[:values])
end

# ╔═╡ 877d8fe0-dc38-4a1c-997d-bb58514b7b86
function value_at_time(trace::T, time::S) where T <: AbstractVector{Tuple{Float64, Float64}} where S <: Number
	
	time_before, value_before = last([(t, v) 
		for (t, v) in trace if t <= time])

	time_after, value_after = first([(t, v) 
		for (t, v) in trace if t >= time])

	Δt = time_after - time_before
	if Δt == 0 # Happens if there is an exact match
		return value_after
	end
	fraction = (time - time_before)/Δt
	return value_before + fraction*(value_after - value_before)
end

# ╔═╡ a833bed7-8f79-4163-b601-6370a337c944
function at_regular_intervals(trace::T, interval::S) where T <: 
		AbstractVector{Tuple{Float64, Float64}} where S <: Number

	t_max = trace[end][1]
	return [ value_at_time(trace, i) for i in 0:interval:prevfloat(t_max) ]
end

# ╔═╡ 43f7a7c3-df0d-4862-9f3c-9fba64c1eb2f
function parse_trace(output, Δt, keywords...)
	result = Dict{String, Vector{Float64}}()
	
	for keyword in keywords
		trace = get_pairs_for_array(output, keyword)
		result[keyword] = at_regular_intervals(trace, Δt)
	end
	result
end

# ╔═╡ 7ee88060-754f-488d-9341-633bf54c2318
@bind model_file TextField(80, 
	default=pwd() ⨝ "Uppaal model with grid-shield/Bouncing Ball Shielded.xml")

# ╔═╡ 38489318-eadb-4cff-bd54-c91d76fa8800
function copy_and_replace(input_file, outfile, replacements)
	file = input_file |> read |> String
	open(outfile, "w") do io
		for line in split(file, "\n")
			line′ = replace(line, replacements...)
			println(io, line′)
		end
	end
	outfile
end

# ╔═╡ e2959cd8-f3da-4a81-a539-671a59b4ddb5
begin
	model_file′ = copy_and_replace(model_file, target_dir ⨝ "BB.xml", 
		[r"/\*capture 1\*/.*/\* end 1\*/" => "import \"$shield_so\" "])
end

# ╔═╡ 80f08a4c-ba20-4408-853e-a694df474a02
@bind query TextField((95, 6), default="""
	Pr[<=100;1000] ([] number_deaths < 1)

	E[<=100;100] (max:LearnerPlayer.fired)

	simulate[<=100;1] { v, p, allowed.nohit }
""")

# ╔═╡ eadf88fa-1290-4e08-9af4-53e0f613dc17
function remove_single_line_breaks(str)
	line_break_placeholder = "¤NEWLINE¤"
	str= replace(str, r"\n\s*\n" => line_break_placeholder)
	str = replace(str, r"\n\s*" => " ")
	str = replace(str, line_break_placeholder => "\n")
end

# ╔═╡ ecfea12e-c933-49a0-b117-88658fd5723a
query |> remove_single_line_breaks |> multiline

# ╔═╡ 9733a667-964a-4007-9e51-f656c3838954
@bind verifyta TextField(80, default=homedir() ⨝ "opt/uppaal-5.0.0-linux64/bin/verifyta")

# ╔═╡ 572dac94-aded-44af-8855-8f2c9811f94d
# 0.01 is the default.
discretization = 0.005

# ╔═╡ 7efc9685-1031-4882-8202-2fcd83a728d8
@bind working_dir TextField(80, default=mktempdir())

# ╔═╡ b703f5a3-e869-4051-83a7-64d53f446000
query_file = let
	query_file = working_dir ⨝ "queries.q"
	write(query_file, remove_single_line_breaks(query))
	query_file
end

# ╔═╡ c2fad765-52d9-479e-a494-faf38736d58c
if isfile(query_file) && isfile(model_file′)
	output = Cmd(String[
		verifyta,
		"-s",
		split("--discretization $discretization --truncation-error $discretization --truncation-time-error $discretization", " ")...,
		model_file′,
		query_file
	]) |> read |> String;
else
	@info "not found" isfile(query_file) isfile(model_file′)
end

# ╔═╡ 1fee975e-294c-407e-ad1a-65f91c33fbd4
md"""
The first query should say `(1000/1000 runs) Pr([] …) in [0.996318,1] (95% CI)`, if the shield is safe.
"""

# ╔═╡ ab63a35f-b674-4de1-bb6c-04f45f034a1d
output[1:min(10000, length(output))] |> multiline

# ╔═╡ 6f1a76ea-5eb7-4908-9edf-c47b7595f913
tt = parse_trace(output, 0.1, "v", "p", "allowed.nohit")

# ╔═╡ 465f5546-e9a1-4ea6-98e9-75436240dd86
@bind ii NumberField(1:10000)

# ╔═╡ 0038f63a-d7fc-4e48-a487-6aa97df100c5
let
	plot(shield_plot_old_statespace)
	plot!(tt["v"][2:end], tt["p"][2:end])
	scatter!([tt["v"][ii]], [tt["p"][ii]])
end

# ╔═╡ 78f6195a-c9fc-4d42-9e96-58a79d4fa06e
tt["p"][ii]

# ╔═╡ bc9a89bd-6b93-4ff9-8326-b24f8d3c54d1
tt["allowed.nohit"][ii]

# ╔═╡ b4849f98-2174-4e55-9070-61b5bded79a4
vv, pp = tt["v"][ii], tt["p"][ii]

# ╔═╡ eb6cb789-6d91-4f96-bedb-b92ba5d1d69a
get_value(box(shield, π(vv, pp)))

# ╔═╡ 2896b900-df7c-4fee-a62c-8ea64f9ffddc
BB.simulate_point(m, (vv, pp), BB.nohit)

# ╔═╡ 1cd28cd7-c3d7-4599-9f7b-b1d68bc094a0
tt["v"][ii + 1], tt["p"][ii + 1]

# ╔═╡ Cell order:
# ╟─c663a860-4562-4de0-9b08-edc041cde9e6
# ╠═9c8abfbc-a5f0-11ec-3a9b-9bfd0b447638
# ╠═bffbac67-8a3b-4155-a665-0c39f93d3dd7
# ╟─43b8a393-b906-4c46-b8fb-d277954510b3
# ╠═a0eeeee8-04b3-4cfd-a75f-03c53a6a91e0
# ╠═0ddc1964-874a-4b47-ad13-0a070d42a51b
# ╠═28d80bda-7ced-4ac1-b538-9f1c5187a4a4
# ╟─e59ef411-1779-4ef6-8d82-37892f1387e8
# ╠═6fee7dcf-a0ee-431a-a5b7-d31c54ffa1a6
# ╠═67d30df1-b60a-4835-a331-94957908ae4a
# ╠═19281cd7-e79d-4535-9c53-9e9de9882eb0
# ╠═7802329e-9ef1-40a5-8d5f-79010fa6ac1f
# ╠═dc918da4-5ab4-4795-a220-67ffbccb97d1
# ╠═5a251063-64e8-4ced-8e28-34cb4812f931
# ╟─0a785d98-a03f-4cb1-8d5a-1f61fda1ab4b
# ╠═d30c6363-75e6-42f8-b3a0-4a032df063ef
# ╟─61fa5c6f-2f91-4429-a479-10266f6332c8
# ╠═f93a65f3-9bdf-493b-994c-a26f34818a96
# ╠═a3af719b-3b92-4c39-a95e-478d5b3179a2
# ╠═b775a061-3279-4121-806c-e99d211c36b0
# ╟─3544f929-e518-485f-bdec-eaf1506f3226
# ╠═b60a9495-7d59-4faa-a399-ac83a83d934d
# ╠═c1a7ffdd-767d-418d-96af-f13b357e980e
# ╟─cad96c13-e9fa-45ae-b046-f976ae2ee901
# ╟─490abcb1-80ea-4bbe-9b4f-b8133d22d9dd
# ╠═2556f5de-5e22-4f88-b4bf-f3f4c87d06be
# ╠═3cfe65d2-7f6b-47c9-9c5f-ebc09229a2e6
# ╠═0ad2d3c9-3a81-4582-b7dc-52225e0c99e9
# ╠═72d3376b-c91a-43c6-b237-53e83970fd4f
# ╠═a92f10f7-46c7-4e66-b406-094e66da8cfe
# ╟─9cbb038a-62d3-43bf-9281-e672db8fb19e
# ╠═78cb48d3-bedf-48e9-9479-8c71bcc10f6f
# ╠═236497fd-670b-45b1-8ca3-41b204a4d287
# ╟─b527f190-ff38-48d3-97ae-aeeed8fdd273
# ╠═ff60b015-12cf-478b-9a60-93a9b93d0f5f
# ╠═d0dd5ad2-97b6-4d7a-a97b-cb33b29230e6
# ╟─87651747-c606-4f15-b335-649492faedd9
# ╟─937afb55-7775-482d-8674-260c8de29614
# ╟─aad4b9e6-2fbb-46a9-9311-f9e534a17002
# ╠═f363e7ad-ad45-4fca-83c3-7b04ffdf48eb
# ╠═8f0f7850-c149-4735-a2d5-f58182251d34
# ╠═26092473-69d3-4777-9890-48fa928ccc94
# ╠═5c3ed2b7-81c0-42e5-b157-d65e25537791
# ╠═fc8619b5-8dbc-47b3-b66e-24ceeeb45f7f
# ╠═d9e38cab-da45-4367-ad44-0c02ced097b7
# ╟─98663ba0-e134-45eb-b700-8fb78dcc6971
# ╠═2408a96c-8634-4fe9-91aa-af32ac2c7dec
# ╠═814e17fe-4824-410d-a46f-da73729d6e8c
# ╠═3e00e758-2e2e-42da-9152-fff188f75875
# ╠═cc239362-e3a9-4e7b-bef7-737233e2d338
# ╟─670639a2-dc12-45af-bb38-5d197ff41fd4
# ╟─7b3c0c32-6845-46a6-9507-a4351c14dd47
# ╟─c2c5207f-ee2e-46fa-877a-18a9bc691e11
# ╟─1f1c79cb-d4d4-4e1b-9a34-b958ed864a7d
# ╟─443301cb-ef1c-40b3-a552-f86e46e0cbe8
# ╟─898fcb77-2f6d-42b9-93c7-dce396664174
# ╠═206a65db-e953-4216-9689-31966739c88d
# ╠═c2d118ff-daaa-4649-8937-76f6f4de684b
# ╠═f0612487-06c4-4330-a0f0-fc4dd367d083
# ╠═5bf69f54-8ec2-4561-b696-7199ce83c839
# ╟─7a76b511-0462-4889-bfe6-939a89ca4702
# ╠═f4364c08-d09b-4dcc-89ea-e3a58490d901
# ╠═0335457d-5081-4f34-b086-7f597413c9f7
# ╠═966304ab-8d5e-452b-9d47-c234a14626e6
# ╠═3fdb6a5a-81e6-43ab-b3f5-4118fe2275c7
# ╠═7f4b10fe-bed4-4f0a-bc4e-0a6f0d0ca8f1
# ╠═104f1f24-44c8-4ea8-9d6a-732984a96e91
# ╟─6327ed76-cf69-4389-8ce2-e0e9c42eb11f
# ╠═01190c0f-b8bb-403f-8eed-57d683ad302a
# ╠═c98583c9-3105-46b3-80b4-06b84d6e1db6
# ╟─f351c1ed-89d0-495c-8720-7f1ffa9ddd93
# ╟─f9fbc97a-3b7d-45f4-a784-45d2cf515fa1
# ╟─94ced2a5-7ad8-49e3-b1d1-e1d5b8ee9868
# ╠═7dad96ac-3c70-4b75-86a1-3ab374d631fa
# ╠═77750ddb-f774-4c95-8963-a4fd45806bb6
# ╠═e762cebe-cea0-48ea-952b-55d14fbba5bb
# ╠═af696d4b-aa09-4339-b471-d9c91f065364
# ╠═24350838-772a-4357-b4fd-5275d6a70393
# ╠═a3e566e8-6b31-4d07-a2b9-b3b90f178d63
# ╠═3f4d6c5b-b4a9-42c2-909e-569061590af7
# ╠═60401048-7e4a-45c8-a0aa-4fb9338714ab
# ╠═a31a8a05-c145-43a9-b844-ccfaf9f49645
# ╠═8790b998-d96e-4437-b9bb-d77571d4bd1b
# ╠═fd928206-accf-44fc-8762-599fe34c26b6
# ╠═22d05a23-bcad-4281-8303-5082a3d8e785
# ╠═2a4c1d40-bd6d-4e83-94d8-c6a3cfa8aee0
# ╠═74c833b3-f3b9-4f3e-b579-b34ff2b17220
# ╠═021e2fb4-1760-4421-916b-fb2ef306cb13
# ╟─2565bc31-e784-4c58-b741-1ccf41813dca
# ╟─8c03490e-0c85-4287-b8a3-03c82ad05f07
# ╠═06ab9ae5-fe2a-4734-8929-ca5babe0d9f5
# ╟─ccbdaed6-0e4d-4bef-941c-9ad831b52e61
# ╠═a566b33b-7005-43c3-afce-b8793447f615
# ╠═e247dfa7-6000-4df1-8a28-328463e32c49
# ╠═702172e9-59d7-4a77-b663-a89f66132a1f
# ╠═5b65f23f-ecd1-4911-98e8-57a582cdb4d3
# ╠═080a4374-104e-4c30-b946-313475fb0c11
# ╠═3961c068-f268-48c5-926c-99cd5c501018
# ╟─e494556c-1106-49ce-85b4-729136b9b0b3
# ╠═efef17e1-8cd7-4d5b-a805-3d4a7345cf9d
# ╠═f5bd346f-ba38-42c5-8920-7ec127f8c547
# ╠═087cbfb4-9f42-4f9a-85cd-e92ff2004cc8
# ╠═d4cbae79-3a44-4f1f-839e-3b652bf83a42
# ╠═92f2f097-02e7-4c7b-a8f0-9d0be416444f
# ╟─d7b1d3d3-4ced-47f0-918a-3c3aa8cae5ed
# ╟─cf85b021-ab85-4926-86c4-16854fbbe545
# ╠═cb82c3f3-a4c7-435d-af89-7dc8785d4c51
# ╠═568bbecc-0726-43d2-ba8e-cc2c468c44b2
# ╠═2a7637b1-4a00-43e5-b601-e1f1bf7144bb
# ╠═e4dd2de5-f959-4bcf-b67e-a12fe3c9d188
# ╠═05b5e4d4-9bea-49b5-ae51-0daa2fb8478d
# ╠═10782c73-e1ef-4576-b139-485b559f6ba6
# ╠═777c1274-e8f8-441d-a79c-17248c85dd39
# ╠═4f238ba9-035c-4046-ba45-604bf514be67
# ╠═b097e128-a1df-44f0-8fb7-347d9317abfc
# ╠═62a83575-a88a-4ce3-8d03-0a88a4864762
# ╟─048a3a18-64c5-4ef3-9b25-23306e100dd2
# ╠═9b4e88ad-2de3-4b8d-8a69-bbb8660cc293
# ╠═1a3ebfb8-47b8-41a3-b63d-875b03187a4e
# ╟─5789cc7e-4a58-4a83-9f3a-87c982028c59
# ╠═d97db0c1-6bba-4c92-b8a6-0f60a9ea6357
# ╠═2813bdf7-9530-40a0-bdb5-ab1213f54b31
# ╠═15e42007-54fd-4af2-b7ea-b86ca62a9f19
# ╠═fbe0e11d-a06b-41fe-b349-ccbcc66ffd3f
# ╠═75b00a19-f89b-4b9b-b85e-f84a7c5f2e9b
# ╠═1cce35be-253e-4e75-8f0f-fdf1aed9799d
# ╟─b62de837-53d8-4e61-97a3-d629d9388165
# ╟─b55d3379-4a1d-40e4-a663-ec07c119df33
# ╠═e1cc2426-a7fe-41c9-bfdc-6bbd64126123
# ╠═1bd93c52-5484-4364-82e3-0407fb6d0779
# ╟─700c196c-dafe-4116-bac8-1024acee9642
# ╠═aefb7e86-276d-4536-9ea0-33487a5015a8
# ╠═3071c7f7-4d77-45fb-866f-a27f76270284
# ╠═a29a9de0-0fee-4887-bf9a-3780725e418f
# ╠═a683c891-05bf-4119-af89-002f6def9ef9
# ╠═8a7697e1-ef70-40b0-9432-4fe2bc302bad
# ╠═4f368fec-9cca-4236-99fe-ffd35cc9335d
# ╠═877d8fe0-dc38-4a1c-997d-bb58514b7b86
# ╠═a833bed7-8f79-4163-b601-6370a337c944
# ╠═43f7a7c3-df0d-4862-9f3c-9fba64c1eb2f
# ╠═7ee88060-754f-488d-9341-633bf54c2318
# ╠═38489318-eadb-4cff-bd54-c91d76fa8800
# ╠═e2959cd8-f3da-4a81-a539-671a59b4ddb5
# ╠═80f08a4c-ba20-4408-853e-a694df474a02
# ╠═eadf88fa-1290-4e08-9af4-53e0f613dc17
# ╠═ecfea12e-c933-49a0-b117-88658fd5723a
# ╠═b703f5a3-e869-4051-83a7-64d53f446000
# ╠═9733a667-964a-4007-9e51-f656c3838954
# ╠═572dac94-aded-44af-8855-8f2c9811f94d
# ╠═c2fad765-52d9-479e-a494-faf38736d58c
# ╠═7efc9685-1031-4882-8202-2fcd83a728d8
# ╟─1fee975e-294c-407e-ad1a-65f91c33fbd4
# ╠═ab63a35f-b674-4de1-bb6c-04f45f034a1d
# ╠═6f1a76ea-5eb7-4908-9edf-c47b7595f913
# ╠═0038f63a-d7fc-4e48-a487-6aa97df100c5
# ╠═465f5546-e9a1-4ea6-98e9-75436240dd86
# ╠═78f6195a-c9fc-4d42-9e96-58a79d4fa06e
# ╠═bc9a89bd-6b93-4ff9-8326-b24f8d3c54d1
# ╠═b4849f98-2174-4e55-9070-61b5bded79a4
# ╠═eb6cb789-6d91-4f96-bedb-b92ba5d1d69a
# ╠═2896b900-df7c-4fee-a62c-8ea64f9ffddc
# ╠═1cd28cd7-c3d7-4599-9f7b-b1d68bc094a0
