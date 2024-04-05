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

# ╔═╡ cb6e988a-f263-11ee-1f3f-53192cebcad4
begin
	using Pkg
	Pkg.activate(".")
	Pkg.develop("GridShielding")
	using GridShielding
	using OrdinaryDiffEq
	using TaylorIntegration
	using Plots
	using PlutoUI
	using Distributions
	using ProgressLogging
	using StaticArrays
	include("Shared Code/FlatUI.jl")
end

# ╔═╡ 2767663f-3ef8-44f5-81a2-8e480158266e
md"""
# Cart Pole Problem
"""

# ╔═╡ 3115801b-0a07-4a44-a6b3-d1ab2b9c0775
md"""
## Preliminaries
"""

# ╔═╡ cd2df9dc-af72-4b37-b1ef-ff8a0dcb9e0f
Pkg.add("StaticArrays")

# ╔═╡ a8aff15c-255d-498f-97dd-4c9c953ec662
begin
	function multi_field(names, types, defaults=nothing)

		if defaults == nothing
			defaults = zeros(length(names))
		end
		return PlutoUI.combine() do Field
			fields = []
	
			for (n, t, d) in zip(names, types, defaults)
				field = "`Unsupported type`"
				if t<:Number
					field = md" $n = $(Field(n, NumberField(-10000.:0.0001:10000., default=d)))"
				elseif t == String
					field = md" $n = $(Field(n, TextField(80), default=d))"
				end
				push!(fields, field)
			end
			
			md"$fields"
		end
	end

	function multi_field(typ::Type)
		multi_field(fieldnames(typ), fieldtypes(typ))
	end

	function multi_field(elem)
		defaults = [getfield(elem, f) for f in fieldnames(typeof(elem))]
		names = fieldnames(typeof(elem))
		types = fieldtypes(typeof(elem))
		multi_field(names, types, defaults)
	end
end

# ╔═╡ 7dd5c185-2b95-4297-a36c-4e2ca38952ea
md"""
## Simulating Cartpole

With $x$ being the position, $\theta$ being the angle and $a$ being the action, the state vector is

$(x, \dot x, \theta, \dot \theta,  a)^\top$



And the (frictionless) system is controlled by [1]

$\displaystyle\ddot{\theta}=\frac{{{g} \sin{\theta}+ \cos{\theta}\cdot{\left(\frac{{-{F}-{m}_{{p}}{l}\dot{\theta}^{2} \sin{\theta}}}{{{m}_{{c}}+{m}_{{p}}}}\right)}}}{{{l}{\left(\frac{4}{{3}}-\frac{{{m}_{{p}}{{\cos}^{2}\theta}}}{{{m}_{{c}}+{m}_{{p}}}}\right)}}}$

and 

$\displaystyle\ddot{{x}}=\frac{{{F}+{m}_{{p}}{l}{\left(\dot{\theta}^{2} \sin{\theta}-\ddot{\theta} \cos{\theta}\right)}}}{{{m}_{{c}}+{m}_{{p}}}}$


Cart-pole code is adapted from [2].

[1] R. V. Florian, “Correct equations for the dynamics of the cart-pole system,” 2005.

[2] C. Schilling, A. Lukina, E. Demirović, and K. Larsen, “Safety Verification of Decision-Tree Policies in Continuous Time,” Advances in Neural Information Processing Systems, vol. 36, pp. 14750–14769, Dec. 2023.

"""

# ╔═╡ 86ab3957-129c-4c71-911d-a54fd80cfd2d
begin
	struct Action
	    name::String
	    value::Int
	end
	
	function Base.isless(A1::Action, A2::Action)
	    return A1.name < A2.name
	end

	function Base.zero(::Action)
		Action("", 0)
	end
	
	
	# The task of `apply(A, X)` is to update the control mode of `X` to `A.value`.
	# By convention, the control mode is stored in the last dimension.
	
	
	function apply(A::Action, X)
	    return apply!(A, deepcopy(X))
	end
	
	function apply(A::Action, X::NTuple{N,T}) where {N, T<:Number}
	    return Tuple(i == N ? T(A.value) : x for (i, x) in enumerate(X))
	end
	
	
	function apply!(A::Action, x::Vector)
	    x[end] = A.value
	    return x
	end
end

# ╔═╡ b52604cc-e8bc-4b53-84ad-79cf019c1667
const left = Action("left", 0); const right = Action("right", 1);

# ╔═╡ e1f217af-759a-4868-b0b2-6ce08de324ea
const left_value = left.value, const right_value = right.value

# ╔═╡ 274e6f36-bcae-4232-9bb6-732f954aa4e5
const actions = [left, right]

# ╔═╡ 788e089e-ded6-4ad7-9951-1d10a32a8295
function get_action(action_value::Int)
	for action in actions
		if action.value == action_value
			return action
		end
	end
	error("Unexpected value: $((;action_value))")
end

# ╔═╡ 35605d87-4c3a-49a9-93d1-a5fceede3653
@kwdef struct CartPoleMechanics
	gravity = 9.8
	cart_mass = 1.0
	pole_mass = 0.1
	total_mass = pole_mass + cart_mass
	pole_length = 0.5  # actually half the pole's length
	polemass_length = pole_mass * pole_length
	force_mag = 10.0
	τ = 0.02  # control cycle length
end

# ╔═╡ 822b3f06-613a-4992-8baf-6450a405d961
const m = CartPoleMechanics()

# ╔═╡ d385483b-df9c-4e19-9071-815df436f7bc
const (;gravity, cart_mass, pole_mass, total_mass, pole_length, polemass_length, force_mag) = m

# ╔═╡ 0b7509e7-433b-41c7-a971-bfdb164c44a1
function cartpole!(ds, s, p, t)
	
	# x = s[1]  # cart position
    x_vel = s[2]  # cart velocity
    θ = s[3]  # pendulum angle
    θ_vel = s[4]  # pendulum angular velocity
    action = s[5]  # action
    if action == left.value
        force = -force_mag
	elseif action == right.value
        force = force_mag
	else
		error("Unexpected action value $action")
    end
	
    cosθ = cos(θ)
    sinθ = sin(θ)

    temp = (force + polemass_length * θ_vel^2 * sinθ) / total_mass
    θ_acc = (gravity * sinθ - cosθ * temp) / (pole_length * (4.0/3.0 - pole_mass * cosθ^2 / total_mass))
    x_acc = temp - polemass_length * θ_acc * cosθ / total_mass

    ds[1] = x_vel
    ds[2] = x_acc
    ds[3] = θ_vel
    ds[4] = θ_acc
    ds[5] = zero(action)
	nothing
end

# ╔═╡ 8d4c5990-4fc1-47cc-a10b-18682331957e
const time_test_ds, time_test_s = 
	Float64[0, 0, 0, 0, 0], Float64[0, 0, 0, 0, left.value]

# ╔═╡ 1d690e0a-2748-4e52-b23e-66cac2e83411
@time cartpole!(time_test_ds, time_test_s, :_, m.τ)

# ╔═╡ bc682eeb-0120-427e-9bfb-b306a7d82d0c
s0() = Float64[0, 0, rand(Uniform(-0.1, 0.1)), 0, left.value]

# ╔═╡ 0a67aaa8-189b-4e19-be48-6693b05bc54d
# Terminal state
st = Float64[-1, -1, -1, -1, -1]

# ╔═╡ 2f21c16b-1dc5-4960-a93d-18d13ce01f0e
solver = TaylorMethod(8)

# ╔═╡ 78bfb1a5-0f9c-4e8b-ae9b-a149f2acd953
solve(ODEProblem(cartpole!, s0(), m.τ), solver)

# ╔═╡ 187344d0-b0c1-4f63-b257-9254f5dab869
function simulate_point(m::CartPoleMechanics, s0, a::Action)
	prob = ODEProblem(cartpole!, apply(a, s0), m.τ)
	result = solve(prob, solver)
	result[2]
end

# ╔═╡ c97daca9-3e71-4c38-af24-19a40d326467
CartPoleState = Vector{Float64}

# ╔═╡ 5306c76f-f0df-427c-88d3-a738f9048721
@kwdef struct CartPoleTrace
	times=Float64[]
	actions=Action[]
	states=CartPoleState[]
end

# ╔═╡ b345bb20-15ee-4db5-865e-d99952a64f02
random_policy = _ -> rand((left, right))

# ╔═╡ a30b08af-8584-481e-930c-793db7ee97d4
picked_right_last = false

# ╔═╡ 5e1614e9-1dca-4a32-b169-5200d708396c
alternating_policy = _ -> if picked_right_last
	picked_right_last = false
	left
else
	picked_right_last = true
	right
end

# ╔═╡ e7108be6-8280-4591-b8e3-3e6d4207aed9
function cart(x)
	w = 0.2 # Half width actually
	h = 0.1 # same with height
	Shape([x - w, x - w, x + w, x + w],
		  [-h, +h, +h, -h])
end

# ╔═╡ ed4eb3ec-3947-49c6-b1a0-d0bcafb281ff
cart(0)

# ╔═╡ 251427d4-0aae-4a1f-a31f-71832877c749
function animate_sequence(trace::CartPoleTrace; speed=1)
	xs = [s[1] for s in trace.states]
	xlims = max(abs(min(xs...)), max(xs...), 1)
	xlims = (-xlims, xlims)
	🎥 = @animate for (i, s) in enumerate(trace.states)
		x = s[1]
		θ = π/2 - s[3]
		
		plot(;
			xlims,
			ylims=(-m.pole_length, m.pole_length*3),
			yticks=nothing,
			aspectratio=:equal)

		# Action #
		a = trace.actions[i]
		if a.value == left.value
			plot!([x, x - 0.3], [0, 0],
				label=nothing,
				marker=:ltriangle,
				markerstrokewidth=0,
				color=colors.EMERALD,)
		elseif a.value == right.value
				plot!([x, x + 0.3], [0, 0],
				label=nothing,
				marker=:rtriangle,
				markerstrokewidth=0,
				color=colors.EMERALD,)
		else
			error("Unexpected action value $a")
		end
		
		# Cart #
		plot!(cart(x), 
			color=colors.WET_ASPHALT,
			linewidth=0,
			label=nothing)

		# Pole #
		pole_start = (x, 0)

		pole_end = (x + 2*m.pole_length*cos(θ), 2*m.pole_length*sin(θ))

		plot!([pole_start, pole_end],
			color=colors.ASBESTOS,
			label=nothing,
			linewidth=4)
	end
	gif(🎥, fps=1/m.τ*speed, show_msg=false)
end

# ╔═╡ 7f1f08f2-bdf7-4563-887c-a8e380549f94
md"""
## Q-learning for fun
"""

# ╔═╡ 1ebddb2c-428e-4bca-bbe1-f5f2189b5418
Q_granularity = [0.5, 0.5, 0.05, 0.2]

# ╔═╡ 0f9da8fe-7a26-4b21-9a64-906e908d835e
Q_bounds = Bounds([-2.2, -4.8, -0.4, -0.836],
	[2.2, 4.8, 0.4, 0.836])

# ╔═╡ 5cfb472e-8d08-48a7-bd27-c23f67a065d5
# reward
function r(s)
	if s ∉ Q_bounds
		-10
	else
		-abs(s[3])
	end
end

# ╔═╡ 05cdc837-58e1-4112-922a-e8344bc4ee66
get_size(Q_granularity, Q_bounds), get_size(Q_granularity, Q_bounds) |> prod

# ╔═╡ 0a9c51e1-0cc3-4bb0-9fa5-0f3962aae605
@bind episodes NumberField(0:typemax(Int64), default=5)

# ╔═╡ 398af1f5-5ffb-4667-b8cc-ee1920de2997
@bind γ NumberField(0.0001:0.0001:1, default=0.9)

# ╔═╡ 4cb85ba8-86f6-436c-9479-8e162c8b7d54
@bind ϵ_base NumberField(0.0001:0.0001:1, default=0.8)

# ╔═╡ 5d2064c4-a987-4296-ae30-fed483057eff
@bind α_base NumberField(0.0001:0.0001:1, default=0.01)

# ╔═╡ 2ef0dc55-796c-4c5f-89be-96872cbc3c50
function α(t; episodes=episodes)
	if t < episodes/2
		α_base
	else
		α_base/(1 + 0.2*(t - episodes/2))
	end
end

# ╔═╡ 742280a3-8f06-49de-aea6-79462cfeb2f8
function ϵ(t; episodes=episodes)
	if t < episodes/2
		ϵ_base
	else
		ϵ_base/(1 + 0.2*(t - episodes/2))
	end
end

# ╔═╡ 534b9c41-52dd-43c5-bb5d-f44dcde0f70f
# ϵ-greedy choice from Q.
function ϵ_greedy(ϵ::Number, Q, s)
	if rand(Uniform(0, 1)) < ϵ
		return rand((left, right))
	else
		return argmax((a) -> get_value(box(Q[a], s)), (left, right))
	end
end

# ╔═╡ 33f338e1-4792-4a68-802c-c814a6c89fd6
begin
	p1 = plot(xlabel="t")
	
	plot!(y -> ϵ(y; episodes), 
		xlim=(0, episodes), 
		label="ϵ", 
		color=colors.ALIZARIN)
	
	hline!([0], line=:black, label=nothing)
	
	p2 = plot(xlabel="t")
	
	plot!(y -> α(y; episodes), 
		xlim=(0, episodes), 
		label="α", 
		color=colors.PETER_RIVER)
	
	hline!([0], line=:black, label=nothing)
	plot(p1, p2, size=(600, 200))
end

# ╔═╡ 57bad1a9-9cea-4e2f-903d-004bccafffb3
s0()

# ╔═╡ 6c84e1aa-f45a-453c-8a78-f3976c605385
md"""
# Making a shield
"""

# ╔═╡ 731d3595-1746-46cd-9851-e285673d6a1b
md"""
### Safety
"""

# ╔═╡ e5371840-2989-4f7f-b15c-4b07e8e96a3b
# Bounds of the state space
# The episode terminates if the cart is outside the [-2.4, 2.4] range
# or the angle is outside [-0.418, 0.418]. I assume the 
cart_pole_bounds = Bounds(
	[-2.4, -2.4/m.τ/4, -0.418, -0.418/m.τ/4],
	[2.4, 2.4/m.τ/4, 0.418, 0.418/m.τ/4])

# ╔═╡ 16060176-987c-4450-819e-9b06dc23a051
function simulate_sequence(m::CartPoleMechanics, s0, policy, duration)
	trace = CartPoleTrace()
	push!(trace.times, 0)
	push!(trace.states, s0)
	for t in m.τ:m.τ:duration
		s = trace.states[end]
		action = policy(s)
		s′ = simulate_point(m, s, action)
		push!(trace.times, t)
		push!(trace.actions, action)
		if s′ ∉ cart_pole_bounds
			break
		else
			push!(trace.states, s′)
		end
			
	end
	trace
end

# ╔═╡ 4a18c6f8-b31c-487b-8a22-a8b6ed0b3b46
trace = simulate_sequence(m, s0(), random_policy, 4)

# ╔═╡ ea84c513-b4ca-41df-96ec-1c230fde9f3d
animate_sequence(trace; speed=1)

# ╔═╡ 40550db6-a7b9-496c-b321-b61cf1239e18
function Q_learn()
	Q = Dict(left => Grid(Q_granularity, Q_bounds, data_type=Float64), 
			right => Grid(Q_granularity, Q_bounds, data_type=Float64))

	# Discourage exploration; we want to stay near s0
	for partition in Q[left]
		set_value!(partition, -1)
	end
	for partition in Q[right]
		set_value!(partition, -1)
	end
	
	@progress for i ∈ 1:episodes
		Sₜ = s0()
		Aₜ = rand((left, right))
		for t ∈ 0:m.τ:10
			Sₜ₊₁ = simulate_point(m, Sₜ, Aₜ)
			if Sₜ₊₁ ∉ cart_pole_bounds break end
			if Sₜ₊₁ ∉ Q_bounds continue end
			Q_Sₜ_Aₜ = box(Q[Aₜ], Sₜ)
			set_value!(Q_Sₜ_Aₜ, 
				get_value(Q_Sₜ_Aₜ) + 
				α(t)*(r(Sₜ) + γ*max([get_value(box(Q[a′], Sₜ₊₁)) 
				for a′ in (left, right)]...) - get_value(Q_Sₜ_Aₜ)))
			
			Aₜ₊₁ = ϵ_greedy(ϵ(t), Q, Sₜ)
			Sₜ, Aₜ = Sₜ₊₁, Aₜ₊₁
		end
	end

	return Q
end

# ╔═╡ 54afb933-b75d-4fb5-8b26-396c123f09ca
Q = Q_learn()

# ╔═╡ 787677e2-bdf3-43e3-ac71-563a483ef8dc
Q_policy = s -> begin
	if s ∉ Q[left] 
		return rand((left, right))
	end
	if get_value(box(Q[left], s)) > get_value(box(Q[right], s))
		return left
	else
		return right
	end
end

# ╔═╡ 78c4dab0-3d95-462a-a214-9578f84b6cb4
Q_trace = simulate_sequence(m, s0(), Q_policy, 4)

# ╔═╡ 2ac120e4-f380-4b02-bc7f-a1d5e84d7c36
animate_sequence(Q_trace)

# ╔═╡ d712571e-ced8-4f06-8b44-6874fcd3e15d
length(Q[left].array |> unique),
length(Q[right].array |> unique)

# ╔═╡ 453f70e4-c22f-4c75-b786-d523c8e4bf9c
begin
	function is_safe(s::Vector{Float64})
		s ∈ cart_pole_bounds
	end

	function is_safe(bounds::Bounds)
		for s in SupportingPoints([2, 1, 2, 1, 1], bounds)
			if !is_safe([s..., 0])
				return false
			end
		end
		return true
	end
end

# ╔═╡ 16315ce1-9dca-4284-936b-32a204b56108
is_safe(s0())

# ╔═╡ e0892cae-9ef0-4e57-9a1c-91bf34043956
md"""
### The Grid
"""

# ╔═╡ 668f4592-75fd-445e-a0fa-56ee02a03f2d
no_action = actions_to_int([])

# ╔═╡ 611ba5df-6af9-413b-8e8a-b3da0c825d3e
any_action = actions_to_int([left.value, right.value])

# ╔═╡ 26cfc8ec-1351-468f-b9dc-e76acec6e777
function initializer(partition)
	is_safe(partition) ? any_action : no_action
end

# ╔═╡ 5dde6492-564f-46cb-848d-8a28ea2adb5f
granularity = Float64[0.2, 2.0, 0.04, 0.4]

# ╔═╡ bcbf4a16-ce8f-451e-b58b-0bf9d8d0d872
get_size(granularity, cart_pole_bounds)

# ╔═╡ 38bc7025-9e1f-4101-a53d-a3a7ff802aa7
grid_bounds = Bounds(
	cart_pole_bounds.lower .- granularity, 
	cart_pole_bounds.upper .+ granularity)

# ╔═╡ 0610d08b-020e-4ec8-9815-1d0a4c592899
get_size(granularity, cart_pole_bounds)

# ╔═╡ b966dc17-050f-40dc-adee-b6f9e79b4b0c
begin
	grid = Grid(granularity, grid_bounds)
	GridShielding.initialize!(grid, initializer)
end

# ╔═╡ f723aa48-e30b-4666-ad70-c20ae10fb4bb
is_safe(Bounds(box(grid, grid.bounds.lower)))

# ╔═╡ d618a32d-370b-4b45-a691-8cd3ce09eff3
length(grid)

# ╔═╡ 3fd479d1-c43a-4c6f-95f8-0c74a9ffbf18
begin
	shieldcolors = [colors.WET_ASPHALT, colors.AMETHYST, colors.SUNFLOWER, colors.CLOUDS]
	shieldlabels = [a for a in 0:3]
	shieldlabels = [[a for a in int_to_actions(Int, a)] for a in shieldlabels]
	shieldlabels = [[get_action(aa).name for aa in a] for a in shieldlabels]
	shieldlabels = [join(a, "," ) for a in shieldlabels]
	shieldlabels = ["{$a}" for a in shieldlabels]

	zip(shieldcolors, shieldlabels) |> collect
end

# ╔═╡ 6de525db-e339-435f-9f87-620fed817839
@bind s_input multi_field(["x", "x_vel", "θ", "θ_vel"], 
	[Float64, Float64, Float64, Float64])

# ╔═╡ 39c77465-5d24-4d5f-808f-90b18ac18446
@bind action Select([a => a.name for a in (left, right)])

# ╔═╡ a5a3b815-23c4-4acd-90d5-a12c721e7866
action

# ╔═╡ 1b882558-e83e-4679-8d51-3dc54040cdf1
s = apply(action, [s_input.x, s_input.x_vel, s_input.θ, s_input.θ_vel, 0])

# ╔═╡ 65f2972d-d1c1-4e21-95ed-f73f2047093d
@time simulate_point(m, s, action)

# ╔═╡ 915e0813-216c-4f01-8341-75c57198dc44
partition = box(grid, s)

# ╔═╡ ee7fdda2-732f-42be-926c-bcdbf0634299
let
	slice = Any[partition.indices...]
	slice[1] = slice[2] = Colon()
	draw(grid, slice,
		colors=shieldcolors,
		color_labels=shieldlabels,
		xlabel="x",
		ylabel="θ")
end

# ╔═╡ 1d8867a6-a0b9-4ee1-8cf0-dabfa4678937
bounds = Bounds(partition)

# ╔═╡ f4c3e866-50aa-440b-a141-65de2daf9c4c
is_safe(bounds)

# ╔═╡ ace4eb92-b880-4fe6-9391-ad5bc586b802
md"""
### Reachability
"""

# ╔═╡ dee84d76-f3a2-45a6-b1db-3ca865877de1
samples_per_axis = [3, 3, 3, 3]

# ╔═╡ 5376f447-716b-4011-bd3a-4b60db5ed110
function reachability_function(partition, action)::Vector{Vector{Int64}}
	result = Vector{Int64}[]
	grid = partition.grid
	for s in SupportingPoints(samples_per_axis, partition)
		if s ∉ Bounds(partition) continue end
		#for r in SupportingPoints(samples_per_random_axis, Bounds((-1,), (1,)))
		s′ = simulate_point(m, [s..., 0.0], action)
		clamp(s′, grid.bounds)
		partition′ = box(grid, s′)
		if partition′.indices ∈ result
			continue
		end
		push!(result, partition′.indices)
		#end
	end
	result
end

# ╔═╡ b1de4c47-e90e-45c3-8c60-340516b42f8e
reachability_function(partition, action)

# ╔═╡ 5c6c1a0d-1442-4d39-bc4d-f7c943e14d97
md"""
### Mainmatter
"""

# ╔═╡ b5453b51-a878-4073-8152-c69d85d30ec1
# ╠═╡ disabled = true
#=╠═╡
reachability_function_precomputed = get_transitions(reachability_function, 
	[0, 1], 
	grid);
  ╠═╡ =#

# ╔═╡ e77abd9a-de23-40b8-a442-b0971339f903
#=╠═╡
reachability_function_precomputed[left][partition.indices]
  ╠═╡ =#

# ╔═╡ b79619f1-aa55-4a5c-851f-7387d411d8eb
@bind max_steps NumberField(0:1000)

# ╔═╡ 0c28089a-1547-47f4-a411-e3a57cac6a6d
#=╠═╡
shield, max_steps_reached = 
	make_shield(reachability_function_precomputed, actions, grid; max_steps)
  ╠═╡ =#

# ╔═╡ c7a4e65c-a907-468e-b31c-ce05393d41d5
#=╠═╡
let
	slice = Any[partition.indices...]
	slice[1] = slice[2] = Colon()
	draw(shield, slice,
		colors=shieldcolors,
		color_labels=shieldlabels,
		xlabel="x",
		ylabel="θ")
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─2767663f-3ef8-44f5-81a2-8e480158266e
# ╟─3115801b-0a07-4a44-a6b3-d1ab2b9c0775
# ╠═cb6e988a-f263-11ee-1f3f-53192cebcad4
# ╠═cd2df9dc-af72-4b37-b1ef-ff8a0dcb9e0f
# ╟─a8aff15c-255d-498f-97dd-4c9c953ec662
# ╟─7dd5c185-2b95-4297-a36c-4e2ca38952ea
# ╠═b52604cc-e8bc-4b53-84ad-79cf019c1667
# ╠═e1f217af-759a-4868-b0b2-6ce08de324ea
# ╠═274e6f36-bcae-4232-9bb6-732f954aa4e5
# ╠═788e089e-ded6-4ad7-9951-1d10a32a8295
# ╠═86ab3957-129c-4c71-911d-a54fd80cfd2d
# ╠═35605d87-4c3a-49a9-93d1-a5fceede3653
# ╠═822b3f06-613a-4992-8baf-6450a405d961
# ╠═d385483b-df9c-4e19-9071-815df436f7bc
# ╠═0b7509e7-433b-41c7-a971-bfdb164c44a1
# ╠═8d4c5990-4fc1-47cc-a10b-18682331957e
# ╠═1d690e0a-2748-4e52-b23e-66cac2e83411
# ╠═bc682eeb-0120-427e-9bfb-b306a7d82d0c
# ╠═0a67aaa8-189b-4e19-be48-6693b05bc54d
# ╠═2f21c16b-1dc5-4960-a93d-18d13ce01f0e
# ╠═78bfb1a5-0f9c-4e8b-ae9b-a149f2acd953
# ╠═187344d0-b0c1-4f63-b257-9254f5dab869
# ╠═65f2972d-d1c1-4e21-95ed-f73f2047093d
# ╠═c97daca9-3e71-4c38-af24-19a40d326467
# ╠═5306c76f-f0df-427c-88d3-a738f9048721
# ╠═16060176-987c-4450-819e-9b06dc23a051
# ╠═b345bb20-15ee-4db5-865e-d99952a64f02
# ╠═a30b08af-8584-481e-930c-793db7ee97d4
# ╠═5e1614e9-1dca-4a32-b169-5200d708396c
# ╠═e7108be6-8280-4591-b8e3-3e6d4207aed9
# ╠═ed4eb3ec-3947-49c6-b1a0-d0bcafb281ff
# ╟─251427d4-0aae-4a1f-a31f-71832877c749
# ╠═4a18c6f8-b31c-487b-8a22-a8b6ed0b3b46
# ╠═ea84c513-b4ca-41df-96ec-1c230fde9f3d
# ╟─7f1f08f2-bdf7-4563-887c-a8e380549f94
# ╠═bcbf4a16-ce8f-451e-b58b-0bf9d8d0d872
# ╠═5cfb472e-8d08-48a7-bd27-c23f67a065d5
# ╠═1ebddb2c-428e-4bca-bbe1-f5f2189b5418
# ╠═0f9da8fe-7a26-4b21-9a64-906e908d835e
# ╠═05cdc837-58e1-4112-922a-e8344bc4ee66
# ╠═40550db6-a7b9-496c-b321-b61cf1239e18
# ╠═0a9c51e1-0cc3-4bb0-9fa5-0f3962aae605
# ╠═398af1f5-5ffb-4667-b8cc-ee1920de2997
# ╠═4cb85ba8-86f6-436c-9479-8e162c8b7d54
# ╠═5d2064c4-a987-4296-ae30-fed483057eff
# ╠═2ef0dc55-796c-4c5f-89be-96872cbc3c50
# ╠═742280a3-8f06-49de-aea6-79462cfeb2f8
# ╠═534b9c41-52dd-43c5-bb5d-f44dcde0f70f
# ╟─33f338e1-4792-4a68-802c-c814a6c89fd6
# ╠═57bad1a9-9cea-4e2f-903d-004bccafffb3
# ╠═a5a3b815-23c4-4acd-90d5-a12c721e7866
# ╠═54afb933-b75d-4fb5-8b26-396c123f09ca
# ╠═787677e2-bdf3-43e3-ac71-563a483ef8dc
# ╠═d712571e-ced8-4f06-8b44-6874fcd3e15d
# ╠═78c4dab0-3d95-462a-a214-9578f84b6cb4
# ╠═2ac120e4-f380-4b02-bc7f-a1d5e84d7c36
# ╟─6c84e1aa-f45a-453c-8a78-f3976c605385
# ╟─731d3595-1746-46cd-9851-e285673d6a1b
# ╠═e5371840-2989-4f7f-b15c-4b07e8e96a3b
# ╠═453f70e4-c22f-4c75-b786-d523c8e4bf9c
# ╠═16315ce1-9dca-4284-936b-32a204b56108
# ╠═f4c3e866-50aa-440b-a141-65de2daf9c4c
# ╠═f723aa48-e30b-4666-ad70-c20ae10fb4bb
# ╟─e0892cae-9ef0-4e57-9a1c-91bf34043956
# ╠═668f4592-75fd-445e-a0fa-56ee02a03f2d
# ╠═611ba5df-6af9-413b-8e8a-b3da0c825d3e
# ╠═26cfc8ec-1351-468f-b9dc-e76acec6e777
# ╠═5dde6492-564f-46cb-848d-8a28ea2adb5f
# ╠═38bc7025-9e1f-4101-a53d-a3a7ff802aa7
# ╠═0610d08b-020e-4ec8-9815-1d0a4c592899
# ╠═b966dc17-050f-40dc-adee-b6f9e79b4b0c
# ╠═d618a32d-370b-4b45-a691-8cd3ce09eff3
# ╠═3fd479d1-c43a-4c6f-95f8-0c74a9ffbf18
# ╟─ee7fdda2-732f-42be-926c-bcdbf0634299
# ╠═6de525db-e339-435f-9f87-620fed817839
# ╠═39c77465-5d24-4d5f-808f-90b18ac18446
# ╠═1b882558-e83e-4679-8d51-3dc54040cdf1
# ╠═915e0813-216c-4f01-8341-75c57198dc44
# ╠═1d8867a6-a0b9-4ee1-8cf0-dabfa4678937
# ╟─ace4eb92-b880-4fe6-9391-ad5bc586b802
# ╠═dee84d76-f3a2-45a6-b1db-3ca865877de1
# ╠═5376f447-716b-4011-bd3a-4b60db5ed110
# ╠═b1de4c47-e90e-45c3-8c60-340516b42f8e
# ╟─5c6c1a0d-1442-4d39-bc4d-f7c943e14d97
# ╠═b5453b51-a878-4073-8152-c69d85d30ec1
# ╠═e77abd9a-de23-40b8-a442-b0971339f903
# ╠═b79619f1-aa55-4a5c-851f-7387d411d8eb
# ╠═0c28089a-1547-47f4-a411-e3a57cac6a6d
# ╠═c7a4e65c-a907-468e-b31c-ce05393d41d5
