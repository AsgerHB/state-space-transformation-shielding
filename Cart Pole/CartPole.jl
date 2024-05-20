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
	Pkg.activate("..")
	Pkg.develop("GridShielding")
	using GridShielding
	using OrdinaryDiffEq
	using TaylorIntegration
	using Plots
	using PlutoUI
	using Distributions
	using ProgressLogging
	using StaticArrays
	using Polynomials
	using Unzip
	include("../Shared Code/FlatUI.jl")
end

# ╔═╡ 5c808e32-0aa6-49ec-a7b1-d01e44aa48ea
using Measures

# ╔═╡ 2767663f-3ef8-44f5-81a2-8e480158266e
md"""
# Cart Pole Problem
"""

# ╔═╡ 3115801b-0a07-4a44-a6b3-d1ab2b9c0775
md"""
## Preliminaries
"""

# ╔═╡ cd2df9dc-af72-4b37-b1ef-ff8a0dcb9e0f
TableOfContents()

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

# ╔═╡ 319eff14-0f2e-458e-99b3-b3a64eccccee
md"""
### 🛠 `make_paper_friendly_figures` 
Toggle: $(@bind make_paper_friendly_figures CheckBox(default=false))
"""

# ╔═╡ 1d5aced4-a27a-4f3b-ab08-69b1a3b1559f
begin
	default_font = default(:fontfamily)
	default_size = default(:size)
	default_margin = default(:margin)
end;

# ╔═╡ 35c76906-f2dd-4f4d-af43-0fafe69d211b
begin
	paper_font = "serif-roman" 	# https://gr-framework.org/fonts.html
	paper_size = (300, 220)
	paper_margin = 0mm
end;

# ╔═╡ 8172cdb6-b9a4-4afc-844f-245e5d951bb7
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

# ╔═╡ e008fe83-3851-4061-9e5a-35e938a53ec1
theme_type; plot(rand(1:10, 10), xlabel="θαβ ∪⋆", ylabel="asdf")

# ╔═╡ a8eefbb3-a09d-4915-9b9a-76f896636dcc
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
	
	function apply!(A::Action, x::MVector)
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

# ╔═╡ 0a67aaa8-189b-4e19-be48-6693b05bc54d
# Terminal state
st = Float64[-1, -1, -1, -1, -1]

# ╔═╡ 2f21c16b-1dc5-4960-a93d-18d13ce01f0e
solver = Tsit5()

# ╔═╡ c97daca9-3e71-4c38-af24-19a40d326467
const CartPoleState = MVector{5, Float64}

# ╔═╡ bc682eeb-0120-427e-9bfb-b306a7d82d0c
s0()::CartPoleState = Float64[0, 0, rand(Uniform(-0.01, 0.01)), 0, left.value]

# ╔═╡ 78bfb1a5-0f9c-4e8b-ae9b-a149f2acd953
solve(ODEProblem(cartpole!, s0(), m.τ), solver)

# ╔═╡ 187344d0-b0c1-4f63-b257-9254f5dab869
function simulate_point(m::CartPoleMechanics, s0, a::Action)::CartPoleState
	prob = ODEProblem(cartpole!, apply!(a, s0), m.τ)
	result = solve(prob, solver)
	result[end]
end

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

# ╔═╡ 7f1f08f2-bdf7-4563-887c-a8e380549f94
md"""
## Q-learning for fun
"""

# ╔═╡ 1ebddb2c-428e-4bca-bbe1-f5f2189b5418
Q_granularity = Float64[0.5, 1, 0.1, 0.2]

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

# ╔═╡ 398af1f5-5ffb-4667-b8cc-ee1920de2997
@bind γ NumberField(0.0001:0.0001:1, default=0.9)

# ╔═╡ 4cb85ba8-86f6-436c-9479-8e162c8b7d54
@bind ϵ_base NumberField(0.0001:0.0001:1, default=0.8)

# ╔═╡ 5d2064c4-a987-4296-ae30-fed483057eff
@bind α_base NumberField(0.0001:0.0001:1, default=0.01)

# ╔═╡ 534b9c41-52dd-43c5-bb5d-f44dcde0f70f
# ϵ-greedy choice from Q.
function ϵ_greedy(ϵ::Number, Q, s)
	if rand(Uniform(0, 1)) < ϵ
		return rand((left, right))
	else
		return argmax((a) -> get_value(box(Q[a], s)), (left, right))
	end
end

# ╔═╡ 57bad1a9-9cea-4e2f-903d-004bccafffb3
s0()

# ╔═╡ 6c84e1aa-f45a-453c-8a78-f3976c605385
md"""
# Making a shield
"""

# ╔═╡ a7296808-4eb8-4f10-8683-adc4963b21ce
md"""
## Altered State Space
"""

# ╔═╡ 9da1c57f-c922-420c-8c46-0d6842553f8c
md"""
### 🛠 Toggle altered state space

`enable_altered_state_space =` $(@bind enable_altered_state_space CheckBox(default=false))
"""

# ╔═╡ 338e5d40-6251-429c-9b0d-ef92460a7e52
const AlteredState = MVector{4, Float64}

# ╔═╡ 3aeb0922-e5f8-4311-b66d-dd42d61f18f3
altered_state_axes, state_axes_simple = if enable_altered_state_space
	([ "\$x\$", "\$p(x, v)\$",
	  "\$θ\$", "\$p(\\theta, \\omega)\$", ],
	[ "cart_pos", "poly1",
	  "pole_ang", "poly2", ])
else
	([ "\$x\$", "\$v\$",
	  "\$θ\$", "\$ω\$", ],
	[ "cart_pos", "cart_vel",
	  "pole_ang", "pole_vel", ])
end

# ╔═╡ b4ddbac3-33e9-42e7-bac7-a2ec32093678
function P1(x, x_vel)
	# Learned for x-safety only (high-resolution shield)
	# 1st degree polynomial.
	x_vel - (1.1576428742879741e-16 - 1.3323889308715684*x) 
end

# ╔═╡ 64be4f5a-97f6-49bc-b848-450506d9ceb1
P1⁻¹(θ, P1_s) = P1_s - P1(θ, 0)

# ╔═╡ 5645d7a7-0f23-4a02-ab30-49a6cda9ce17
# Naming confusion: p1, p2 ... are plots. P1, P2 are polynomials
function P2(θ, θ_vel)

	# Learned from 3 steps of a 20x20 shield computation
	return θ_vel - (- 4.550831135117032*θ - 141.6953270125445*θ^3)

	# Learned from a fully compuetd 100x100 shield
	θ_vel - (-1.7440907420745881e-16 - 5.829601684910079*θ - 7.54337188727965e-16*θ^2 - 26.782941687007654*θ^3)
end

# ╔═╡ d90ea6c4-316d-46b1-8688-23d95f3cca61
P2⁻¹(θ, P2_s) = P2_s - P2(θ, 0)

# ╔═╡ 6eedcd3c-1aa5-4d8b-9e5f-8465e608f6c7
if enable_altered_state_space

	function f(s::CartPoleState)::AlteredState
		x = s[1]
		x_vel = s[2]
		θ = s[3]
		θ_vel = s[4]
		AlteredState(x, P1(x, x_vel), θ, P2(θ, θ_vel))
	end

else
	
	function f(s::CartPoleState)::AlteredState
		x = s[1]
		x_vel = s[2]
		θ = s[3]
		θ_vel = s[4]
		AlteredState(x, x_vel, θ, θ_vel)
	end
	
end

# ╔═╡ 1eb460d5-3de2-44b0-9a75-378a35e8cf6f
if enable_altered_state_space

	function f⁻¹(s::AlteredState, samples=4)::Vector{CartPoleState}
		x = s[1]
		P1_s = s[2]
		θ = s[3]
		P2_s = s[4]
		#= [
			CartPoleState(x, x_vel_1, θ, θ_vel, 0),
			CartPoleState(x, x_vel_2, θ, θ_vel, 0),
		] =#
		[CartPoleState(x, P1⁻¹(x, P1_s), θ, P2⁻¹(θ, P2_s), 0)]
	end

else

	function f⁻¹(s::AlteredState, samples=4)::Vector{CartPoleState}
		x = s[1]
		x_vel = s[2]
		θ = s[3]
		θ_vel = s[4]
		[CartPoleState(x, x_vel, θ, θ_vel, 0)]
	end
	
end

# ╔═╡ f52cb596-8f54-4721-879c-65d8f206a224
s0

# ╔═╡ edfae222-4bd5-4411-bfa4-a6e9fc5ff725
# ╠═╡ disabled = true
#=╠═╡
let
	fs = f(s)
	f⁻¹s = f⁻¹(fs, 10)
	
	scatter([s[3]],
		xlim=(cart_pole_bounds.lower[3], cart_pole_bounds.upper[3]),
		ylim=(cart_pole_bounds.lower[4], cart_pole_bounds.upper[4]),
		[s[4]],
		label="s",
		markersize=6,
		markerstrokewidth=0,
		color=colors.EMERALD)

	plot!([], [],
		linewidth=0,
		label="f(s)=$(fs[1])")
	
	scatter!([s[3] for s in f⁻¹s],
		[s[4] for s in f⁻¹s],
		label="f⁻¹(f(s))",
		markersize=3,
		markerstrokewidth=0,
		color=colors.PETER_RIVER)
end
  ╠═╡ =#

# ╔═╡ 731d3595-1746-46cd-9851-e285673d6a1b
md"""
## Safety
"""

# ╔═╡ 2f2a6934-b633-47f1-9a26-88b7ece30f46
md"""
### 🛠 safety constraints
`concerned_with_position =` $(@bind concerned_with_position CheckBox(default=false))

`concerned_with_angle =` $(@bind concerned_with_angle CheckBox(default=true))
"""

# ╔═╡ e5371840-2989-4f7f-b15c-4b07e8e96a3b
# Bounds of the state space
# The episode terminates if the cart is outside the [-2.4, 2.4] range
# or the angle is outside [-0.418, 0.418].
# The theoretical upper bound for the velocities are the values 
# that will cause a violation in the next step.
# However, these are hideously large, so I've done some empirical observations
# and set my own bounds
const cart_pole_bounds = Bounds(
	[-2.4, -6, -0.418, -3],
	[2.4, 6, 0.418, 3])

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

		# Wrap back if leaving frame
		if !concerned_with_position && s′[1] < cart_pole_bounds.lower[1]
			s′[1] = cart_pole_bounds.upper[1]
			push!(trace.states, s′)
		elseif !concerned_with_position && s′[1] >= cart_pole_bounds.upper[1]
			s′[1] = cart_pole_bounds.lower[1]
			push!(trace.states, s′)
		
			
		# Common case
		else
			push!(trace.states, s′)
		end
			
	end
	trace
end

# ╔═╡ 4a18c6f8-b31c-487b-8a22-a8b6ed0b3b46
trace = simulate_sequence(m, s0(), random_policy, 4)

# ╔═╡ 58458e09-c9cf-4d75-b153-17d807b783f3
begin
	function plot_state!(state::CartPoleState, action::Union{Action, Nothing}; 
			plotargs...)
		
		xlims = (cart_pole_bounds.lower[1] - 0.3, cart_pole_bounds.upper[1] + 0.3)
		x = state[1]
		θ = π/2 - state[3]
		
		plot!(;
			xlims,
			ylims=(-m.pole_length, m.pole_length*3),
			yticks=nothing,
			aspectratio=:equal,
			plotargs...)
		
		# Action 
		if action == nothing
		elseif action.value == left.value
			plot!([x, x - 0.3], [0, 0],
				label=nothing,
				marker=:ltriangle,
				markerstrokewidth=0,
				color=colors.EMERALD,)
		elseif action.value == right.value
				plot!([x, x + 0.3], [0, 0],
				label=nothing,
				marker=:rtriangle,
				markerstrokewidth=0,
				color=colors.EMERALD,)
		else
			error("Unexpected action value $action")
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
	
	function plot_state(state::CartPoleState, action::Union{Action, Nothing};
			plotargs...)
		
		plot()
		plot_state!(state, action; plotargs...)
	end
end

# ╔═╡ bd006c20-2291-4e26-af86-2aaa4eeb0640
theme_type; plot_state(rand(trace.states), rand((left, right)), xlim=(-1, 1))

# ╔═╡ 251427d4-0aae-4a1f-a31f-71832877c749
# ╠═╡ disabled = true
#=╠═╡
function animate_sequence(trace::CartPoleTrace; speed=1)
	🎥 = @animate for (i, s) in enumerate(trace.states)

		a = get(trace.actions, i, nothing)
		plot_state(trace.states[i], a)
	end
	gif(🎥, fps=1/m.τ*speed, show_msg=false)
end
  ╠═╡ =#

# ╔═╡ ea84c513-b4ca-41df-96ec-1c230fde9f3d
#=╠═╡
animate_sequence(trace; speed=1)
  ╠═╡ =#

# ╔═╡ 227f0131-fc76-4382-902e-18874ce66104
cart_pole_bounds

# ╔═╡ dfdaf4cc-3490-4b93-a8d9-9f4d01c39c09
let
	all_good = true
	for i in 1:100
		s = [rand(cart_pole_bounds.lower[i]:0.0001:cart_pole_bounds.upper[i])
			for i in 1:4]
	
		s = CartPoleState(s..., 0)

		(s′,) = f⁻¹(f(s))
		if (concerned_with_position && (s[1] ≉ s′[1] || s[2] ≉ s′[2]) ||
			concerned_with_angle && (s[3] ≉ s′[3] || s[4] ≉ s′[4])) &&
			(concerned_with_position && (s[1] ≉ s″[1] || s[2] ≉ s″[2]) ||
				concerned_with_angle && (s[3] ≉ s″[3] || s[4] ≉ s″[4]))
			@error "s ≉ f⁻¹(f(s))" s f(s) f⁻¹(f(s))
			all_good = false
			break
		end
	end
	if all_good
		"The inverse function seems to work 👍"
	else
		"Found an example where the inverse function doesn't work :-("
	end
end

# ╔═╡ e0892cae-9ef0-4e57-9a1c-91bf34043956
md"""
## The Grid
"""

# ╔═╡ 38bc7025-9e1f-4101-a53d-a3a7ff802aa7
grid_bounds = let
	if enable_altered_state_space
		lower = Float64[-2.4, -5, -0.2095, -3]
		upper = Float64[ 2.4,  5,  0.2095,  3]
	else
		lower = Float64[-2.4, -10, -0.2095, -3]
		upper = Float64[ 2.4,  10,  0.2095,  3]
	end
	Bounds(lower, upper)
end

# ╔═╡ 453f70e4-c22f-4c75-b786-d523c8e4bf9c
begin
	function is_safe(s::AlteredState)
		for (i, l) in enumerate(grid_bounds.lower)
			if !concerned_with_position && (i == 1 || i == 2)
				 continue
			elseif !concerned_with_angle && (i == 3 || i == 4)
				continue
			end
			if l ≈ s[i]
				return false
			end
		end
		for (i, u) in enumerate(grid_bounds.upper)
			if !concerned_with_position && (i == 1 || i == 2)
				 continue
			elseif !concerned_with_angle && (i == 3 || i == 4)
				continue
			end
			if u ≈ s[i]
				return false
			end
		end
		return true
	end

	function is_safe(bounds::Bounds)
		for (i, l) in enumerate(bounds.lower)
			if !concerned_with_position && (i == 1 || i == 2)
				 continue
			elseif !concerned_with_angle && (i == 3 || i == 4)
				continue
			end
			if l ≈ grid_bounds.lower[i]
				return false
			end
		end
		for (i, u) in enumerate(bounds.upper)
			if !concerned_with_position && (i == 1 || i == 2)
				 continue
			elseif !concerned_with_angle && (i == 3 || i == 4)
				continue
			end
			if u ≈ grid_bounds.upper[i]
				return false
			end
		end
		return true
	end
end

# ╔═╡ 16315ce1-9dca-4284-936b-32a204b56108
is_safe(f(s0()))

# ╔═╡ 1aa6c36a-7884-48d1-9c5f-1dc87f1278e3
@bind resolution PlutoUI.combine() do field
	fields = [field("$i", NumberField(1:1000, default=20)) for i in 1:4]
	md"""
	### 🛠 `resolution` 
	
	$fields"""
end

# ╔═╡ 5dde6492-564f-46cb-848d-8a28ea2adb5f
granularity = let
	l, u = grid_bounds.lower, grid_bounds.upper
	span = u .- l

	if concerned_with_angle && concerned_with_position
		Float64[
			span[1]/resolution[1],
			span[2]/resolution[2],
			span[3]/resolution[3],
			span[4]/resolution[4]]
	elseif concerned_with_angle
		Float64[
			span[1],
			span[2],
			span[3]/resolution[3],
			span[4]/resolution[4]]
	elseif concerned_with_position
		Float64[
			span[1]/resolution[1],
			span[2]/resolution[2],
			span[3],
			span[4]]
	else # Not really valid
		Float64[
			span[1],
			span[2],
			span[3],
			span[4]]
	end
end

# ╔═╡ bcbf4a16-ce8f-451e-b58b-0bf9d8d0d872
get_size(granularity, cart_pole_bounds)

# ╔═╡ 50506cd4-bdf7-4efc-a12d-6ffc9a0f70a7
# ╠═╡ disabled = true
#=╠═╡
granularity = let
	l, u = grid_bounds.lower, grid_bounds.upper
	span = u .- l

	[span[1], span[2], 0.005, 0.006]
end
  ╠═╡ =#

# ╔═╡ 0610d08b-020e-4ec8-9815-1d0a4c592899
get_size(granularity, grid_bounds)

# ╔═╡ 6f4c7279-43f4-4c08-8eb1-45e3d70db8cd
# Add European-style thousands separator. Stolen from StackOverflow.
function commas(num::Integer)
	str = string(num)
	return replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

# ╔═╡ 299658d1-c3df-48a2-b992-02ef94c1bb59
get_size(granularity, grid_bounds) |> prod |> commas

# ╔═╡ 668f4592-75fd-445e-a0fa-56ee02a03f2d
no_action = actions_to_int([])

# ╔═╡ 611ba5df-6af9-413b-8e8a-b3da0c825d3e
any_action = actions_to_int([left.value, right.value])

# ╔═╡ 26cfc8ec-1351-468f-b9dc-e76acec6e777
function initializer(bounds::Bounds)
	is_safe(bounds) ? any_action : no_action
end

# ╔═╡ b966dc17-050f-40dc-adee-b6f9e79b4b0c
grid = let
	grid = Grid(granularity, grid_bounds)
	GridShielding.initialize!(grid, initializer)
	grid
end

# ╔═╡ f723aa48-e30b-4666-ad70-c20ae10fb4bb
is_safe(Bounds(box(grid, grid.bounds.lower)))

# ╔═╡ 030e5edc-0c1c-460f-b6b4-602640f5caf7
size(grid)

# ╔═╡ 7b1d742c-2b70-4aab-a032-cca6814d0d98
length(grid)

# ╔═╡ ace4eb92-b880-4fe6-9391-ad5bc586b802
md"""
## Reachability
"""

# ╔═╡ dee84d76-f3a2-45a6-b1db-3ca865877de1
const samples_per_axis = [2, 2, 2, 2]

# ╔═╡ 6b1c4273-a6f4-4101-b85e-cd63c308f8cc
crude_clamp!(x, bounds::Bounds) = begin
	for i in 1:get_dim(bounds)
		x[i] = clamp(x[i], bounds.lower[i], bounds.upper[i] - 0.0001)
	end
	x
end

# ╔═╡ 5376f447-716b-4011-bd3a-4b60db5ed110
function reachability_function(partition::Partition, action)::Vector{Vector{Int64}}
	result = Vector{Int64}[]
	grid = partition.grid
	for s::AlteredState in SupportingPoints(samples_per_axis, partition)
		#for r in SupportingPoints(samples_per_random_axis, Bounds((-1,), (1,)))
		for f⁻¹s in f⁻¹(s)
			s′ = simulate_point(m, f⁻¹s, get_action(action))
			s′ = f(s′)
			crude_clamp!(s′, grid.bounds)
			partition′ = box(grid, s′)
			if partition′.indices ∈ result
				continue
			end
			push!(result, partition′.indices)
		end
		#end
	end
	result
end

# ╔═╡ 5c6c1a0d-1442-4d39-bc4d-f7c943e14d97
md"""
## Mainmatter
"""

# ╔═╡ dac58385-3443-436d-acf4-dc15ce28c4af
begin reachability_function, grid, concerned_with_position, concerned_with_angle
	@bind do_it_button CounterButton("Do it.")
end

# ╔═╡ b5453b51-a878-4073-8152-c69d85d30ec1
if do_it_button > 0 || true
	reachability_function_precomputed = get_transitions(reachability_function, 
		[0, 1],
		grid);
end

# ╔═╡ a17770c5-de14-4c29-8efb-3d49f96a2950
simulate_point(m, [-2.4000000000000004, 2.0, -0.41800000000000004, -1.5, 0], left
) |> f

# ╔═╡ fa3ae0ca-d4b1-4961-afaa-c0d98174e0d2
size(grid)

# ╔═╡ b7980c62-7947-40e8-a8c5-195e93f4eb18
length(grid)

# ╔═╡ b79619f1-aa55-4a5c-851f-7387d411d8eb
@bind max_steps NumberField(0:1000, default=1000)

# ╔═╡ 0c28089a-1547-47f4-a411-e3a57cac6a6d
if @isdefined reachability_function_precomputed
	shield, max_steps_reached = 
		make_shield(reachability_function_precomputed, [0, 1], grid; max_steps)
else
	shield, max_steps_reached = grid, true
end

# ╔═╡ abaa6617-7932-4a10-a355-b2218bad4103
md"""

### 🛠 slice
`slice_axis_1 =` $(@bind slice_axis_1 Select(
	[i => n for (i, n) in enumerate(altered_state_axes)], 
	default=(concerned_with_position ? 1 : 3)))

`slice_axis_2 =` $(@bind slice_axis_2 Select(
	[i => n for (i, n) in enumerate(altered_state_axes)],
	default=(concerned_with_position ? 2 : 4)))
"""

# ╔═╡ 551c4f5c-615e-4f99-9549-eb9926bf6450
if enable_altered_state_space let
	theme_type # reactivity
	
	if slice_axis_2 < slice_axis_1
		sa1, sa2 = slice_axis_2, slice_axis_1
	else
		sa1, sa2 = slice_axis_1, slice_axis_2
	end
	xlabel=altered_state_axes[1]
	ylabel=altered_state_axes[2]
	l, u = shield.bounds.lower, shield.bounds.upper
	
	draw_function(
		s -> let
			s′ = f(CartPoleState(0, 0, s..., 4))
			if s′ ∉ shield
				return 0
			else
				return get_value(box(shield, s′))
			end
		end,
		-0.2, 0.2, -3, 3, (make_paper_friendly_figures ? 0.001 : 0.005);
		color=cgrad([colors.WET_ASPHALT, colors.AMETHYST, colors.SUNFLOWER, colors.CLOUDS], 10, categorical=true),
		xlabel="\$\\theta\$",
		ylabel="\$\\omega\$",
		colorbar=nothing)
end end

# ╔═╡ 6de525db-e339-435f-9f87-620fed817839
md"""
### 🛠 `s`

`x =`
$(@bind x NumberField(cart_pole_bounds.lower[1]:round(granularity[1], digits=3):cart_pole_bounds.upper[1], default=0))

`x_vel =`
$(@bind x_vel NumberField(cart_pole_bounds.lower[2]:round(granularity[2], digits=3):cart_pole_bounds.upper[2], default=0))

`θ =`
$(@bind θ NumberField(cart_pole_bounds.lower[3]:round(granularity[3], digits=3):cart_pole_bounds.upper[3], default=0))

`θ_vel =`
$(@bind θ_vel NumberField(cart_pole_bounds.lower[4]:round(granularity[4], digits=3):cart_pole_bounds.upper[4], default=0))

`action =`
$(@bind action Select([a => a.name for a in (left, right)]))
"""

# ╔═╡ a5a3b815-23c4-4acd-90d5-a12c721e7866
action

# ╔═╡ 1b882558-e83e-4679-8d51-3dc54040cdf1
s = CartPoleState(apply(action, [x, x_vel, θ, θ_vel, 0]))

# ╔═╡ 65f2972d-d1c1-4e21-95ed-f73f2047093d
@time simulate_point(m, s, action)

# ╔═╡ 23317253-20ab-44fe-b6f3-82d40307f5be
s

# ╔═╡ a2917851-62f0-471f-b57a-f64c14526f56
f(s)

# ╔═╡ 6127c10d-87c5-4123-87d8-7986a4f3a311
f⁻¹(f(s))

# ╔═╡ e83d4cda-ac7b-4efa-8200-b4a3f3fca38f
f⁻¹(f(s))

# ╔═╡ 19ef2f43-8338-4661-9375-b4836b81b647
f(s)

# ╔═╡ 1b797c47-31dd-45df-94f0-07bd376b57b7
s == f⁻¹(f(s))

# ╔═╡ 3ae14c03-d786-4e79-8744-3c52a8f4266d
f(s)

# ╔═╡ 915e0813-216c-4f01-8341-75c57198dc44
partition = box(grid, f(s))

# ╔═╡ 1d8867a6-a0b9-4ee1-8cf0-dabfa4678937
bounds = Bounds(partition)

# ╔═╡ f4c3e866-50aa-440b-a141-65de2daf9c4c
is_safe(bounds)

# ╔═╡ 641cc511-cb53-4d08-81f1-43a94b3fbb1c
is_safe(bounds)

# ╔═╡ f908b62b-4183-4ee8-9dbc-cab4a8164e70
begin
	slice = Any[partition.indices...]
	slice[slice_axis_1] = slice[slice_axis_2] = Colon()
	slice
end

# ╔═╡ b1de4c47-e90e-45c3-8c60-340516b42f8e
@time reachability_function(partition, action.value)

# ╔═╡ e77abd9a-de23-40b8-a442-b0971339f903
reachability_function_precomputed[action.value][partition.indices...]

# ╔═╡ 891b2c11-c79b-4fc2-8188-cf7a1097bb6d
@bind show_reachability CheckBox(default=false)

# ╔═╡ 1161cbd5-7e47-4358-9b7e-139dcd6740a1
@bind zoom CheckBox(default=false)

# ╔═╡ 4bae6425-1730-4d6a-8d79-a7534f9d131a
bounds

# ╔═╡ 481c90d7-dbb1-4ddb-aebd-66d018c27d92
s, f(s)

# ╔═╡ a05e58c6-9bf1-4865-9052-a1a4a231f3b2
show_grid = size(grid)[slice_axis_1] < 50 && size(grid)[slice_axis_2] < 50

# ╔═╡ c7a4e65c-a907-468e-b31c-ce05393d41d5
p3 = let
	theme_type # reactivity
	
	if slice_axis_2 < slice_axis_1
		sa1, sa2 = slice_axis_2, slice_axis_1
	else
		sa1, sa2 = slice_axis_1, slice_axis_2
	end
	
	
	xlabel=altered_state_axes[sa1]
	ylabel=altered_state_axes[sa2]
	
	draw(shield, slice;
		show_grid,
		colors=shieldcolors,
		#color_labels=shieldlabels,
		clims=(0, 3),
		xlabel,
		ylabel,
		legend=:top)

	if zoom
		plot!(
		xlim=(s[sa1] - 1, s[sa1] + 1),
		ylim=(s[sa2] - 1, s[sa2] + 1),)
	end
	if show_reachability
		sp = SupportingPoints(samples_per_axis, bounds)
		xs = [s[sa1] for s in sp]
		ys = [s[sa2] for s in sp]
		
		scatter!(xs, ys,
			color=colors.EMERALD,
			markerstrokewidth=0,
			label="initial")

		sp = [simulate_point(m, s′, action) 
			for s in sp
			for s′ in f⁻¹(AlteredState(s))
			]
				
		xs = [s[sa1] for s in sp]
		ys = [s[sa2] for s in sp]
		
		scatter!(xs, ys,
			color=colors.BELIZE_HOLE,
			markerstrokewidth=0,
			markersize=2,
			label="reachable")

		#=
		reachable = reachability_function_precomputed[action.value][partition.indices...]
		reachable = [Partition(shield, is) for is in reachable]
		reachable = [Bounds(p) for p in reachable]
		
		reachable = [
			Bounds([b.lower[sa1], b.lower[sa2]], [b.upper[sa1], b.upper[sa2]]) 
			for b in reachable]

		reachable = reachable |> unique
		
		for r in reachable
			plot!(r,
				label=nothing,
				linewidth=1,
				color=colors.PETER_RIVER)
		end
		=#
	end
	plot!()
end

# ╔═╡ 0871379f-8cf8-4949-a033-d64c9e3e633d
if max_steps_reached
	md"""
	!!! warning "Shield not done"
		Either synthesis has not even been started, or the `max_steps` variable controlling the number of iterations has been set too low."""
end

# ╔═╡ cc196ed4-c434-4244-9284-353b48686a33
round_8(n) = round(n, digits=8)

# ╔═╡ de02fd2c-cfd0-429a-9602-0e3da5418c70
prettyprint(b::Bounds) = 
		"Bounds($(round_8.(b.lower)), $(round_8.(b.upper)))"

# ╔═╡ 8f424501-1641-4779-bf60-0204f6ea3efc
shieldcolors[get_value(partition) + 1], prettyprint(bounds)

# ╔═╡ 5b373e7a-6254-4fe4-bead-eababbd8f065
# Reachability from s
let
	
	reachable = reachability_function_precomputed[action.value][partition.indices...]
	reachable = [Partition(shield, indices) for indices in reachable]
	reachable = [(get_value(partition), partition) for partition in reachable]
	reachable = [(shieldcolors[v+1], Bounds(p)) for (v, p) in reachable]
	reachable = [(c, prettyprint(b)) for (c, b) in reachable]
end

# ╔═╡ 973fba84-d206-454c-a743-0d9eae296c28
md"""
## Fit a Polynomial
"""

# ╔═╡ acbdbfe1-66e1-45d9-81d6-96a059aafb6f
@bind polynomial_degree NumberField(1:24, default=3)

# ╔═╡ a54f79bb-a240-47e5-9f04-8b2674ac9be1
@bind fit_to Select(["Averaged", "Upper", "Lower"])

# ╔═╡ 05f7c96a-dd50-41a7-8916-293938c03b40
# ╠═╡ disabled = true
#=╠═╡
let
	if slice_axis_2 < slice_axis_1
		sa1, sa2 = slice_axis_2, slice_axis_1
	else
		sa1, sa2 = slice_axis_1, slice_axis_2
	end
	
	plot(p3)
	lower_border = border_points(shield, 2, 3, slice) |> sort
	upper_border = border_points(shield, 1, 3, slice) |> sort
	
	averaged_border = [(x1, (l + u)/2) 
			for ((x1, l), (x2, u)) in zip(upper_border, lower_border)
			if x1 == x2]

	pol1 = Polynomials.fit((averaged_border[1:15] |> unzip)..., 1)
	pol2 = Polynomials.fit((averaged_border[16:85] |> unzip)..., 0)
	pol3 = Polynomials.fit((averaged_border[86:end] |> unzip)..., 1)
	@show pol1
	@show pol2
	@show averaged_border[15]
	@show pol3
	@show averaged_border[85]
	
	piecewise(x) =  
			x < averaged_border[15][1] ? pol1(x) : 
			x < averaged_border[85][1] ? pol2(x) : pol3(x)

	f(x, x_vel) = piecewise(x) + x_vel

	plot!([(x, f(x, x_vel)) for x in grid_bounds.lower[sa1]:0.1:grid_bounds.upper[sa1]],
		color=colors.SILVER,
		label="p")
	
	scatter!(averaged_border[1:end],
		color=colors.ASBESTOS,
		markersize=2,
		markerstrokewidth=0,
		label="average")

	scatter!([x], [x_vel], 
		color=colors.EMERALD,
		markersize=2,
		markerstrokewidth=0,
		label="s")

end
  ╠═╡ =#

# ╔═╡ 442a2427-e8a7-4088-8221-ec7a5dc9f1c2
function border_points(grid::Grid, value_1, value_2, slice)
	result = Tuple{Float64, Float64}[]
	sa1, sa2 = indexof((==)(Colon()), slice)
	for partition in grid
		if get_value(partition) != value_1
			continue
		elseif !all([partition.indices[i] == index 
				for (i, index) in enumerate(slice) 
				if index != Colon()])
			continue
		end
		above = copy(partition.indices)
		above[sa2] += 1
		above = Partition(grid, above)
		below = copy(partition.indices)
		below[sa2] -= 1
		below = Partition(grid, below)
		if get_value(above) == value_2 || get_value(below) == value_2
			bounds = Bounds(partition)
			# Middle of the bounds
			x = (bounds.upper[sa1] - bounds.lower[sa1])/2 + bounds.lower[sa1]
			y = (bounds.upper[sa2] - bounds.lower[sa2])/2 + bounds.lower[sa2]
			push!(result, (x, y))
		end
	end
	result
end

# ╔═╡ d73c2ee5-e8bf-4cc4-8855-d239224ba843
p4 = let
	theme_type # reactivity
	
	if slice_axis_2 < slice_axis_1
		sa1, sa2 = slice_axis_2, slice_axis_1
	else
		sa1, sa2 = slice_axis_1, slice_axis_2
	end
	
	plot(p3)
	lower_border = border_points(shield, 1, 3, slice)
	upper_border = border_points(shield, 2, 3, slice)
	
	averaged_border = [(x1, (l + u)/2) 
			for (x1, l) in lower_border 
			for (x2, u) in upper_border
			if x1 == x2]

	to_fit_to = fit_to == "Upper" ? upper_border :
				fit_to == "Lower" ? lower_border :
				fit_to == "Averaged" ? averaged_border : error("whu? wha?")

	p = Polynomials.fit((to_fit_to |> unzip)..., polynomial_degree)
	@show p
	
	scatter!(upper_border,
		legend=:topright,
		color=colors.CONCRETE,
		markersize=3,
		marker=:utriangle,
		markerstrokewidth=0,
		label="Upper")
	
	scatter!(averaged_border,
		color=colors.CONCRETE,
		markersize=3,
		marker=:diamond,
		markerstrokewidth=0,
		label="Average")
	
	scatter!(lower_border,
		color=colors.CONCRETE,
		markersize=3,
		marker=:dtriangle,
		markerstrokewidth=0,
		label="Lower")

	plot!([(x, p(x)) 
			for x in grid_bounds.lower[sa1]:granularity[sa1]/2:grid_bounds.upper[sa1]],
		color=colors.SILVER,
		xlim=(shield.bounds.lower[sa1], shield.bounds.upper[sa1]),
		ylim=(shield.bounds.lower[sa2], shield.bounds.upper[sa2]),
		label="p")
end

# ╔═╡ 2498792a-a7b9-4295-bfb9-7e9068a02d7d
random_policy

# ╔═╡ 0f5f2fea-3b84-4d1d-80bc-08715f947661
md"""
## Check Safety
"""

# ╔═╡ 9fda178a-0fcd-49ad-b0b8-025523995691
# ╠═╡ disabled = true
#=╠═╡
animate_sequence(shielded_trace)
  ╠═╡ =#

# ╔═╡ e2db38df-8347-4bf4-be27-d9ad19c96823
function get_allowed(s::CartPoleState)
	f_s = f(s)
	f_s = crude_clamp!(f_s, shield.bounds)
	allowed = int_to_actions(Int, get_value((box(shield, f_s))))
	allowed = [get_action(a) for a in allowed]
end

# ╔═╡ d08f05d8-6227-4bbf-aab2-744152726107
function shield_policy(policy)
	return s -> begin
		a = policy(s)
		#if f(s) ∉ shield
		#	error("Outside grid: s=$s f(s)=$(f(s))")
		#end
		allowed = get_allowed(s)
		if a ∈ allowed
			return a
		elseif length(allowed) > 0
			return rand(allowed)
		else
			error("Unsafe state reached: $s")
		end
	end
end

# ╔═╡ fbaa39a7-d1d2-4ad8-a475-c712cdafe35d
shielded_random_policy = shield_policy(random_policy)

# ╔═╡ d6ee45cf-765a-41d6-8bc7-b74662ac9243
s0_const = s0()

# ╔═╡ d3f51f26-92da-4e23-9316-13a249079100
get_allowed(s0_const)

# ╔═╡ 2dbb749a-cd7e-4092-b662-519b10d9552d
runs = 100

# ╔═╡ 08d65223-892b-4921-9d09-af959524bb7a
function check_safety(m::CartPoleMechanics, policy, duration; runs=1000)
	deaths = 0
	example_trace = nothing
	@progress for run in 1:runs
		trace = simulate_sequence(m, s0(), policy, duration)
		for s in trace.states
			if (concerned_with_position &&
					!(cart_pole_bounds.lower[1] < s[1] < cart_pole_bounds.upper[1])
				) ||
				(concerned_with_angle &&
					!(cart_pole_bounds.lower[3] < s[3] < cart_pole_bounds.upper[3])
				)
				
				deaths += 1
				example_trace = trace
				break
			end
		end
		example_trace = something(example_trace, trace)
	end
	deaths, example_trace
end

# ╔═╡ 84734786-a79c-484a-95a9-5de041436c2f
cart_pole_bounds

# ╔═╡ 9bad8bb4-bfa1-47a3-821c-dd3448c3f534
deaths, shielded_trace = check_safety(m, shielded_random_policy, 10; runs)

# ╔═╡ 3d63e8c3-6218-4eec-86de-698bb61d8f96
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

# ╔═╡ 65f542c6-d5f3-40e0-be5f-ab66786eaf72
shielded_trace.states[end]

# ╔═╡ 834ab28c-a08a-4d92-8ab5-5284196cc2db
md"""
# Exporting the Shield
"""

# ╔═╡ 15952e74-05fb-4ffe-b6b6-97e1bd5fc815
@bind target_dir TextField(95, default=mktempdir())

# ╔═╡ 8f78d419-7c9d-4d92-9d66-a9696a203b07
target_dir; @bind open_folder_button CounterButton("Open Folder")

# ╔═╡ d7915d1f-7bd3-463e-a5bb-d31b6fe9fe23
if open_folder_button > 0
	run(`nautilus $target_dir`, wait=false)
end; "This cell opens `$target_dir` in nautilus"

# ╔═╡ 4e27eb77-69b4-4477-9f21-2e8a9702898c
md"""
### Export as serialized julia-tuple

Easy export and import between julia code.
"""

# ╔═╡ 22df2512-8b9f-47c0-aaaf-ee32f3675d3d
let
	if enable_altered_state_space
		state_space = "Altered State Space"
	else
		state_space = "Standard State Space"
	end
	filename = "Cart Pole - $state_space.shield"
	robust_grid_serialization(joinpath(target_dir, filename), shield)
	"Exported `'$filename'`." |> Markdown.parse
end

# ╔═╡ 6052d4f5-6874-4bc1-b17f-d323a07d6c54
md"""
### Export as a function in a shared-object library

Use this library to access the shield from C and C++ code.

The shield is compiled into a shared-object binary, which exports the function `int get_value(double v, double p)`. It takes the state-variables as input and returns the bit-encoded list of allowed actions. (See `int_to_actions`.)
"""

# ╔═╡ 747c5cee-f701-4c28-b7df-04efa3740d61
let
	if enable_altered_state_space
		state_space = "altered_state_space"
	else
		state_space = "standard_state_space"
	end
	shield_so = "cart_pole_shield_$state_space.so"
	shield_so = joinpath(target_dir, shield_so)
	
	get_libshield(shield; destination=shield_so, force=true)
	
	"Exported `'$shield_so'`." |> Markdown.parse
end

# ╔═╡ d2217f57-024c-4338-bcdb-361990f94d37
md"""
### Export to Numpy

Exports a zip-file containing a serialized numpy-array along with a JSON file with details on how to read it.
"""

# ╔═╡ 2864b91e-7654-4c54-9b7a-411ef982d01d
@enum Actions move_left move_right

# ╔═╡ 0a9c51e1-0cc3-4bb0-9fa5-0f3962aae605
@bind episodes NumberField(0:typemax(Int64), default=5)

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

# ╔═╡ 40550db6-a7b9-496c-b321-b61cf1239e18
function Q_learn()
	Q = Dict(left => Grid(Q_granularity, Q_bounds, data_type=Float64), 
			right => Grid(Q_granularity, Q_bounds, data_type=Float64))

	# Discourage exploration; we want to stay near s0
	for partition in Q[left]
		set_value!(partition, -2)
	end
	for partition in Q[right]
		set_value!(partition, -2)
	end
	
	@progress for i ∈ 1:episodes
		Sₜ = s0()
		Aₜ = rand((left, right))
		for t ∈ 0:m.τ:10
			Sₜ₊₁ = simulate_point(m, Sₜ, Aₜ)
			if Sₜ₊₁ ∉ Q_bounds continue end
			Q_Sₜ_Aₜ = box(Q[Aₜ], Sₜ)
			set_value!(Q_Sₜ_Aₜ, 
				get_value(Q_Sₜ_Aₜ) + 
				α(t)*(r(Sₜ) + γ*max([get_value(box(Q[a′], Sₜ₊₁)) 
				for a′ in (left, right)]...) - get_value(Q_Sₜ_Aₜ)))
			
			Aₜ₊₁ = ϵ_greedy(ϵ(t), Q, Sₜ)
			Sₜ, Aₜ = Sₜ₊₁, Aₜ₊₁
			if Sₜ₊₁ ∉ cart_pole_bounds break end
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

# ╔═╡ 45785f69-c79a-4172-b8b4-9009ee08e613
Q_trace.states[end] ∈ cart_pole_bounds

# ╔═╡ eb53b3db-db53-4842-ad48-4272a667b7cf
max([abs(x_vel) for (x, x_vel, θ, θ_vel, _) in Q_trace.states]...)

# ╔═╡ 5bfe3632-dba7-4e12-ba21-823b9803b9df
max([abs(θ_vel) for (x, x_vel, θ, θ_vel, _) in Q_trace.states]...)

# ╔═╡ 2ac120e4-f380-4b02-bc7f-a1d5e84d7c36
#=╠═╡
animate_sequence(Q_trace)
  ╠═╡ =#

# ╔═╡ d712571e-ced8-4f06-8b44-6874fcd3e15d
length(Q[left].array |> unique),
length(Q[right].array |> unique)

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

# ╔═╡ 618ab712-95b6-413b-973f-8309405b7da9
let
	if enable_altered_state_space
		state_space = "Altered State Space"
	else
		state_space = "Standard State Space"
	end

	meta_info = (;variables=state_axes_simple, 
		actions=Actions,
		env_id="Cart Pole")
	
	filename = "Cart Pole Shield - $state_space.zip"
	
	numpy_zip_file(shield, joinpath(target_dir, filename); meta_info...)
	
	"Exported `'$filename'`." |> Markdown.parse, meta_info
end

# ╔═╡ 3c578986-4a73-42bd-a202-fbf2bf534151
get_allowed(CartPoleState(0, 0, -0.16, -0.99, 0))

# ╔═╡ Cell order:
# ╟─2767663f-3ef8-44f5-81a2-8e480158266e
# ╟─3115801b-0a07-4a44-a6b3-d1ab2b9c0775
# ╠═cb6e988a-f263-11ee-1f3f-53192cebcad4
# ╠═cd2df9dc-af72-4b37-b1ef-ff8a0dcb9e0f
# ╠═5c808e32-0aa6-49ec-a7b1-d01e44aa48ea
# ╟─a8aff15c-255d-498f-97dd-4c9c953ec662
# ╟─319eff14-0f2e-458e-99b3-b3a64eccccee
# ╠═1d5aced4-a27a-4f3b-ab08-69b1a3b1559f
# ╠═35c76906-f2dd-4f4d-af43-0fafe69d211b
# ╠═8172cdb6-b9a4-4afc-844f-245e5d951bb7
# ╠═e008fe83-3851-4061-9e5a-35e938a53ec1
# ╠═a8eefbb3-a09d-4915-9b9a-76f896636dcc
# ╟─3fd479d1-c43a-4c6f-95f8-0c74a9ffbf18
# ╟─7dd5c185-2b95-4297-a36c-4e2ca38952ea
# ╠═35605d87-4c3a-49a9-93d1-a5fceede3653
# ╠═b52604cc-e8bc-4b53-84ad-79cf019c1667
# ╠═e1f217af-759a-4868-b0b2-6ce08de324ea
# ╠═274e6f36-bcae-4232-9bb6-732f954aa4e5
# ╠═788e089e-ded6-4ad7-9951-1d10a32a8295
# ╠═86ab3957-129c-4c71-911d-a54fd80cfd2d
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
# ╠═58458e09-c9cf-4d75-b153-17d807b783f3
# ╠═bd006c20-2291-4e26-af86-2aaa4eeb0640
# ╠═251427d4-0aae-4a1f-a31f-71832877c749
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
# ╠═45785f69-c79a-4172-b8b4-9009ee08e613
# ╠═227f0131-fc76-4382-902e-18874ce66104
# ╠═eb53b3db-db53-4842-ad48-4272a667b7cf
# ╠═5bfe3632-dba7-4e12-ba21-823b9803b9df
# ╠═2ac120e4-f380-4b02-bc7f-a1d5e84d7c36
# ╟─6c84e1aa-f45a-453c-8a78-f3976c605385
# ╟─a7296808-4eb8-4f10-8683-adc4963b21ce
# ╟─9da1c57f-c922-420c-8c46-0d6842553f8c
# ╠═338e5d40-6251-429c-9b0d-ef92460a7e52
# ╠═3aeb0922-e5f8-4311-b66d-dd42d61f18f3
# ╠═b4ddbac3-33e9-42e7-bac7-a2ec32093678
# ╠═64be4f5a-97f6-49bc-b848-450506d9ceb1
# ╠═5645d7a7-0f23-4a02-ab30-49a6cda9ce17
# ╠═d90ea6c4-316d-46b1-8688-23d95f3cca61
# ╠═6eedcd3c-1aa5-4d8b-9e5f-8465e608f6c7
# ╠═1eb460d5-3de2-44b0-9a75-378a35e8cf6f
# ╠═f52cb596-8f54-4721-879c-65d8f206a224
# ╠═23317253-20ab-44fe-b6f3-82d40307f5be
# ╠═a2917851-62f0-471f-b57a-f64c14526f56
# ╠═6127c10d-87c5-4123-87d8-7986a4f3a311
# ╠═e83d4cda-ac7b-4efa-8200-b4a3f3fca38f
# ╟─edfae222-4bd5-4411-bfa4-a6e9fc5ff725
# ╠═19ef2f43-8338-4661-9375-b4836b81b647
# ╠═1b797c47-31dd-45df-94f0-07bd376b57b7
# ╠═dfdaf4cc-3490-4b93-a8d9-9f4d01c39c09
# ╟─731d3595-1746-46cd-9851-e285673d6a1b
# ╟─2f2a6934-b633-47f1-9a26-88b7ece30f46
# ╠═e5371840-2989-4f7f-b15c-4b07e8e96a3b
# ╠═453f70e4-c22f-4c75-b786-d523c8e4bf9c
# ╠═16315ce1-9dca-4284-936b-32a204b56108
# ╠═f4c3e866-50aa-440b-a141-65de2daf9c4c
# ╠═f723aa48-e30b-4666-ad70-c20ae10fb4bb
# ╟─e0892cae-9ef0-4e57-9a1c-91bf34043956
# ╠═38bc7025-9e1f-4101-a53d-a3a7ff802aa7
# ╟─1aa6c36a-7884-48d1-9c5f-1dc87f1278e3
# ╟─5dde6492-564f-46cb-848d-8a28ea2adb5f
# ╠═50506cd4-bdf7-4efc-a12d-6ffc9a0f70a7
# ╠═0610d08b-020e-4ec8-9815-1d0a4c592899
# ╠═6f4c7279-43f4-4c08-8eb1-45e3d70db8cd
# ╠═299658d1-c3df-48a2-b992-02ef94c1bb59
# ╠═668f4592-75fd-445e-a0fa-56ee02a03f2d
# ╠═611ba5df-6af9-413b-8e8a-b3da0c825d3e
# ╠═26cfc8ec-1351-468f-b9dc-e76acec6e777
# ╠═b966dc17-050f-40dc-adee-b6f9e79b4b0c
# ╠═030e5edc-0c1c-460f-b6b4-602640f5caf7
# ╠═7b1d742c-2b70-4aab-a032-cca6814d0d98
# ╠═1b882558-e83e-4679-8d51-3dc54040cdf1
# ╠═3ae14c03-d786-4e79-8744-3c52a8f4266d
# ╠═915e0813-216c-4f01-8341-75c57198dc44
# ╠═1d8867a6-a0b9-4ee1-8cf0-dabfa4678937
# ╠═641cc511-cb53-4d08-81f1-43a94b3fbb1c
# ╟─ace4eb92-b880-4fe6-9391-ad5bc586b802
# ╠═dee84d76-f3a2-45a6-b1db-3ca865877de1
# ╠═6b1c4273-a6f4-4101-b85e-cd63c308f8cc
# ╠═5376f447-716b-4011-bd3a-4b60db5ed110
# ╠═b1de4c47-e90e-45c3-8c60-340516b42f8e
# ╟─5c6c1a0d-1442-4d39-bc4d-f7c943e14d97
# ╠═dac58385-3443-436d-acf4-dc15ce28c4af
# ╠═b5453b51-a878-4073-8152-c69d85d30ec1
# ╠═a17770c5-de14-4c29-8efb-3d49f96a2950
# ╠═e77abd9a-de23-40b8-a442-b0971339f903
# ╠═fa3ae0ca-d4b1-4961-afaa-c0d98174e0d2
# ╠═b7980c62-7947-40e8-a8c5-195e93f4eb18
# ╠═0c28089a-1547-47f4-a411-e3a57cac6a6d
# ╠═b79619f1-aa55-4a5c-851f-7387d411d8eb
# ╟─c7a4e65c-a907-468e-b31c-ce05393d41d5
# ╟─551c4f5c-615e-4f99-9549-eb9926bf6450
# ╟─abaa6617-7932-4a10-a355-b2218bad4103
# ╟─f908b62b-4183-4ee8-9dbc-cab4a8164e70
# ╟─6de525db-e339-435f-9f87-620fed817839
# ╠═891b2c11-c79b-4fc2-8188-cf7a1097bb6d
# ╠═1161cbd5-7e47-4358-9b7e-139dcd6740a1
# ╠═4bae6425-1730-4d6a-8d79-a7534f9d131a
# ╠═481c90d7-dbb1-4ddb-aebd-66d018c27d92
# ╠═a05e58c6-9bf1-4865-9052-a1a4a231f3b2
# ╟─0871379f-8cf8-4949-a033-d64c9e3e633d
# ╠═cc196ed4-c434-4244-9284-353b48686a33
# ╠═de02fd2c-cfd0-429a-9602-0e3da5418c70
# ╠═8f424501-1641-4779-bf60-0204f6ea3efc
# ╠═5b373e7a-6254-4fe4-bead-eababbd8f065
# ╟─973fba84-d206-454c-a743-0d9eae296c28
# ╠═acbdbfe1-66e1-45d9-81d6-96a059aafb6f
# ╠═a54f79bb-a240-47e5-9f04-8b2674ac9be1
# ╠═d73c2ee5-e8bf-4cc4-8855-d239224ba843
# ╠═05f7c96a-dd50-41a7-8916-293938c03b40
# ╠═442a2427-e8a7-4088-8221-ec7a5dc9f1c2
# ╠═2498792a-a7b9-4295-bfb9-7e9068a02d7d
# ╟─0f5f2fea-3b84-4d1d-80bc-08715f947661
# ╠═9fda178a-0fcd-49ad-b0b8-025523995691
# ╠═d08f05d8-6227-4bbf-aab2-744152726107
# ╠═e2db38df-8347-4bf4-be27-d9ad19c96823
# ╠═fbaa39a7-d1d2-4ad8-a475-c712cdafe35d
# ╠═d6ee45cf-765a-41d6-8bc7-b74662ac9243
# ╠═d3f51f26-92da-4e23-9316-13a249079100
# ╠═2dbb749a-cd7e-4092-b662-519b10d9552d
# ╠═08d65223-892b-4921-9d09-af959524bb7a
# ╠═84734786-a79c-484a-95a9-5de041436c2f
# ╠═9bad8bb4-bfa1-47a3-821c-dd3448c3f534
# ╟─3d63e8c3-6218-4eec-86de-698bb61d8f96
# ╠═65f542c6-d5f3-40e0-be5f-ab66786eaf72
# ╟─834ab28c-a08a-4d92-8ab5-5284196cc2db
# ╠═15952e74-05fb-4ffe-b6b6-97e1bd5fc815
# ╟─8f78d419-7c9d-4d92-9d66-a9696a203b07
# ╟─d7915d1f-7bd3-463e-a5bb-d31b6fe9fe23
# ╟─4e27eb77-69b4-4477-9f21-2e8a9702898c
# ╠═22df2512-8b9f-47c0-aaaf-ee32f3675d3d
# ╟─6052d4f5-6874-4bc1-b17f-d323a07d6c54
# ╠═747c5cee-f701-4c28-b7df-04efa3740d61
# ╟─d2217f57-024c-4338-bcdb-361990f94d37
# ╠═2864b91e-7654-4c54-9b7a-411ef982d01d
# ╠═618ab712-95b6-413b-973f-8309405b7da9
# ╠═3c578986-4a73-42bd-a202-fbf2bf534151
