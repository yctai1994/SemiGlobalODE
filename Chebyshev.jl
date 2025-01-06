### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 79d44ed1-3968-4372-a96c-38dc8f87440a
begin
	import Pkg
	Pkg.activate()

	import CairoMakie as cm
end

# ╔═╡ 91814e80-eaf7-4eab-92d7-e32f1ab1fbaf
test_func(t::Real) = 0.1 * abs2(sech((t − 500) / 140)) * cos(0.06 * (t − 500))

# ╔═╡ 192519a8-c953-4049-be72-2e90ae695e1d
global_test_xarr = collect(range(start = 0.0, step = 1e-2, stop = 1e3))

# ╔═╡ 71cc131e-8f37-4104-9123-e15ddb2856ab
global_test_yarr = let yarr = similar(global_test_xarr)
	@inbounds for i in eachindex(yarr)
		yarr[i] = test_func(global_test_xarr[i])
	end
	yarr
end

# ╔═╡ d1f77238-dc23-4441-84bf-aec4b7586c34
md"---"

# ╔═╡ f05169ec-cc4f-11ef-2aed-3d58cadc3837
function num_recipes_chebyshev_coeff(
	func::Function, a::Real, b::Real, n:: Int
)
	bpa = 0.5 * (b + a) # right_bound + left_bound
	bma = 0.5 * (b - a) # right_bound - left_bound

	x = Vector{Float64}(undef, n)
	f = Vector{Float64}(undef, n)
	c = Vector{Float64}(undef, n)

	# for k in 0:n-1
	# 	x[k+1] = cos(pi * (k + 0.5) / n) * bma + bpa
	# 	f[k+1] = func(x[k+1])
	# end

	for k in 1:n
		x[k] = cos(pi * (2k - 1) / (2n)) * bma + bpa
		f[k] = func(x[k])
	end

	fac = 2.0 / n
	for j in 0:n-1
		tot = 0.0
		for k in 0:n-1
			tot += f[k+1] * cos(pi * j * (k + 0.5) / n)
		end
		c[j+1] = fac * tot
	end

	return c
end

# ╔═╡ d5c055d8-fa44-4f87-9d04-30e789f343ff
num_recipes_c = num_recipes_chebyshev_coeff(
	test_func, 0.0, global_test_xarr[end], 101
)

# ╔═╡ f82eaceb-821f-41df-8d1e-6f32cc3c09ec
function num_recipes_chebyshev_eval(
	x::Real, a::Real, b::Real, coeff::Vector{Float64}
)
	if (x - a) * (x - b) > 0.0
		error("x not in range in Chebyshev::eval")
	end

	bpa = b + a
	bma = b - a

	y = (2.0 * x - bpa) / bma
	y2 = 2.0 * y

	j = length(coeff)
	dⱼ₋₁ = 0.0
	dⱼ = 0.0
	while j > 1
		sv = dⱼ₋₁
		dⱼ₋₁ = y2 * dⱼ₋₁ - dⱼ + coeff[j]
		dⱼ = sv
		j -= 1
	end

	return y * dⱼ₋₁ - dⱼ + 0.5 * coeff[1]
end

# ╔═╡ 02af5007-404e-4bef-891b-a7ec0c8a35ef
md"---"

# ╔═╡ c72c4e0a-ccc2-4bbf-ab91-3b1f6e80aea6
num_recipes_x = collect(range(start = 0.0, step = 1e1, stop = 1e3))

# ╔═╡ 84dbe091-49b6-496d-ab24-1ddfd932eaf3
num_recipes_y = let yarr = similar(num_recipes_x)
	@inbounds for i in eachindex(yarr)
		yarr[i] = num_recipes_chebyshev_eval(
			num_recipes_x[i],
			global_test_xarr[1],
			global_test_xarr[end],
			num_recipes_c
		)
	end
	yarr
end

# ╔═╡ c61fc524-01af-4f49-b622-bff9211dc2e8
md"---"

# ╔═╡ 37efcab2-4827-48c2-9311-f2c1efed8297
let fig = cm.Figure()
	ax1 = cm.Axis(fig[1,1];
		limits = (global_test_xarr[1], global_test_xarr[end], -0.105, 0.105,),
		xticks = collect(range(
			start = global_test_xarr[1],
			step = 100,
			stop = global_test_xarr[end],
		)),
	)
	cm.lines!(ax1,
		global_test_xarr, global_test_yarr;
		color = (:steelblue, 0.5),
	)
	cm.scatterlines!(ax1,
		num_recipes_x, num_recipes_y;
		color = (:firebrick, 0.5),
		markercolor = (:firebrick, 0.5),
	)
	fig
end

# ╔═╡ Cell order:
# ╠═79d44ed1-3968-4372-a96c-38dc8f87440a
# ╟─91814e80-eaf7-4eab-92d7-e32f1ab1fbaf
# ╟─192519a8-c953-4049-be72-2e90ae695e1d
# ╟─71cc131e-8f37-4104-9123-e15ddb2856ab
# ╟─d1f77238-dc23-4441-84bf-aec4b7586c34
# ╟─f05169ec-cc4f-11ef-2aed-3d58cadc3837
# ╠═d5c055d8-fa44-4f87-9d04-30e789f343ff
# ╠═f82eaceb-821f-41df-8d1e-6f32cc3c09ec
# ╟─02af5007-404e-4bef-891b-a7ec0c8a35ef
# ╟─c72c4e0a-ccc2-4bbf-ab91-3b1f6e80aea6
# ╠═84dbe091-49b6-496d-ab24-1ddfd932eaf3
# ╟─c61fc524-01af-4f49-b622-bff9211dc2e8
# ╠═37efcab2-4827-48c2-9311-f2c1efed8297
