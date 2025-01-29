### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a5881ab9-053c-47ac-b23c-da051ecc171d
import LinearAlgebra as la

# ╔═╡ 109c3c7c-ceb2-11ef-3483-89b230552d1f
const A0 = ComplexF64[
	complex(1.0, 2.0) complex(-5.0, 1.0) complex(2.0, 1.0) complex(-3.0, -1.0);
	complex(3.0, 1.0) complex(-2.0, 3.0) complex(6.0, 7.0) complex(-3.0, -4.0);
	complex(0.0, 0.0) complex(-1.0, 5.0) complex(1.0, 2.0) complex(-1.0, -3.0);
	complex(0.0, 0.0) complex( 0.0, 0.0) complex(7.0, 6.0) complex(-4.0, -3.0);
]

# ╔═╡ 8e2a472e-881f-428b-b7fe-c128b3fee75e
function qr!(A::Matrix{T}, R::Matrix{T}, Qt::Matrix{T}) where T<:Complex
	N = size(A, 2)

	@inbounds for j in 1:N
		for i in 1:N
			Qt[i,j] = complex(0.0, 0.0)
		end
		@inbounds Qt[j,j] = complex(1.0, 0.0)
	end

	copyto!(R, A)

	for k in 1:N-1
		scale = maximum(abs, @view(R[k:N,k]))
		@inbounds for i in k:N
			R[i,k] /= scale
		end

		sigma = la.norm(@view(R[k:N,k]))
		phase = cis(angle(R[k,k]))
		beta = inv(sigma * (sigma + abs(R[k,k])))
		R[k,k] += sigma * phase

		for j in k+1:N
			tmp = la.dot(@view(R[k:N,k]), @view(R[k:N,j]))
			tmp *= beta
			for i in k:N
				R[i,j] -= tmp * R[i,k]
			end
		end

		for j in 1:N
			tmp = la.dot(@view(R[k:N,k]), @view(Qt[k:N,j]))
			tmp *= beta
			for i in k:N
				Qt[i,j] -= tmp * R[i,k]
			end
		end

		@inbounds R[k,k] = -scale * sigma * phase
		@inbounds for i in k+1:N
			R[i,k] = complex(0.0, 0.0)
		end
	end
end

# ╔═╡ 5e2b2de8-2530-480f-b911-ea8ce569bb0f
function rq!(A::Matrix{T}, R::Matrix{T}, Qt::Matrix{T}) where T<:Complex
	N = size(A, 2)

	@inbounds for i in 1:N, j in 1:N
		tmp = complex(0.0, 0.0)
		for k in i:N
			tmp += R[i,k] * conj(Qt[j,k])
		end
		A[i,j] = tmp
	end
end

# ╔═╡ ee08c69d-e40c-4596-92b9-5a5cd3a294b5
function wilkinson(A::AbstractMatrix{T}) where T<:Complex
	m = 0.5 * (A[1,1] + A[2,2])
	p = A[1,1] * A[2,2] - A[1,2] * A[2,1]

	z2 = m * m - p
	r = abs(z2)
	s = sqrt(r)

	θ = angle(z2)

	z_plus = s * cis(0.5θ)
	z_minus = s * cis(0.5θ + π)

	eigen_minus = m + z_minus
	eigen_plus  = m + z_plus

	if abs(eigen_minus - A[2,2]) < abs(eigen_plus - A[2,2])
		return eigen_minus
	else
		return eigen_plus
	end
end

# ╔═╡ 39ce686e-e072-44b8-b35a-5c9218a2d976
la.eigen(copy(A0))

# ╔═╡ b1d7f708-8ebe-4372-a117-f4ecf21d3e88
let _A_ = copy(A0), _R_ = similar(_A_), _Qt_ = similar(_A_)
	N = size(_A_, 2)
	for _ in 1:16
		eigen_guess = wilkinson(@view(_A_[N-1:N, N-1:N]))
		@inbounds for i in 1:N
			_A_[i,i] -= eigen_guess
		end

		qr!(_A_, _R_, _Qt_)
		rq!(_A_, _R_, _Qt_)

		@inbounds for i in 1:N
			_A_[i,i] += eigen_guess
		end
	end
	print(_A_)
end

# ╔═╡ df96a77d-7dfd-4468-8690-c56bfc0fdf70
md"---"

# ╔═╡ 9bd88489-110a-49fd-891e-2bad132d4b70
const B0 = ComplexF64[
	complex(1.0, 2.0) complex(-5.0, 1.0);
	complex(6.0, 1.0) complex(-6.0, 7.0);
]

# ╔═╡ bae57ded-9b4f-4d3d-a735-3729d8d24cf5
la.eigen(copy(B0))

# ╔═╡ 709d278d-18e6-4255-9d7c-e5142af42054
let (m, p) = (
	0.5 * (B0[1,1] + B0[2,2]),
	B0[1,1] * B0[2,2] - B0[1,2] * B0[2,1],
)
	z2 = m * m - p
	r = abs(z2)
	s = sqrt(r)

	θ = angle(z2)

	z_plus = s * cis(0.5θ)
	z_minus = s * cis(0.5θ + π)

	m + z_minus, m + z_plus
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═a5881ab9-053c-47ac-b23c-da051ecc171d
# ╟─109c3c7c-ceb2-11ef-3483-89b230552d1f
# ╟─8e2a472e-881f-428b-b7fe-c128b3fee75e
# ╟─5e2b2de8-2530-480f-b911-ea8ce569bb0f
# ╟─ee08c69d-e40c-4596-92b9-5a5cd3a294b5
# ╠═39ce686e-e072-44b8-b35a-5c9218a2d976
# ╠═b1d7f708-8ebe-4372-a117-f4ecf21d3e88
# ╟─df96a77d-7dfd-4468-8690-c56bfc0fdf70
# ╠═9bd88489-110a-49fd-891e-2bad132d4b70
# ╠═bae57ded-9b4f-4d3d-a735-3729d8d24cf5
# ╠═709d278d-18e6-4255-9d7c-e5142af42054
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
