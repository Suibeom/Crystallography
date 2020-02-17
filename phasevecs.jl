struct PhaseVec{T}
    x::T
    y::T
end

struct Hamiltonian{T, N}
    H(PhaseVec{T})::N
end
