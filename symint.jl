using Plots

struct PhaseVec
    x::Array{Float32, 1}
    p::Array{Float32, 1}
    m::Float32
    c::Float32
end

function ΣF(q, W, F)
    Σ = q.x .* 0.0
    for i in 1:length(W)
        Σ = Σ + F(q, W[i])
    end
    return Σ
end

function F_gc(u, v, G, k_e)
    r_sq = sum((u.x - v.x).^2)
    dx = u.x - v.x
    r_sq == 0 && return 0.0*dx
    dx = dx/sqrt(sum(dx.^2))
    f = (G*u.m*v.m + k_e*u.c*v.c)
    return dx.*f/sum((u.x - v.x).^2)
end

function cdstep(u::PhaseVec, F, t, c, d)
    # x = x + c*t*v
    # p = p + d*t*F
    return PhaseVec(u.x + c*t*u.p/u.m, u.p + d*t*F, u.m, u.c)
end

function verletstep(W, F, t)
    #verlet with c1 = 0, c2 = 1
    #and d1 = d2 = 1/2
    for i in 1:length(W)
        u = W[i]
        Fu = ΣFV(u, W, F)
        #c = 0, d1 = 1/2
        W[i] = cdstep(u, Fu, t, 0, 0.5)
    end
    for i in 1:length(W)
        u = W[i]
        Fu = ΣFV(u, W, F)
        #c = 1, d = 1/2
        W[i] = cdstep(u, Fu, t, 1.0, 0.5)
    end
    return W
end

function V_gc(u, v, G, k_e)
    r_sq = sum((u.x - v.x).^2)
    r_sq == 0 && return 0.0
    f = (G*u.m*v.m + k_e*u.c*v.c)
    return f/sqrt(sum((u.x - v.x).^2))
end


function potential(W, V)
    Σ = 0.0
    for i in 1:length(W)
        u = W[i]
        for j in 1:i
            v = W[j]
            Σ = Σ + V(u, v)
        end
    end
    return Σ
end

function verletstep(W, F, t, c)
    #verlet with c1 = 0, c2 = 1
    #and d1 = d2 = 1/2
    for i in 1:length(W)
        u = W[i]
        Fu = ΣF(u, W, F)
        Fu += u.x .* c*(1-xsq(u))
        #c = 0, d1 = 1/2
        W[i] = cdstep(u, Fu, t, 0, 0.5)
    end
    for i in 1:length(W)
        u = W[i]
        Fu = ΣF(u, W, F)
        Fu += u.x .* c*(1-xsq(u))
        #c = 1, d = 1/2
        W[i] = cdstep(u, Fu, t, 1.0, 0.5)
    end
    return W
end


function vsq(u)
    return sum((u.p./u.m).^2)
end

function xsq(u)
    return sum(u.x.^2)
end

function drag(W, λ)
    return [PhaseVec(u.x,u.p - c*vsq(u).*u.p, u.m, u.c) for u in W]
end

function decay(W, λ)
    return [PhaseVec(u.x, λ .*u.p, u.m, u.c) for u in W]
end

function x(u::PhaseVec)
    return u.x
end

function xyz(u::PhaseVec)
    return [u.x[1], u.x[2], u.x[3]]
end

function constrain(u::PhaseVec)
    x = u.x
    p = u.p
    x = x ./ sqrt(x'*x)
    p = p - x.*(x'*p)/(x'*x)
    return PhaseVec(x, p, u.m, u.c)
end



function main(n)
    G_0 = -1.0
    k_e0 = 1.0
    F(u, v) = F_gc(u, v, G_0, k_e0)
    W = []
    z = PhaseVec(0.0*rand(3), 0.0* rand(3), 10000.0, 0.0)
    push!(W,z)
    for i in 1:n
        push!(W, PhaseVec(2rand(3).-0.5, 2rand(3) .- 0.5, 1, 1))
    end
    pts = [xyz(a) for a in W]
    for i in 1:1000
        W = verletstep(W, F, 0.01)
        #W = decay(W, 0.99)
    end
    return W
end



function dynam(n)
    G_0 = -0.001
    k_e0 = 1.0
    F(u, v) = F_gc(u, v, G_0, k_e0)
    V(u, v) = V_gc(u, v, G_0, k_e0)
    W = []
    for i in 1:n
        push!(W, PhaseVec(2rand(3).-0.5, 2rand(3) .- 0.5, 1, 1))
    end
    W = [constrain(u) for u in W]
    for i in 1:100000
        W = verletstep(W, F, 0.01, 0.1)
        W = decay(W, 0.999)
        if i % 10 == 0
            println(potential(W, V))
        end
    end
    return W
end

v(θ) = [[cos(θ) -sin(θ) 0]; [sin(θ) cos(θ) 0]; [0 0 1]]

function dynmplot(n)
    W = dynam(n)
    pts = [xyz(u) for u in W]
    anim = @animate for k in 1:314
        scatter3d([a[1] for a in pts], [a[2] for a in pts], [a[3] for a in pts])
        pts = [v(0.05) *pt for pt in pts]
    end
    gif(anim, "/tmp/anim_llfps15.gif", fps = 15)

end
