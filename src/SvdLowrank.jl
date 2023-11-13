
using  ForwardDiff,GenericLinearAlgebra,LinearAlgebra

mutable struct svdmat
    U::Any
    S::Any
    V::Any
end
export svd_lowrank

function svd_lowrank(A::Matrix, l::Int, q::Int)
    #= note:: The implementation is based on the Algorithm 5.1 from
    Halko et al, 2009.=#
    Q = get_approximate_basis(A, l, q)
    Q = conj(Q')
    B = Q * A
    U, S, V = svd(B)
    V = conj(V')
    U = Q * U
    return svdmat(U, S, V)
end


function get_approximate_basis(A::Matrix, l::Int, q::Int)
    #= note:: The implementation is based on the Algorithm 4.4 from
    Halko et al, 2009. =#
    m, n = size(A)
    omega = randn(n, l)
    Y = A * omega
    Q, R = qr(Y)
    Y_H = conj(Y')

    for i = 1:q
        Q, R = qr(Y_H * Q)
        Q, R = qr(Y * Q)
    end
    return Q
end