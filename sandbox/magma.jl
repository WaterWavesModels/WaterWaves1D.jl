using MAGMA
import LinearAlgebra.LAPACK: gels!, gesv!, getrs!, getri!
using Test
using CuArrays
using TimerOutputs

T = ComplexF64

# randomly generate a 2 by 2 matrix for testing
size_test = 1024

A = rand(T, size_test, size_test)
B = rand(T, size_test, size_test)

dA = cu(A)
dB = cu(B)

A_copy = copy(A)
B_copy = copy(B)

A_test = dA
B_test = dB 

right_answer = gesv!(A, B)

magma_init()
result = magma_gesv!(A_test, B_test)
magma_finalize()

for i in 1:length(result)
    @test Array(result[i]) ≈ right_answer[i]
end
