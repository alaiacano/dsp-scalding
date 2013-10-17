x = [1 2 3]

# julia> toeplitz(x)
# 3x3 Array{Int64,2}:
#  1  2  3
#  2  1  2
#  3  2  1
function toeplitz{T}(x::Array{T})
    n = length(x)
    A = zeros(T, n, n)
    
    for i in 1:n
        A[i,i:] = x[1:n - i + 1]
        A[i:,i] = x[1:n - i + 1]
    end
    A
end

# modified version. It still creates a toeplitz matrix, but is useful for convolution.
# julia> toeplitz_filter(x)
# 3x3 Array{Int64,2}:
#  1  2  3
#  0  1  2
#  0  0  1
function toeplitz_filter{T}(x::Array{T})
    n = length(x)
    A = zeros(T, n, n)
    for i in 1:n
        A[i,i:] = x[1:n - i + 1]
    end
    A
end

# to add some padding for more proper convolution:
# julia> toeplitz_filter([1 2 3 0 0 0 0 0 0 0 0 0 0])
# 13x13 Array{Int64,2}:
#  1  2  3  0  0  0  0  0  0  0  0  0  0
#  0  1  2  3  0  0  0  0  0  0  0  0  0
#  0  0  1  2  3  0  0  0  0  0  0  0  0
#  0  0  0  1  2  3  0  0  0  0  0  0  0
#  0  0  0  0  1  2  3  0  0  0  0  0  0
#  0  0  0  0  0  1  2  3  0  0  0  0  0
#  0  0  0  0  0  0  1  2  3  0  0  0  0
#  0  0  0  0  0  0  0  1  2  3  0  0  0
#  0  0  0  0  0  0  0  0  1  2  3  0  0
#  0  0  0  0  0  0  0  0  0  1  2  3  0
#  0  0  0  0  0  0  0  0  0  0  1  2  3
#  0  0  0  0  0  0  0  0  0  0  0  1  2
#  0  0  0  0  0  0  0  0  0  0  0  0  1


# use in a linear filter:
t = linspace(0, 6pi, 1000)

# input: square wave and its components
X = [sin(2*pi*t) sin(3*2*pi*t)/3 sin(5*2*pi*t)/5]
X = [X sum(X,2)]

# build the filter
H = toeplitz_filter([repeat([1/7], outer=[1,7]) repeat([0],outer=[1,length(t)-7])])

# filtered output
Y = transpose(transpose(X) * H)
