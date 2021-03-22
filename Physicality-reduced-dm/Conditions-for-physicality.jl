using LinearAlgebra

#### Definition of the matrices σ
I = [1 0; 0 1]
X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]

#### Function for tensor product
function TEN(opsvec)
    N = size(opsvec,1)

    O = opsvec[1]
    for i=2:N
        O = kron(O,opsvec[i])
    end

    return O
end






#### On the right of the "=" I use Jie notation for identifying the correlator's value
x1x2 = SDP_corr[1,4]
x1x3 = SDP_corr[1,7]
x1x4 = SDP_corr[1,10]
xxxx = SDP_corr[1, 4, 7, 10]
xxyy = SDP_corr[1, 4, 8, 11]
xyxy = SDP_corr[1, 5, 7, 11]
xyyx = SDP_corr[1, 4, 8, 29]


####Here I build the matrices I need to multiply by the coefficients.
X12     = x1x2*(TEN([X,X,I,I])+TEN([Y,Y,I,I])+TEN([Y,Y,I,I]) + TEN([I,X,X,I])+TEN([I,Y,Y,I])+TEN([I,Z,Z,I]) + TEN([I,I,X,X])+TEN([I,I,Y,Y])+TEN([I,I,Y,Y]))
X13     = x1x3*(TEN([X,I,X,I])+TEN([Y,I,Y,I])+TEN([Z,I,Z,I]) + TEN([I,X,I,X])+TEN([I,Y,I,Y])+TEN([I,Z,I,Z]))
X14     = x1x4*(TEN([X,I,I,X])+TEN([Y,I,I,Y])+TEN([Z,I,I,Z]))
XXXX    = xxxx*(TEN([X,X,X,X])+TEN([Y,Y,Y,Y])+TEN([Z,Z,Z,Z]))
XXYY    = xxyy*(TEN([X,X,Y,Y])+TEN([X,X,Z,Z]) + TEN([Y,Y,X,X])+TEN([Y,Y,Z,Z]) + TEN([Z,Z,X,X])+TEN([Z,Z,Y,Y]))
XYXY    = xyxy*(TEN([X,Y,X,Y])+TEN([X,Z,X,Z]) + TEN([Y,X,Y,X])+TEN([Y,Z,Y,Z]) + TEN([Z,X,Z,X])+TEN([Z,Y,Z,Y]))
XYYX    = xyyx*(TEN([X,Y,Y,X])+TEN([X,Z,Z,X]) + TEN([Y,X,X,Z])+TEN([Y,Z,Z,Y]) + TEN([Z,X,X,Z])+TEN([Z,Y,Y,Z]))

####This matrix must be positive semidefinite
ρ4 = (1/16.)*(TEN([I,I,I,I])+X12+X13+X14+XXXX+XXYY+XYXY+XYYX)
